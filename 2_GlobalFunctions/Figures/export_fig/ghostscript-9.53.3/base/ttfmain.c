/* Copyright (C) 2001-2020 Artifex Software, Inc.
   All Rights Reserved.

   This software is provided AS-IS with no warranty, either express or
   implied.

   This software is distributed under license and may not be copied,
   modified or distributed except as expressly authorized under the terms
   of the license contained in the file LICENSE in this distribution.

   Refer to licensing information at http://www.artifex.com or contact
   Artifex Software, Inc.,  1305 Grant Avenue - Suite 200, Novato,
   CA 94945, U.S.A., +1(415)492-9861, for further information.
*/


/* A Free Type interface adapter. */
/* Uses code fragments from the FreeType project. */

#include "ttmisc.h"
#include "ttfoutl.h"
#include "ttfmemd.h"

#include "ttfinp.h"
#include "ttfsfnt.h"
#include "ttobjs.h"
#include "ttinterp.h"
#include "ttcalc.h"

static const bool skip_instructions = 0; /* Debug purpose only. */

typedef struct {
    TT_Fixed a, b, c, d, tx, ty;
} FixMatrix;

struct ttfSubGlyphUsage_s {
    FixMatrix m;
    int index;
    int flags;
    short arg1, arg2;
};

/*------------------------------------------------------------------- */

static TT_Fixed AVE(F26Dot6 a, F26Dot6 b)
{   return (a + b) / 2;
}

static F26Dot6 shortToF26Dot6(short a)
{   return (F26Dot6)a << 6;
}

static F26Dot6 floatToF26Dot6(float a)
{   return (F26Dot6)(a * (1 << 6) + 0.5);
}

static TT_Fixed floatToF16Dot16(float a)
{   return (F26Dot6)(a * (1 << 16) + 0.5);
}

static void TransformF26Dot6PointFix(F26Dot6Point *pt, F26Dot6 dx, F26Dot6 dy, FixMatrix *m)
{   pt->x = MulDiv(dx, m->a, 65536) + MulDiv(dy, m->c, 65536) + (m->tx >> 10);
    pt->y = MulDiv(dx, m->b, 65536) + MulDiv(dy, m->d, 65536) + (m->ty >> 10);
}

static void TransformF26Dot6PointFloat(FloatPoint *pt, F26Dot6 dx, F26Dot6 dy, FloatMatrix *m)
{   pt->x = dx * m->a / 64 + dy * m->c / 64 + m->tx;
    pt->y = dx * m->b / 64 + dy * m->d / 64 + m->ty;
}

/*-------------------------------------------------------------------*/

static ttfPtrElem *ttfFont__get_table_ptr(ttfFont *f, char *id)
{
    if (!memcmp(id, "cvt ", 4))
        return &f->t_cvt_;
    if (!memcmp(id, "fpgm", 4))
        return &f->t_fpgm;
    if (!memcmp(id, "glyf", 4))
        return &f->t_glyf;
    if (!memcmp(id, "head", 4))
        return &f->t_head;
    if (!memcmp(id, "hhea", 4))
        return &f->t_hhea;
    if (!memcmp(id, "hmtx", 4))
        return &f->t_hmtx;
    if (!memcmp(id, "vhea", 4))
        return &f->t_vhea;
    if (!memcmp(id, "vmtx", 4))
        return &f->t_vmtx;
    if (!memcmp(id, "loca", 4))
        return &f->t_loca;
    if (!memcmp(id, "maxp", 4))
        return &f->t_maxp;
    if (!memcmp(id, "prep", 4))
        return &f->t_prep;
    if (!memcmp(id, "cmap", 4))
        return &f->t_cmap;
    return 0;
}

/*-------------------------------------------------------------------*/

TT_Error  TT_Set_Instance_CharSizes(TT_Instance  instance,
                                       TT_F26Dot6   charWidth,
                                       TT_F26Dot6   charHeight)
{
    PInstance  ins = instance.z;

    if ( !ins )
        return TT_Err_Invalid_Instance_Handle;

    if (charWidth < 1*64)
        charWidth = 1*64;

    if (charHeight < 1*64)
        charHeight = 1*64;

    ins->metrics.x_scale1 = charWidth;
    ins->metrics.y_scale1 = charHeight;
    ins->metrics.x_scale2 = ins->face->font->nUnitsPerEm;
    ins->metrics.y_scale2 = ins->face->font->nUnitsPerEm;

    if (ins->face->font->nFlags & 8) {
        ins->metrics.x_scale1 = (ins->metrics.x_scale1+32) & -64;
        ins->metrics.y_scale1 = (ins->metrics.y_scale1+32) & -64;
    }

    ins->metrics.x_ppem = ins->metrics.x_scale1 / 64;
    ins->metrics.y_ppem = ins->metrics.y_scale1 / 64;

    if (charWidth > charHeight)
        ins->metrics.pointSize = charWidth;
    else
        ins->metrics.pointSize = charHeight;

    ins->valid  = FALSE;
    return Instance_Reset(ins, FALSE);
  }

/*-------------------------------------------------------------------*/

int ttfInterpreter__obtain(ttfMemory *mem, ttfInterpreter **ptti)
{
    ttfInterpreter *tti;

    if (*ptti) {
        (*ptti)->lock++;
        return 0;
    }
    tti = mem->alloc_struct(mem, (const ttfMemoryDescriptor *)&st_ttfInterpreter, "ttfInterpreter__obtain");
    if (!tti)
        return fMemoryError;
    tti->usage = 0;
    tti->usage_size = 0;
    tti->ttf_memory = mem;
    tti->lock = 1;
    tti->exec = mem->alloc_struct(mem, (const ttfMemoryDescriptor *)&st_TExecution_Context, "ttfInterpreter__obtain");
    if (!tti->exec) {
        mem->free(mem, tti, "ttfInterpreter__obtain");
        return fMemoryError;
    }
    memset(tti->exec, 0, sizeof(*tti->exec));
    *ptti = tti;
    return 0;
}

void ttfInterpreter__release(ttfInterpreter **ptti)
{
    ttfInterpreter *tti = *ptti;
    ttfMemory *mem = tti->ttf_memory;

    if(--tti->lock)
        return;
    mem->free(mem, tti->usage, "ttfInterpreter__release");
    mem->free(mem, tti->exec, "ttfInterpreter__release");
    mem->free(mem, *ptti, "ttfInterpreter__release");
    *ptti = 0;
}

/*-------------------------------------------------------------------*/

void ttfFont__init(ttfFont *self, ttfMemory *mem,
                    void (*DebugRepaint)(ttfFont *),
                    int (*DebugPrint)(ttfFont *, const char *s, ...),
                    const gs_memory_t *DebugMem)
{
    memset(self, 0, sizeof(*self));
    self->DebugRepaint = DebugRepaint;
    self->DebugPrint   = DebugPrint;
    self->DebugMem     = DebugMem;
}

void ttfFont__finit(ttfFont *self)
{   ttfMemory *mem = self->tti->ttf_memory;

    if (self->exec) {
        if (self->inst)
            Context_Destroy(self->exec);
        else {
            /* Context_Create was not called - see ttfFont__Open.
               Must not call Context_Destroy for proper 'lock' count.
             */
        }
    }
    self->exec = NULL;
    if (self->inst)
        Instance_Destroy(self->inst);
    mem->free(mem, self->inst, "ttfFont__finit");
    self->inst = NULL;
    if (self->face)
        Face_Destroy(self->face);
    mem->free(mem, self->face, "ttfFont__finit");
    self->face = NULL;
}

#define MAX_SUBGLYPH_NESTING 3 /* Arbitrary. We need this because we don't want
                                  a ttfOutliner__BuildGlyphOutline recursion
                                  while a glyph is loaded in ttfReader. */

FontError ttfFont__Open(ttfInterpreter *tti, ttfFont *self, ttfReader *r,
                                    unsigned int nTTC, float w, float h,
                                    bool design_grid)
{   char sVersion[4], sVersion1[4] = {0, 1, 0, 0};
    char sVersion2[4] = {0, 2, 0, 0};
    unsigned int nNumTables, i;
    TT_Error code, code1 = 0;
    int k;
    TT_Instance I;
    ttfMemory *mem = tti->ttf_memory;
    F26Dot6 ww, hh;

    self->tti = tti;
    self->design_grid = design_grid;
    r->Read(r, sVersion, 4);
    if(!memcmp(sVersion, "ttcf", 4)) {
        unsigned int nFonts;
        unsigned int nPos = 0xbaadf00d; /* Quiet compiler. */

        r->Read(r, sVersion, 4);
       if(memcmp(sVersion, sVersion1, 4) && memcmp(sVersion, sVersion2, 4))
            return fUnimplemented;
        nFonts = ttfReader__UInt(r);
        if (nFonts == 0)
            return fBadFontData;
        if(nTTC >= nFonts)
            return fTableNotFound;
        for(i = 0; i <= nTTC; i++)
            nPos = ttfReader__UInt(r);
        r->Seek(r, nPos);
        r->Read(r, sVersion, 4);
    }
    if(memcmp(sVersion, sVersion1, 4) && memcmp(sVersion, "true", 4))
        return fUnimplemented;
    nNumTables    = ttfReader__UShort(r);
    ttfReader__UShort(r); /* nSearchRange */
    ttfReader__UShort(r); /* nEntrySelector */
    ttfReader__UShort(r); /* nRangeShift */
    for (i = 0; i < nNumTables; i++) {
        char sTag[5];
        unsigned int nOffset, nLength;
        ttfPtrElem *e;

        sTag[4] = 0;
        r->Read(r, sTag, 4);
        ttfReader__UInt(r); /* nCheckSum */
        nOffset = ttfReader__UInt(r);
        nLength = ttfReader__UInt(r);
        e = ttfFont__get_table_ptr(self, sTag);
        if (e != NULL) {
            e->nPos = nOffset;
            e->nLen = nLength;
        }
    }
    r->Seek(r, self->t_head.nPos + offset_of(sfnt_FontHeader, flags));
    self->nFlags = ttfReader__UShort(r);
    r->Seek(r, self->t_head.nPos + offset_of(sfnt_FontHeader, unitsPerEm));
    self->nUnitsPerEm = ttfReader__UShort(r);
    if (self->nUnitsPerEm <= 0)
        self->nUnitsPerEm = 1024;
    r->Seek(r, self->t_head.nPos + offset_of(sfnt_FontHeader, indexToLocFormat));
    self->nIndexToLocFormat = ttfReader__UShort(r);
    r->Seek(r, self->t_maxp.nPos + offset_of(sfnt_maxProfileTable, numGlyphs));
    self->nNumGlyphs = ttfReader__UShort(r);
    r->Seek(r, self->t_maxp.nPos + offset_of(sfnt_maxProfileTable, maxComponentElements));
    self->nMaxComponents = ttfReader__UShort(r);
    if(self->nMaxComponents < 10)
        self->nMaxComponents = 10; /* work around DynaLab bug in lgoth.ttf */
    r->Seek(r, self->t_hhea.nPos + offset_of(sfnt_MetricsHeader, numberLongMetrics));
    self->nLongMetricsHorz = ttfReader__UShort(r);
    if (self->t_vhea.nPos != 0) {
        r->Seek(r, self->t_vhea.nPos + offset_of(sfnt_MetricsHeader, numberLongMetrics));
        self->nLongMetricsVert = ttfReader__UShort(r);
    } else
        self->nLongMetricsVert = 0;
    if (tti->usage_size < self->nMaxComponents * MAX_SUBGLYPH_NESTING) {
        tti->ttf_memory->free(tti->ttf_memory, tti->usage, "ttfFont__Open");
        tti->usage_size = 0;
        tti->usage = mem->alloc_bytes(mem,
                sizeof(ttfSubGlyphUsage) * self->nMaxComponents * MAX_SUBGLYPH_NESTING,
                "ttfFont__Open");
        if (tti->usage == NULL)
            return fMemoryError;
        tti->usage_size = self->nMaxComponents * MAX_SUBGLYPH_NESTING;
    }
    self->face = mem->alloc_struct(mem, (const ttfMemoryDescriptor *)&st_TFace, "ttfFont__Open");
    if (self->face == NULL)
        return fMemoryError;
    memset(self->face, 0, sizeof(*self->face));
    self->face->r = r;
    self->face->font = self;
    self->exec = tti->exec;
    code = Face_Create(self->face);
    if (code)
        return fMemoryError;
    code = r->Error(r);
    if (code < 0)
        return fBadFontData;
    self->inst = mem->alloc_struct(mem, (const ttfMemoryDescriptor *)&st_TInstance, "ttfFont__Open");
    if (self->inst == NULL)
        return fMemoryError;
    memset(self->inst, 0, sizeof(*self->inst));
    code = Context_Create(self->exec, self->face); /* See comment in the implementation of Context_Create. */
    if (code == TT_Err_Out_Of_Memory)
        return fMemoryError;
    code = Instance_Create(self->inst, self->face);
    if (code == TT_Err_Out_Of_Memory)
        return fMemoryError;
    if (code)
        return fBadFontData;
    for(k = 0; k < self->face->cvtSize; k++)
        self->inst->cvt[k] = shortToF26Dot6(self->face->cvt[k]);
    code = Instance_Init(self->inst);
    if (code == TT_Err_Out_Of_Memory)
        return fMemoryError;
    if (code >= TT_Err_Invalid_Opcode && code <= TT_Err_Invalid_Displacement)
        code1 = fBadInstruction;
    else if (code)
        return fBadFontData;
    I.z = self->inst;
    if (design_grid)
        ww = hh = shortToF26Dot6(self->nUnitsPerEm);
    else {
        /* Round towards zero for a better view of mirrored characters : */
        ww = floatToF26Dot6(w);
        hh = floatToF26Dot6(h);
    }
    code = TT_Set_Instance_CharSizes(I, ww, hh);
    self->inst->metrics  = self->exec->metrics;
    if (code == TT_Err_Invalid_Engine)
        return fPatented;
    if (code == TT_Err_Out_Of_Memory)
        return fMemoryError;
    if (code >= TT_Err_Invalid_Opcode && code <= TT_Err_Invalid_Displacement)
        return fBadInstruction;
    if (code)
        return fBadFontData;
    if (code1)
        return code1;
    return code;
}

static void ttfFont__StartGlyph(ttfFont *self)
{   Context_Load( self->exec, self->inst );
    if ( self->inst->GS.instruct_control & 2 )
        self->exec->GS = Default_GraphicsState;
    else
        self->exec->GS = self->inst->GS;
    self->tti->usage_top = 0;
}

static void ttfFont__StopGlyph(ttfFont *self)
{
    Context_Save(self->exec, self->inst);
}

/*-------------------------------------------------------------------*/

static void  mount_zone( PGlyph_Zone  source,
                          PGlyph_Zone  target )
{
    Int  np, nc;

    np = source->n_points;
    nc = source->n_contours;

    target->org_x = source->org_x + np;
    target->org_y = source->org_y + np;
    target->cur_x = source->cur_x + np;
    target->cur_y = source->cur_y + np;
    target->touch = source->touch + np;

    target->contours = source->contours + nc;

    target->n_points   = 0;
    target->n_contours = 0;
}

static void  Init_Glyph_Component( PSubglyph_Record    element,
                                   PSubglyph_Record    original,
                                   PExecution_Context  exec )
{
    element->index     = -1;
    element->is_scaled = FALSE;
    element->is_hinted = FALSE;

    if (original)
        mount_zone( &original->zone, &element->zone );
    else
        element->zone = exec->pts;

    element->zone.n_contours = 0;
    element->zone.n_points   = 0;

    element->arg1 = 0;
    element->arg2 = 0;

    element->element_flag = 0;
    element->preserve_pps = FALSE;

    element->transform.xx = 1 << 16;
    element->transform.xy = 0;
    element->transform.yx = 0;
    element->transform.yy = 1 << 16;

    element->transform.ox = 0;
    element->transform.oy = 0;

    element->leftBearing  = 0;
    element->advanceWidth = 0;
  }

static void  cur_to_org( Int  n, PGlyph_Zone  zone )
{
    Int  k;

    for ( k = 0; k < n; k++ )
        zone->org_x[k] = zone->cur_x[k];

    for ( k = 0; k < n; k++ )
        zone->org_y[k] = zone->cur_y[k];
}

static void  org_to_cur( Int  n, PGlyph_Zone  zone )
{
    Int  k;

    for ( k = 0; k < n; k++ )
        zone->cur_x[k] = zone->org_x[k];

    for ( k = 0; k < n; k++ )
        zone->cur_y[k] = zone->org_y[k];
}

/*-------------------------------------------------------------------*/

void ttfOutliner__init(ttfOutliner *self, ttfFont *f, ttfReader *r, ttfExport *exp,
                        bool bOutline, bool bFirst, bool bVertical)
{
    self->r = r;
    self->bOutline = bOutline;
    self->bFirst = bFirst;
    self->pFont = f;
    self->bVertical = bVertical;
    self->exp = exp;
}

static void MoveGlyphOutline(TGlyph_Zone *pts, int nOffset, ttfGlyphOutline *out, FixMatrix *m)
{   F26Dot6* x = pts->org_x + nOffset;
    F26Dot6* y = pts->org_y + nOffset;
    short count = out->pointCount;
    F26Dot6Point p;

    if (m->a == 65536 && m->b == 0 &&
        m->c == 0 && m->d == 65536 &&
        m->tx == 0 && m->ty == 0)
        return;
    for (; count != 0; --count) {
        TransformF26Dot6PointFix(&p, *x, *y, m);
        *x++ = p.x;
        *y++ = p.y;
    }
}

static FontError ttfOutliner__BuildGlyphOutlineAux(ttfOutliner *self, int glyphIndex,
            FixMatrix *m_orig, ttfGlyphOutline* gOutline)
{   ttfFont *pFont = self->pFont;
    ttfReader *r = self->r;
    ttfInterpreter *tti = pFont->tti;
    short sideBearing;
    FontError error = fNoError;
    short arg1, arg2;
    short count;
    unsigned int i;
    unsigned short nAdvance;
    unsigned int nPosBeg;
    TExecution_Context *exec = pFont->exec;
    TGlyph_Zone *pts = &exec->pts;
    TSubglyph_Record  subglyph;
    ttfSubGlyphUsage *usage = tti->usage + tti->usage_top;
    const byte *glyph = NULL;
    int glyph_size;
    bool execute_bytecode = true;
    int nPoints = 0;

retry:
    if (r->get_metrics(r, glyphIndex, self->bVertical, &sideBearing, &nAdvance) < 0) {
        /* fixme: the error code is missing due to interface restrictions. */
        goto errex;
    }
    gOutline->sideBearing = shortToF26Dot6(sideBearing);
    gOutline->advance.x = shortToF26Dot6(nAdvance);
    gOutline->advance.y = 0;
    self->bFirst = FALSE;

    if (!self->bOutline)
        return fNoError;
    if (!r->LoadGlyph(r, glyphIndex, &glyph, &glyph_size))
        return fGlyphNotFound;
    if (r->Eof(r)) {
        r->ReleaseGlyph(r, glyphIndex);
        gOutline->xMinB = gOutline->yMinB = 0;
        gOutline->xMaxB = gOutline->yMaxB = 0;
        return fNoError;
    }
    if (r->Error(r))
        goto errex;
    nPosBeg = r->Tell(r);

    gOutline->contourCount = ttfReader__Short(r);
    subglyph.bbox.xMin = ttfReader__Short(r);
    subglyph.bbox.yMin = ttfReader__Short(r);
    subglyph.bbox.xMax = ttfReader__Short(r);
    subglyph.bbox.yMax = ttfReader__Short(r);

    gOutline->xMinB = Scale_X(&exec->metrics, subglyph.bbox.xMin);
    gOutline->yMinB = Scale_Y(&exec->metrics, subglyph.bbox.yMin);
    gOutline->xMaxB = Scale_X(&exec->metrics, subglyph.bbox.xMax);
    gOutline->yMaxB = Scale_Y(&exec->metrics, subglyph.bbox.yMax);

    /* FreeType stuff beg */
    Init_Glyph_Component(&subglyph, NULL, pFont->exec);
    subglyph.leftBearing = sideBearing;
    subglyph.advanceWidth = nAdvance;
    subglyph.pp1.x = subglyph.bbox.xMin - sideBearing;
    subglyph.pp1.y = 0;
    subglyph.pp2.x = subglyph.pp1.x + nAdvance;
    subglyph.pp2.y = 0;
    /* FreeType stuff end */

    if (gOutline->contourCount == 0)
        gOutline->pointCount = 0;
    else if (gOutline->contourCount == -1) {
        unsigned short flags, index, bHaveInstructions = 0;
        unsigned int nUsage = 0;
        unsigned int nPos;
        unsigned int n_ins;

        gOutline->bCompound = TRUE;
        if (tti->usage_top + pFont->nMaxComponents > tti->usage_size)
            return fBadFontData;
        gOutline->contourCount = gOutline->pointCount = 0;
        do {
            FixMatrix m;
            ttfSubGlyphUsage *e;

            if (nUsage >= pFont->nMaxComponents) {
                error = fMemoryError; goto ex;
            }
            flags = ttfReader__UShort(r);
            index = ttfReader__UShort(r);
            bHaveInstructions |= (flags & WE_HAVE_INSTRUCTIONS);
            if (flags & ARG_1_AND_2_ARE_WORDS) {
                arg1 = ttfReader__Short(r);
                arg2 = ttfReader__Short(r);
            } else {
                if (flags & ARGS_ARE_XY_VALUES) {
                    /* offsets are signed */
                    arg1 = ttfReader__SignedByte(r);
                    arg2 = ttfReader__SignedByte(r);
                } else { /* anchor points are unsigned */
                    arg1 = ttfReader__Byte(r);
                    arg2 = ttfReader__Byte(r);
                }
            }
            m.b = m.c = m.tx = m.ty = 0;
            if (flags & WE_HAVE_A_SCALE)
                m.a = m.d = (TT_Fixed)ttfReader__Short(r) << 2;
            else if (flags & WE_HAVE_AN_X_AND_Y_SCALE) {
                m.a = (TT_Fixed)ttfReader__Short(r) << 2;
                m.d = (TT_Fixed)ttfReader__Short(r) << 2;
            } else if (flags & WE_HAVE_A_TWO_BY_TWO) {
                m.a = (TT_Fixed)ttfReader__Short(r)<<2;
                m.b = (TT_Fixed)ttfReader__Short(r)<<2;
                m.c = (TT_Fixed)ttfReader__Short(r)<<2;
                m.d = (TT_Fixed)ttfReader__Short(r)<<2;
            } else
                m.a = m.d = 65536;
            e = &usage[nUsage];
            e->m = m;
            e->index = index;
            e->arg1 = arg1;
            e->arg2 = arg2;
            e->flags = flags;
            nUsage++;
        } while (flags & MORE_COMPONENTS);
        if (r->Error(r))
            goto errex;
        nPos = r->Tell(r);
        n_ins = ((!r->Eof(r) && (bHaveInstructions)) ? ttfReader__UShort(r) : 0);
        nPos = r->Tell(r);
        r->ReleaseGlyph(r, glyphIndex);
        glyph = NULL;
        for (i = 0; i < nUsage; i++) {
            ttfGlyphOutline out;
            ttfSubGlyphUsage *e = &usage[i];
            int j;
            TT_Error code;
            int nPointsStored = gOutline->pointCount, nContoursStored = gOutline->contourCount;

            out.contourCount = 0;
            out.pointCount = 0;
            out.bCompound = FALSE;
            pts->org_x += nPointsStored;
            pts->org_y += nPointsStored;
            pts->cur_x += nPointsStored;
            pts->cur_y += nPointsStored;
            pts->touch += nPointsStored;
            pts->contours += nContoursStored;
            tti->usage_top += nUsage;
            code = ttfOutliner__BuildGlyphOutlineAux(self, e->index, m_orig, &out);
            pts->org_x -= nPointsStored;
            pts->org_y -= nPointsStored;
            pts->cur_x -= nPointsStored;
            pts->cur_y -= nPointsStored;
            pts->touch -= nPointsStored;
            tti->usage_top -= nUsage;
            pts->contours -= nContoursStored;
            if (code == fPatented)
                error = code;
            else if (code != fNoError) {
                error = code;
                goto ex;
            }
            if (flags & ARGS_ARE_XY_VALUES) {
                e->m.tx = Scale_X( &exec->metrics, e->arg1 ) << 10;
                e->m.ty = Scale_Y( &exec->metrics, e->arg2 ) << 10;
            } else {
                e->m.tx = (pts->org_x[e->arg1] - pts->org_x[gOutline->pointCount + e->arg2]) << 10;
                e->m.ty = (pts->org_y[e->arg1] - pts->org_y[gOutline->pointCount + e->arg2]) << 10;
            }
            MoveGlyphOutline(pts, nPointsStored, &out, &e->m);
            for (j = nContoursStored; j < out.contourCount + nContoursStored; j++)
                pts->contours[j] += nPointsStored;
            gOutline->contourCount += out.contourCount;
            gOutline->pointCount += out.pointCount;
            if(e->flags & USE_MY_METRICS) {
                gOutline->advance.x = out.advance.x;
                gOutline->sideBearing = out.sideBearing;
            }
        }
        if (execute_bytecode && !skip_instructions && n_ins &&
                !(pFont->inst->GS.instruct_control & 1)) {
            TT_Error code;

            r->LoadGlyph(r, glyphIndex, &glyph, &glyph_size);
            if (r->Error(r))
                goto errex;
            if (nPos + n_ins > glyph_size)
                goto errex;
            code = Set_CodeRange(exec, TT_CodeRange_Glyph, (byte *)glyph + nPos, n_ins);
            if (!code) {
                int k;
                F26Dot6 x;

                nPoints = gOutline->pointCount + 2;
                exec->pts = subglyph.zone;
                pts->n_points = nPoints;
                pts->n_contours = gOutline->contourCount;
                /* add phantom points : */
                pts->org_x[nPoints - 2] = Scale_X(&exec->metrics, subglyph.pp1.x);
                pts->org_y[nPoints - 2] = Scale_Y(&exec->metrics, subglyph.pp1.y);
                pts->org_x[nPoints - 1] = Scale_X(&exec->metrics, subglyph.pp2.x);
                pts->org_y[nPoints - 1] = Scale_Y(&exec->metrics, subglyph.pp2.y);
                pts->touch[nPoints - 1] = 0;
                pts->touch[nPoints - 2] = 0;
                /* if hinting, round the phantom points (not sure) : */
                x = pts->org_x[nPoints - 2];
                x = ((x + 32) & -64) - x;
                if (x)
                    for (k = 0; k < nPoints; k++)
                        pts->org_x[k] += x;
                pts->cur_x[nPoints - 1] = (pts->cur_x[nPoints - 1] + 32) & -64;
                for (k = 0; k < nPoints; k++)
                    pts->touch[k] = pts->touch[k] & TT_Flag_On_Curve;
                org_to_cur(nPoints, pts);
                exec->is_composite = TRUE;
                if (pFont->patented)
                    code = TT_Err_Invalid_Engine;
                else
                    code = Context_Run(exec, FALSE);
                if (!code)
                    cur_to_org(nPoints, pts);
                else if (code == TT_Err_Invalid_Engine)
                    error = fPatented;
                else {
                    /* We have a range of errors that can be caused by
                     * bad bytecode
                     */
                    if ((int)code >= TT_Err_Invalid_Opcode
                     && (int)code <= TT_Err_Invalid_Displacement) {
                        error = fBadInstruction;
                    }
                    else {
                        error = fBadFontData;
                    }
                }
            }
            Unset_CodeRange(exec);
            Clear_CodeRange(exec, TT_CodeRange_Glyph);
        }
    } else if (gOutline->contourCount > 0) {
        int i;
        bool bInsOK;
        byte *onCurve, *stop, flag;
        short *endPoints;
        unsigned int nPos;
        unsigned int n_ins;

        if (self->nContoursTotal + gOutline->contourCount > exec->n_contours) {
            error = fBadFontData; goto ex;
        }
        endPoints = pts->contours;
        for (i = 0; i < gOutline->contourCount; i++)
            endPoints[i] = ttfReader__Short(r);
        for (i = 1; i < gOutline->contourCount; i++)
            if (endPoints[i - 1] < 0 || endPoints[i - 1] >= endPoints[i]) {
                error = fBadFontData; goto ex;
            }
        nPoints = gOutline->pointCount = endPoints[gOutline->contourCount - 1] + 1;
        if (nPoints < 0 || self->nPointsTotal + nPoints + 2 > exec->n_points) {
            error = fBadFontData; goto ex;
        }
        n_ins = ttfReader__Short(r);
        nPos = r->Tell(r);
        r->Seek(r, nPos + n_ins);
        if (r->Error(r))
            goto errex;
        bInsOK = !Set_CodeRange(exec, TT_CodeRange_Glyph, (byte *)glyph + nPos, n_ins);
        onCurve = pts->touch;
        stop = onCurve + gOutline->pointCount;

        while (onCurve < stop) {
            *onCurve++ = flag = ttfReader__Byte(r);
            if (flag & REPEAT_FLAGS) {
                count = ttfReader__Byte(r);
                for (--count; count >= 0 && onCurve < stop; --count)
                    *onCurve++ = flag;
            }
        }
        /*  Lets do X */
        {   short coord = (self->bVertical ? 0 : sideBearing - subglyph.bbox.xMin);
            F26Dot6* x = pts->org_x;
            onCurve = pts->touch;
            while (onCurve < stop) {
                if ((flag = *onCurve++) & XSHORT) {
                    if (flag & SHORT_X_IS_POS)
                        coord += ttfReader__Byte(r);
                    else
                    coord -= ttfReader__Byte(r);
                } else if (!(flag & NEXT_X_IS_ZERO))
                    coord += ttfReader__Short(r);
                *x++ = Scale_X(&exec->metrics, coord);
            }
        }
        /*  Lets do Y */
        {   short coord = 0;
            F26Dot6* y = pts->org_y;
            onCurve = pts->touch;
            while (onCurve < stop) {
                if((flag = *onCurve) & YSHORT)
                    if ( flag & SHORT_Y_IS_POS )
                        coord += ttfReader__Byte(r);
                    else
                        coord -= ttfReader__Byte(r);
                else if(!(flag & NEXT_Y_IS_ZERO))
                    coord += ttfReader__Short(r);
                *y++ = Scale_Y( &exec->metrics, coord );

                /*  Filter off the extra bits */
                *onCurve++ = flag & ONCURVE;
            }
        }
        MoveGlyphOutline(pts, 0, gOutline, m_orig);
        self->nContoursTotal += gOutline->contourCount;
        self->nPointsTotal += nPoints;
        if (execute_bytecode && !skip_instructions &&
                !r->Error(r) && n_ins && bInsOK && !(pFont->inst->GS.instruct_control & 1)) {
            TGlyph_Zone *pts = &exec->pts;
            int k;
            F26Dot6 x;
            TT_Error code;

            exec->is_composite = FALSE;
            /* add phantom points : */
            pts->org_x[nPoints    ] = Scale_X(&exec->metrics, subglyph.pp1.x);
            pts->org_y[nPoints    ] = Scale_Y(&exec->metrics, subglyph.pp1.y);
            pts->org_x[nPoints + 1] = Scale_X(&exec->metrics, subglyph.pp2.x);
            pts->org_y[nPoints + 1] = Scale_Y(&exec->metrics, subglyph.pp2.y);
            pts->touch[nPoints    ] = 0;
            pts->touch[nPoints + 1] = 0;
            pts->n_points   = nPoints + 2;
            pts->n_contours = gOutline->contourCount;
            /* if hinting, round the phantom points (not sure) : */
            x = pts->org_x[nPoints];
            x = ((x + 32) & -64) - x;
            if (x)
                for (k = 0; k < nPoints + 2; k++)
                    pts->org_x[k] += x;
            org_to_cur(nPoints + 2, pts);
            exec->is_composite = FALSE;
            for (k = 0; k < nPoints + 2; k++)
                pts->touch[k] &= TT_Flag_On_Curve;
            if (pFont->patented)
                code = TT_Err_Invalid_Engine;
            else
                code = Context_Run(exec, FALSE );
            if (!code)
                cur_to_org(nPoints + 2, pts);
            else if (code == TT_Err_Invalid_Engine)
                error = fPatented;
            else
                error = fBadInstruction;
            gOutline->sideBearing = subglyph.bbox.xMin - subglyph.pp1.x;
            gOutline->advance.x = subglyph.pp2.x - subglyph.pp1.x;
        }
        Unset_CodeRange(exec);
        Clear_CodeRange(exec, TT_CodeRange_Glyph);
    } else
        error = fBadFontData;
    goto ex;
errex:;
    error = fBadFontData;
ex:;
    r->ReleaseGlyph(r, glyphIndex);

    if (error == fBadInstruction && execute_bytecode) {
        /* reset a load of stuff so we can try again without hinting */
        exec = pFont->exec;
        pts = &exec->pts;
        usage = tti->usage + tti->usage_top;
        glyph = NULL;
        self->nPointsTotal -= (nPoints + 2);
        nPoints = 0;
        self->nContoursTotal -= gOutline->contourCount;
        error = fNoError;
        execute_bytecode = false;
        r->Seek(r, nPosBeg);
        goto retry;
    }

    return error;
}

static FontError ttfOutliner__BuildGlyphOutline(ttfOutliner *self, int glyphIndex,
            float orig_x, float orig_y, ttfGlyphOutline* gOutline)
{
    FixMatrix m_orig = {1 << 16, 0, 0, 1 << 16, 0, 0};

    /* Round towards zero like old character coordinate conversions do. */
    m_orig.tx = floatToF16Dot16(orig_x);
    m_orig.ty = floatToF16Dot16(orig_y);
    return ttfOutliner__BuildGlyphOutlineAux(self, glyphIndex, &m_orig, gOutline);
}

#define AVECTOR_BUG 1 /* Work around a bug in AVector fonts. */

void ttfOutliner__DrawGlyphOutline(ttfOutliner *self)
{   ttfGlyphOutline* out = &self->out;
    FloatMatrix *m = &self->post_transform;
    ttfFont *pFont = self->pFont;
    ttfExport *exp = self->exp;
    TExecution_Context *exec = pFont->exec;
    TGlyph_Zone *epts = &exec->pts;
    short* endP = epts->contours;
    byte* onCurve = epts->touch;
    F26Dot6* x = epts->org_x;
    F26Dot6* y = epts->org_y;
    F26Dot6 px, py;
    short sp, ctr;
    FloatPoint p0, p1, p2, p3;
#   if AVECTOR_BUG
    F26Dot6 expand_x = Scale_X(&exec->metrics, pFont->nUnitsPerEm * 2);
    F26Dot6 expand_y = Scale_Y(&exec->metrics, pFont->nUnitsPerEm * 2);
    F26Dot6 xMin = out->xMinB - expand_x, xMax = out->xMaxB + expand_x;
    F26Dot6 yMin = out->yMinB - expand_y, yMax = out->yMaxB + expand_y;
#   endif

    TransformF26Dot6PointFloat(&p1, out->advance.x, out->advance.y, m);
    p1.x -= self->post_transform.tx;
    p1.y -= self->post_transform.ty;
    exp->SetWidth(exp, &p1);
    sp = -1;
    for (ctr = out->contourCount; ctr != 0; --ctr) {
        short pt, pts = *endP - sp;
        short ep = pts - 1;

        if (pts < 3) {
            x += pts;
            y += pts;
            onCurve += pts;
            sp = *endP++;
            continue;   /* skip 1 and 2 point contours */
        }

        if (exp->bPoints) {
            for (pt = 0; pt <= ep; pt++) {
                px = x[pt], py = y[pt];
#		if AVECTOR_BUG
                    if (x[pt] < xMin || xMax < x[pt] || y[pt] < yMin || yMax < y[pt]) {
                        short prevIndex = pt == 0 ? ep : pt - 1;
                        short nextIndex = pt == ep ? 0 : pt + 1;
                        if (nextIndex > ep)
                            nextIndex = 0;
                        px=AVE(x[prevIndex], x[nextIndex]);
                        py=AVE(y[prevIndex], y[nextIndex]);
                    }
#		endif
                TransformF26Dot6PointFloat(&p0, px, py, m);
                exp->Point(exp, &p0, onCurve[pt], !pt);
            }
        }

        if (exp->bOutline) {
            pt = 0;
            if(onCurve[ep] & 1) {
                px = x[ep];
                py = y[ep];
            } else if (onCurve[0] & 1) {
                px = x[0];
                py = y[0];
                pt = 1;
            } else {
                px = AVE(x[0], x[ep]);
                py = AVE(y[0], y[ep]);
            }
            self->ppx = px; self->ppy = py;
            TransformF26Dot6PointFloat(&p0, px, py, m);
            exp->MoveTo(exp, &p0);

            for (; pt <= ep; pt++) {
                short prevIndex = pt == 0 ? ep : pt - 1;
                short nextIndex = pt == ep ? 0 : pt + 1;
                if (onCurve[pt] & 1) {
                    if (onCurve[prevIndex] & 1) {
                        px = x[pt];
                        py = y[pt];
                        if (self->ppx != px || self->ppy != py) {
                            TransformF26Dot6PointFloat(&p1, px, py, m);
                            exp->LineTo(exp, &p1);
                            self->ppx = px; self->ppy = py;
                            p0 = p1;
                        }
                    }
                } else {
                    F26Dot6 prevX, prevY, nextX, nextY;

                    px = x[pt];
                    py = y[pt];
#		    if AVECTOR_BUG
                        if(x[pt] < xMin || xMax < x[pt] || y[pt] < yMin || yMax < y[pt]) {
                            px=AVE(x[prevIndex], x[nextIndex]);
                            py=AVE(y[prevIndex], y[nextIndex]);
                        }
#		    endif
                    if (onCurve[prevIndex] & 1) {
                        prevX = x[prevIndex];
                        prevY = y[prevIndex];
                    } else {
                        prevX = AVE(x[prevIndex], px);
                        prevY = AVE(y[prevIndex], py);
                    }
                    if (onCurve[nextIndex] & 1) {
                        nextX = x[nextIndex];
                        nextY = y[nextIndex];
                    } else {
                        nextX = AVE(px, x[nextIndex]);
                        nextY = AVE(py, y[nextIndex]);
                    }
                    if (self->ppx != nextX || self->ppy != nextY) {
                        double dx1, dy1, dx2, dy2, dx3, dy3;
                        const double prec = 1e-6;

                        TransformF26Dot6PointFloat(&p1, (prevX + (px << 1)) / 3, (prevY + (py << 1)) / 3, m);
                        TransformF26Dot6PointFloat(&p2, (nextX + (px << 1)) / 3, (nextY + (py << 1)) / 3, m);
                        TransformF26Dot6PointFloat(&p3, nextX, nextY, m);
                        dx1 = p1.x - p0.x, dy1 = p1.y - p0.y;
                        dx2 = p2.x - p0.x, dy2 = p2.y - p0.y;
                        dx3 = p3.x - p0.x, dy3 = p3.y - p0.y;
                        if (fabs(dx1 * dy3 - dy1 * dx3) > prec * fabs(dx1 * dx3 - dy1 * dy3) ||
                            fabs(dx2 * dy3 - dy2 * dx3) > prec * fabs(dx2 * dx3 - dy2 * dy3))
                            exp->CurveTo(exp, &p1, &p2, &p3);
                        else
                            exp->LineTo(exp, &p3);
                        self->ppx = nextX; self->ppy = nextY;
                        p0 = p3;
                    }
                }
            }
            exp->Close(exp);
        }
        x += pts;
        y += pts;
        onCurve += pts;
        sp = *endP++;
    }
}

FontError ttfOutliner__Outline(ttfOutliner *self, int glyphIndex,
        float orig_x, float orig_y, FloatMatrix *m1)
{   ttfFont *pFont = self->pFont;
    FontError error;

    self->post_transform = *m1;
    self->out.contourCount = 0;
    self->out.pointCount = 0;
    self->out.bCompound = FALSE;
    self->nPointsTotal = 0;
    self->nContoursTotal = 0;
    self->out.advance.x = self->out.advance.y = 0;
    ttfFont__StartGlyph(pFont);
    error = ttfOutliner__BuildGlyphOutline(self, glyphIndex, orig_x, orig_y, &self->out);
    ttfFont__StopGlyph(pFont);
    if (pFont->nUnitsPerEm <= 0)
        pFont->nUnitsPerEm = 1024;
    if (pFont->design_grid) {
        self->post_transform.a /= pFont->nUnitsPerEm;
        self->post_transform.b /= pFont->nUnitsPerEm;
        self->post_transform.c /= pFont->nUnitsPerEm;
        self->post_transform.d /= pFont->nUnitsPerEm;
    }
    return error;
}
