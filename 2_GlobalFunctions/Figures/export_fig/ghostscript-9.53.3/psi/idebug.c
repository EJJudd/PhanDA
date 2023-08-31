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


/* Debugging support for Ghostscript interpreter */
/* This file must always be compiled with DEBUG set. */
#undef DEBUG
#define DEBUG
#include "string_.h"
#include "ghost.h"
#include "gxalloc.h"		/* for procs for getting struct type */
#include "idebug.h"		/* defines interface */
#include "idict.h"
#include "iname.h"
#include "ipacked.h"
#include "istack.h"
#include "iutil.h"
#include "ivmspace.h"
#include "opdef.h"

/* Table of type name strings */
static const char *const type_strings[] = {
    REF_TYPE_DEBUG_PRINT_STRINGS
};

/* First unassigned type index */
extern const int tx_next_index;	/* in interp.c */

/* Print a name. */
void
debug_print_name(const gs_memory_t *mem, const ref * pnref)
{
    ref sref;

    name_string_ref(mem, pnref, &sref);
    debug_print_string(mem, sref.value.const_bytes, r_size(&sref));
}
void
debug_print_name_index(const gs_memory_t *mem, name_index_t nidx)
{
    ref nref;

    name_index_ref(mem, nidx, &nref);
    debug_print_name(mem, &nref);
}

/* Print a ref. */
static void
debug_print_full_ref(const gs_memory_t *mem, const ref * pref)
{
    uint size = r_size(pref);
    ref nref;

    dmprintf1(mem, "(%x)", r_type_attrs(pref));
    switch (r_type(pref)) {
        case t_array:
            dmprintf2(mem, "array(%u)"PRI_INTPTR, size, (intptr_t)pref->value.refs);
            break;
        case t_astruct:
            goto strct;
        case t_boolean:
            dmprintf1(mem, "boolean %x", pref->value.boolval);
            break;
        case t_device:
            dmprintf1(mem, "device "PRI_INTPTR, (intptr_t) pref->value.pdevice);
            break;
        case t_dictionary:
            dmprintf3(mem, "dict(%u/%u)"PRI_INTPTR,
                     dict_length(pref), dict_maxlength(pref),
                     (intptr_t)pref->value.pdict);
            break;
        case t_file:
            dmprintf1(mem, "file "PRI_INTPTR, (intptr_t)pref->value.pfile);
            break;
        case t_fontID:
            goto strct;
        case t_integer:
            dmprintf1(mem, "int %"PRIpsint, pref->value.intval);
            break;
        case t_mark:
            dmprintf(mem, "mark");
            break;
        case t_mixedarray:
            dmprintf2(mem, "mixed packedarray(%u)"PRI_INTPTR"", size,
                     (intptr_t)pref->value.packed);
            break;
        case t_name:
            dmprintf2(mem, "name("PRI_INTPTR"#%u)", (intptr_t)pref->value.pname,
                     name_index(mem, pref));
            debug_print_name(mem, pref);
            break;
        case t_null:
            dmprintf(mem, "null");
            break;
        case t_oparray:
            dmprintf2(mem, "op_array(%u)"PRI_INTPTR":", size, (intptr_t) pref->value.const_refs);
            {
                const op_array_table *opt = get_op_array(mem, size);

                name_index_ref(mem, opt->nx_table[size - opt->base_index], &nref);
            }
            debug_print_name(mem, &nref);
            break;
        case t_operator:
            dmprintf1(mem, "op(%u", size);
            if (size > 0 && size < op_def_count)	/* just in case */
                dmprintf1(mem, ":%s", (const char *)(op_index_def(size)->oname + 1));
            dmprintf1(mem, ")"PRI_INTPTR"", (intptr_t)pref->value.opproc);
            break;
        case t_real:
            dmprintf1(mem, "real %f", pref->value.realval);
            break;
        case t_save:
            dmprintf1(mem, "save %lu", pref->value.saveid);
            break;
        case t_shortarray:
            dmprintf2(mem, "short packedarray(%u)"PRI_INTPTR"", size,
                     (intptr_t)pref->value.packed);
            break;
        case t_string:
            dmprintf2(mem, "string(%u)"PRI_INTPTR"", size, (intptr_t)pref->value.bytes);
            break;
        case t_struct:
          strct:{
                obj_header_t *obj = (obj_header_t *) pref->value.pstruct;
                /* HACK: We know this object was allocated with gsalloc.c. */
                gs_memory_type_ptr_t otype =
                    gs_ref_memory_procs.object_type(NULL, obj);

                dmprintf2(mem, "struct %s "PRI_INTPTR"",
                         (r_is_foreign(pref) ? "-foreign-" :
                          gs_struct_type_name_string(otype)),
                         (intptr_t)obj);
            }
            break;
        default:
            dmprintf1(mem, "type 0x%x", r_type(pref));
    }
}
static void
debug_print_packed_ref(const gs_memory_t *mem, const ref_packed *pref)
{
    ushort elt = *pref & packed_value_mask;
    ref nref;

    switch (*pref >> r_packed_type_shift) {
        case pt_executable_operator:
            dmprintf(mem, "<op_name>");
            op_index_ref(mem, elt, &nref);
            debug_print_ref(mem, &nref);
            break;
        case pt_integer:
            dmprintf1(mem, "<int> %d", (int)elt + packed_min_intval);
            break;
        case pt_literal_name:
            dmprintf(mem, "<lit_name>");
            goto ptn;
        case pt_executable_name:
            dmprintf(mem, "<exec_name>");
    ptn:    name_index_ref(mem, elt, &nref);
            dmprintf2(mem, "("PRI_INTPTR"#%u)", (intptr_t) nref.value.pname, elt);
            debug_print_name(mem, &nref);
            break;
        default:
            dmprintf2(mem, "<packed_%d?>0x%x", *pref >> r_packed_type_shift, elt);
    }
}
void
debug_print_ref_packed(const gs_memory_t *mem, const ref_packed *rpp)
{
    if (r_is_packed(rpp))
        debug_print_packed_ref(mem, rpp);
    else
        debug_print_full_ref(mem, (const ref *)rpp);
    dmflush(mem);
}
void
debug_print_ref(const gs_memory_t *mem, const ref * pref)
{
    debug_print_ref_packed(mem, (const ref_packed *)pref);
}

/* Dump one ref. */
static void print_ref_data(const gs_memory_t *mem, const ref *);
void
debug_dump_one_ref(const gs_memory_t *mem, const ref * p)
{
    uint attrs = r_type_attrs(p);
    uint type = r_type(p);
    static const ref_attr_print_mask_t apm[] = {
        REF_ATTR_PRINT_MASKS,
        {0, 0, 0}
    };
    const ref_attr_print_mask_t *ap = apm;

    if (type >= tx_next_index)
        dmprintf1(mem, "0x%02x?? ", type);
    else if (type >= t_next_index)
        dmprintf(mem, "opr* ");
    else
        dmprintf1(mem, "%s ", type_strings[type]);
    for (; ap->mask; ++ap)
        if ((attrs & ap->mask) == ap->value)
            dmputc(mem, ap->print);
    dmprintf2(mem, " 0x%04x 0x%08lx", r_size(p), *(const ulong *)&p->value);
    print_ref_data(mem, p);
    dmflush(mem);
}
static void
print_ref_data(const gs_memory_t *mem, const ref *p)
{
#define BUF_SIZE 30
    byte buf[BUF_SIZE + 1];
    const byte *pchars;
    uint plen;

    if (obj_cvs(mem, p, buf, countof(buf) - 1, &plen, &pchars) >= 0 &&
        pchars == buf &&
        ((buf[plen] = 0), strcmp((char *)buf, "--nostringval--"))
        )
        dmprintf1(mem, " = %s", (char *)buf);
#undef BUF_SIZE
}

/* Dump a region of memory containing refs. */
void
debug_dump_refs(const gs_memory_t *mem, const ref * from,
                uint size, const char *msg)
{
    const ref *p = from;
    uint count = size;

    if (size && msg)
        dmprintf2(mem, "%s at "PRI_INTPTR":\n", msg, (intptr_t)from);
    while (count--) {
        dmprintf2(mem, PRI_INTPTR": 0x%04x ", (intptr_t)p, r_type_attrs(p));
        debug_dump_one_ref(mem, p);
        dmputc(mem, '\n');
        p++;
    }
}

/* Dump a stack. */
void
debug_dump_stack(const gs_memory_t *mem,
                 const ref_stack_t * pstack, const char *msg)
{
    int i;
    const char *m = msg;

    for (i = ref_stack_count(pstack); i > 0;) {
        const ref *p = ref_stack_index(pstack, --i);

        if (m) {
            dmprintf2(mem, "%s at "PRI_INTPTR":\n", m, (intptr_t)pstack);
            m = NULL;
        }
        dmprintf2(mem, PRI_INTPTR": 0x%02x ", (intptr_t)p, r_type(p));
        debug_dump_one_ref(mem, p);
        dmputc(mem, '\n');
    }
}

/* Dump an array. */
void
debug_dump_array(const gs_memory_t *mem, const ref * array)
{
    const ref_packed *pp;
    uint type = r_type(array);
    uint len;

    switch (type) {
        default:
            dmprintf2(mem, "%s at "PRI_INTPTR" isn't an array.\n",
                      (type < countof(type_strings) ?
                       type_strings[type] : "????"),
                      (intptr_t)array);
            return;
        case t_oparray:
            /* This isn't really an array, but we'd like to see */
            /* its contents anyway. */
            debug_dump_array(mem, array->value.const_refs);
            return;
        case t_array:
        case t_mixedarray:
        case t_shortarray:
            ;
    }

    /* This "packed" loop works for all array-types. */
    for (len = r_size(array), pp = array->value.packed;
         len > 0;
         len--, pp = packed_next(pp)) {
        ref temp;

        packed_get(mem, pp, &temp);
        if (r_is_packed(pp)) {
            dmprintf2(mem, PRI_INTPTR"* 0x%04x ", (intptr_t)pp, (uint)*pp);
            print_ref_data(mem, &temp);
        } else {
            dmprintf2(mem, PRI_INTPTR": 0x%02x ", (intptr_t)pp, r_type(&temp));
            debug_dump_one_ref(mem, &temp);
        }
        dmputc(mem, '\n');
    }
}
