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


/* Image operators */
#include "math_.h"
#include "memory_.h"
#include "stat_.h" /* get system header early to avoid name clash on Cygwin */
#include "ghost.h"
#include "oper.h"
#include "gscolor.h"
#include "gscspace.h"
#include "gscolor2.h"
#include "gsmatrix.h"
#include "gsimage.h"
#include "gxfixed.h"
#include "gsstruct.h"
#include "gxiparam.h"
#include "idict.h"
#include "idparam.h"
#include "estack.h"		/* for image[mask] */
#include "ialloc.h"
#include "igstate.h"
#include "ilevel.h"
#include "store.h"
#include "stream.h"
#include "ifilter.h"		/* for stream exception handling */
#include "iimage.h"
#include "gxcspace.h"

/* Forward references */
static int zimage_data_setup(i_ctx_t *i_ctx_p, const gs_pixel_image_t * pim,
                                 gx_image_enum_common_t * pie,
                                 const ref * sources, int npop);
static int image_proc_process(i_ctx_t *);
static int image_file_continue(i_ctx_t *);
static int image_string_continue(i_ctx_t *);
static int image_cleanup(i_ctx_t *);

/* Extract and check the parameters for a gs_data_image_t. */
int
data_image_params(const gs_memory_t *mem,
                  const ref *op, gs_data_image_t *pim,
                  image_params *pip, bool require_DataSource,
                  int num_components, int max_bits_per_component,
                  bool islab)
{
    int code;
    ref *pds;

    check_type(*op, t_dictionary);
    check_dict_read(*op);
    code = dict_int_param(op, "Width", 0, max_int_in_fixed/2, -1, &pim->Width);
    if (code < 0)
    	return code;
    code = dict_int_param(op, "Height", 0, max_int_in_fixed/2, -1, &pim->Height);
    if (code < 0)
    	return code;
    code = dict_matrix_param(mem, op, "ImageMatrix", &pim->ImageMatrix);
    if (code < 0)
    	return code;
    code = dict_bool_param(op, "MultipleDataSources", false, &pip->MultipleDataSources);
    if (code < 0)
    	return code;
    code = dict_int_param(op, "BitsPerComponent", 1, max_bits_per_component, -1, &pim->BitsPerComponent);
    if (code < 0)
    	return code;
    code = dict_bool_param(op, "Interpolate", false, &pim->Interpolate);
    if (code < 0)
    	return code;

    /* Decode size pulled out of here to catch case of Lab color space which
       has a 4 entry range.  We also do NOT want to do Lab decoding IF range
       is the common -128 127 for a and b. Otherwise we end up doing multiple
       decode operations, since ICC flow will expect encoded data.
       That is resolved later.  Also discovered that PDF write will stick
       6 entry range in wich appears to be OK as far as AR is concerned so
       we have to handle that too. */
    if (islab) {
        /* Note that it is possible that only the ab range values are there
           or the lab values.  I have seen both cases.... */
        code = dict_floats_param(mem, op, "Decode", 4,
                                 &pim->Decode[2], NULL);
        if (code < 0) {
            /* Try for all three */
            code = dict_floats_param(mem, op, "Decode", 6,
                                                        &pim->Decode[0], NULL);
        } else {
            /* Set the range on the L */
            pim->Decode[0] = 0;
            pim->Decode[1] = 100.0;
        }
        if (code < 0)
            return code;
    } else {
        code = dict_floats_param(mem, op, "Decode",
                                                    num_components * 2,
                                                    &pim->Decode[0], NULL);
        if (code < 0)
            return code;
    }
    pip->pDecode = &pim->Decode[0];
    /* Extract and check the data sources. */
    if ((code = dict_find_string(op, "DataSource", &pds)) <= 0) {
        if (require_DataSource)
            return (code < 0 ? code : gs_note_error(gs_error_rangecheck));
        return 1;		/* no data source */
    }
    if (pip->MultipleDataSources) {
        ref *ds = pip->DataSource;
        long i;
        if (!r_is_array(pds))
            return_error(gs_error_typecheck);
        if (r_size(pds) != num_components)
            return_error(gs_error_rangecheck);
        for (i = 0; i < num_components; ++i)
            array_get(mem, pds, i, &ds[i]);
        if (r_type(&ds[0]) == t_string) {
            /* We don't have a problem with the strings of different length
             * but Adobe does and CET tast 12-02.ps reports this as an error.
             */
            for (i = 1; i < num_components; ++i) {
                if (r_type(&ds[i]) == t_string && r_size(&ds[i]) != r_size(&ds[0])) {
                    return_error(gs_error_rangecheck);
                }
            }
        }
    } else
        pip->DataSource[0] = *pds;
    return 0;
}

/* Extract and check the parameters for a gs_pixel_image_t. */
int
pixel_image_params(i_ctx_t *i_ctx_p, const ref *op, gs_pixel_image_t *pim,
                   image_params *pip, int max_bits_per_component,
                   gs_color_space *csp)
{
    bool islab = false;
    int num_components =
        gs_color_space_num_components(csp);
    int code;

    if (num_components < 1)
        return_error(gs_error_rangecheck);	/* Pattern space not allowed */
    pim->ColorSpace = csp;

    if (pim->ColorSpace->cmm_icc_profile_data != NULL)
        islab = pim->ColorSpace->cmm_icc_profile_data->islab;

    code = data_image_params(imemory, op, (gs_data_image_t *) pim, pip, true,
                             num_components, max_bits_per_component, islab);
    if (code < 0)
        return code;
    pim->format =
        (pip->MultipleDataSources ? gs_image_format_component_planar :
         gs_image_format_chunky);
    return dict_bool_param(op, "CombineWithColor", false,
                           &pim->CombineWithColor);
}

/* Common setup for all Level 1 and 2 images, and ImageType 4 images. */
int
zimage_setup(i_ctx_t *i_ctx_p, const gs_pixel_image_t * pim,
             const ref * sources, bool uses_color, int npop)
{
    gx_image_enum_common_t *pie;
    int code =
        gs_image_begin_typed((const gs_image_common_t *)pim, igs,
                             uses_color, false, &pie);

    if (code < 0)
        return code;
    return zimage_data_setup(i_ctx_p, (const gs_pixel_image_t *)pim, pie,
                             sources, npop);
}

/* <dict> .image1 - */
static int
zimage1(i_ctx_t *i_ctx_p)
{
    os_ptr          op = osp;
    gs_image_t      image;
    image_params    ip;
    int             code;
    gs_color_space *csp = gs_currentcolorspace(igs);

    /* Adobe interpreters accept sampled images when the current color
     * space is a pattern color space using the base color space instead
     * of the pattern space. CET 12-07a-12
     * If all conditions are not met the pattern color space goes through
     * triggering a rangecheck error.
     */
    if (gs_currentcpsimode(imemory) && gs_color_space_num_components(csp) < 1) {
       gs_color_space *bsp = csp->base_space;
       if (bsp)
         csp = bsp;
    }

    gs_image_t_init(&image, csp);
    code = pixel_image_params( i_ctx_p,
                               op,
                               (gs_pixel_image_t *)&image,
                               &ip,
                               (level2_enabled ? 16 : 8),
                               csp);
    if (code < 0)
        return code;

    image.Alpha = gs_image_alpha_none;
        /* swap Width, Height, and ImageMatrix so that it comes out the same */
        /* This is only for performance, so only do it for non-skew cases */
    if (image.Width == 1 && image.Height > 1 && image.BitsPerComponent == 8 &&
        image.ImageMatrix.xy == 0.0 && image.ImageMatrix.yx == 0.0 &&
        image.ImageMatrix.tx == 0.0) {
        float ftmp;
        int   itemp;

        itemp = image.Width;
        image.Width = image.Height;
        image.Height = itemp;

        image.ImageMatrix.xy = image.ImageMatrix.xx;
        image.ImageMatrix.yx = image.ImageMatrix.yy;
        image.ImageMatrix.xx = 0.0;
        image.ImageMatrix.yy = 0.0;
        ftmp = image.ImageMatrix.tx;
        image.ImageMatrix.tx = image.ImageMatrix.ty;
        image.ImageMatrix.ty = ftmp;
    }
    return zimage_setup( i_ctx_p,
                         (gs_pixel_image_t *)&image,
                         &ip.DataSource[0],
                         image.CombineWithColor,
                         1 );
}

/* <dict> .imagemask1 - */
static int
zimagemask1(i_ctx_t *i_ctx_p)
{
    os_ptr op = osp;
    gs_image_t image;
    image_params ip;
    int code;

    gs_image_t_init_mask_adjust(&image, false,
                                gs_incachedevice(igs) != CACHE_DEVICE_NONE);
    code = data_image_params(imemory, op, (gs_data_image_t *) & image,
                             &ip, true, 1, 1, false);
    if (code < 0)
        return code;
    return zimage_setup(i_ctx_p, (gs_pixel_image_t *)&image, &ip.DataSource[0],
                        true, 1);
}

/* Common setup for all Level 1 and 2 images, and ImageType 3 and 4 images. */
/*
 * We push the following on the estack.
 *      control mark,
 *	num_sources,
 *      for I = num_sources-1 ... 0:
 *          data source I,
 *          aliasing information:
 *              if source is not file, 1, except that the topmost value
 *		  is used for bookkeeping in the procedure case (see below);
 *              if file is referenced by a total of M different sources and
 *                this is the occurrence with the lowest I, M;
 *              otherwise, -J, where J is the lowest I of the same file as
 *                this one;
 *      current plane index,
 *      num_sources,
 *      enumeration structure.
 */
#define NUM_PUSH(nsource) ((nsource) * 2 + 5)
/*
 * We can access these values either from the bottom (esp at control mark - 1,
 * EBOT macros) or the top (esp = enumeration structure, ETOP macros).
 * Note that all macros return pointers.
 */
#define EBOT_NUM_SOURCES(ep) ((ep) + 2)
#define EBOT_SOURCE(ep, i)\
  ((ep) + 3 + (EBOT_NUM_SOURCES(ep)->value.intval - 1 - (i)) * 2)
#define ETOP_SOURCE(ep, i)\
  ((ep) - 4 - (i) * 2)
#define ETOP_PLANE_INDEX(ep) ((ep) - 2)
#define ETOP_NUM_SOURCES(ep) ((ep) - 1)
static int
zimage_data_setup(i_ctx_t *i_ctx_p, const gs_pixel_image_t * pim,
                  gx_image_enum_common_t * pie, const ref * sources, int npop)
{
    int num_sources = pie->num_planes;
    int inumpush = NUM_PUSH(num_sources);
    int code;
    gs_image_enum *penum;
    int px;
    const ref *pp;
    bool string_sources = true;

    check_estack(inumpush + 2);	/* stuff above, + continuation + proc */
    make_int(EBOT_NUM_SOURCES(esp), num_sources);
    /*
     * Note that the data sources may be procedures, strings, or (Level
     * 2 only) files.  (The Level 1 reference manual says that Level 1
     * requires procedures, but Adobe Level 1 interpreters also accept
     * strings.)  The sources must all be of the same type.
     *
     * The Adobe documentation explicitly says that if two or more of the
     * data sources are the same or inter-dependent files, the result is not
     * defined.  We don't have a problem with the bookkeeping for
     * inter-dependent files, since each one has its own buffer, but we do
     * have to be careful if two or more sources are actually the same file.
     * That is the reason for the aliasing information described above.
     */
    for (px = 0, pp = sources; px < num_sources; px++, pp++) {
        es_ptr ep = EBOT_SOURCE(esp, px);

        make_int(ep + 1, 1);	/* default is no aliasing */
        switch (r_type(pp)) {
            case t_file:
                if (!level2_enabled)
                    return_error(gs_error_typecheck);
                /* Check for aliasing. */
                {
                    int pi;

                    for (pi = 0; pi < px; ++pi)
                        if (sources[pi].value.pfile == pp->value.pfile) {
                            /* Record aliasing */
                            make_int(ep + 1, -pi);
                            EBOT_SOURCE(esp, pi)[1].value.intval++;
                            break;
                        }
                }
                string_sources = false;
                /* falls through */
            case t_string:
                if (r_type(pp) != r_type(sources)) {
                    gx_image_end(pie, false);    /* Clean up pie */
                    return_error(gs_error_typecheck);
                }
                check_read(*pp);
                break;
            default:
                if (!r_is_proc(sources)) {
                    static const char ds[] = "DataSource";
                    gx_image_end(pie, false);    /* Clean up pie */
                    gs_errorinfo_put_pair(i_ctx_p, ds, sizeof(ds) - 1, pp);
                    return_error(gs_error_typecheck);
                }
                check_proc(*pp);
                string_sources = false;
        }
        *ep = *pp;
    }
    /* Always place the image enumerator into local memory,
       because pie may have local objects inherited from igs,
       which may be local when the current allocation mode is global.
       Bug 688140. */
    if ((penum = gs_image_enum_alloc(imemory_local, "image_setup")) == 0)
        return_error(gs_error_VMerror);
    code = gs_image_enum_init(penum, pie, (const gs_data_image_t *)pim, igs);
    if (code != 0 || (pie->skipping && string_sources)) {		/* error, or empty image */
        int code1 = gs_image_cleanup_and_free_enum(penum, igs);

        if (code >= 0)		/* empty image */
            pop(npop);
        if (code >= 0 && code1 < 0)
            code = code1;
        return code;
    }
    push_mark_estack(es_other, image_cleanup);
    esp += inumpush - 1;
    make_int(ETOP_PLANE_INDEX(esp), 0);
    make_int(ETOP_NUM_SOURCES(esp), num_sources);
    make_struct(esp, avm_local, penum);
    switch (r_type(sources)) {
        case t_file:
            push_op_estack(image_file_continue);
            break;
        case t_string:
            push_op_estack(image_string_continue);
            break;
        default:		/* procedure */
            push_op_estack(image_proc_process);
            break;
    }
    pop(npop);
    return o_push_estack;
}
/* Pop all the control information off the e-stack. */
static es_ptr
zimage_pop_estack(es_ptr tep)
{
    return tep - NUM_PUSH(ETOP_NUM_SOURCES(tep)->value.intval);
}

/*
 * Continuation for procedure data source.  We use the topmost aliasing slot
 * to remember whether we've just called the procedure (1) or whether we're
 * returning from a RemapColor callout (0).
 */
static int
image_proc_continue(i_ctx_t *i_ctx_p)
{
    os_ptr op = osp;
    gs_image_enum *penum = r_ptr(esp, gs_image_enum);
    int px = ETOP_PLANE_INDEX(esp)->value.intval;
    int num_sources = ETOP_NUM_SOURCES(esp)->value.intval;
    uint size, used[GS_IMAGE_MAX_COMPONENTS];
    gs_const_string plane_data[GS_IMAGE_MAX_COMPONENTS];
    const byte *wanted;
    int i, code;

    if (!r_has_type_attrs(op, t_string, a_read)) {
        check_op(1);
        /* Procedure didn't return a (readable) string.  Quit. */
        esp = zimage_pop_estack(esp);
        image_cleanup(i_ctx_p);
        return_error(!r_has_type(op, t_string) ? gs_error_typecheck : gs_error_invalidaccess);
    }
    size = r_size(op);
    if (size == 0 && ETOP_SOURCE(esp, 0)[1].value.intval == 0)
        code = 1;
    else {
        for (i = 0; i < num_sources; i++)
            plane_data[i].size = 0;
        plane_data[px].data = op->value.bytes;
        plane_data[px].size = size;
        code = gs_image_next_planes(penum, plane_data, used);
        if (code == gs_error_Remap_Color) {
            op->value.bytes += used[px]; /* skip used data */
            r_dec_size(op, used[px]);
            ETOP_SOURCE(esp, 0)[1].value.intval = 0; /* RemapColor callout */
            return code;
        }
    }
    if (code) {			/* Stop now. */
        esp = zimage_pop_estack(esp);
        pop(1);
        image_cleanup(i_ctx_p);
        return (code < 0 ? code : o_pop_estack);
    }
    pop(1);
    wanted = gs_image_planes_wanted(penum);
    do {
        if (++px == num_sources)
            px = 0;
    } while (!wanted[px]);
    ETOP_PLANE_INDEX(esp)->value.intval = px;
    return image_proc_process(i_ctx_p);
}
static int
image_proc_process(i_ctx_t *i_ctx_p)
{
    int px = ETOP_PLANE_INDEX(esp)->value.intval;
    gs_image_enum *penum = r_ptr(esp, gs_image_enum);
    const byte *wanted = gs_image_planes_wanted(penum);
    int num_sources = ETOP_NUM_SOURCES(esp)->value.intval;
    const ref *pp;

    ETOP_SOURCE(esp, 0)[1].value.intval = 0; /* procedure callout */
    while (!wanted[px]) {
        if (++px == num_sources)
            px = 0;
        ETOP_PLANE_INDEX(esp)->value.intval = px;
    }
    pp = ETOP_SOURCE(esp, px);
    push_op_estack(image_proc_continue);
    *++esp = *pp;
    return o_push_estack;
}

/* Continue processing data from an image with file data sources. */
static int
image_file_continue(i_ctx_t *i_ctx_p)
{
    gs_image_enum *penum = r_ptr(esp, gs_image_enum);
    int num_sources = ETOP_NUM_SOURCES(esp)->value.intval;

    for (;;) {
        uint min_avail = max_int;
        gs_const_string plane_data[GS_IMAGE_MAX_COMPONENTS];
        int code;
        int px;
        const ref *pp;
        int at_eof_count = 0;
        int total_used;

        /*
         * Do a first pass through the files to ensure that at least
         * one has data available in its buffer.
         */

        for (px = 0, pp = ETOP_SOURCE(esp, 0); px < num_sources;
             ++px, pp -= 2
            ) {
            int num_aliases = pp[1].value.intval;
            stream *s = pp->value.pfile;
            int min_left;
            uint avail;

            if (num_aliases <= 0)
                num_aliases = ETOP_SOURCE(esp, -num_aliases)[1].value.intval;
            while ((avail = sbufavailable(s)) <=
                   (min_left = sbuf_min_left(s)) + num_aliases - 1) {
                int next = s->end_status;

                switch (next) {
                case 0:
                    s_process_read_buf(s);
                    continue;
                case EOFC:
                    at_eof_count++;
                    break;	/* with no data available */
                case INTC:
                case CALLC:
                    return
                        s_handle_read_exception(i_ctx_p, next, pp,
                                                NULL, 0, image_file_continue);
                default:
                    /* case ERRC: */
                    return_error(gs_error_ioerror);
                }
                break;		/* for EOFC */
            }
            /*
             * Note that in the EOF case, we can get here with no data
             * available.
             */
            if (avail >= min_left)
                avail = (avail - min_left) / num_aliases; /* may be 0 */
            if (avail < min_avail)
                min_avail = avail;
            plane_data[px].data = sbufptr(s);
            plane_data[px].size = avail;
        }

        /*
         * Now pass the available buffered data to the image processor.
         * Even if there is no available data, we must call
         * gs_image_next_planes one more time to finish processing any
         * retained data.
         */

        {
            int pi;
            uint used[GS_IMAGE_MAX_COMPONENTS];

            code = gs_image_next_planes(penum, plane_data, used);
            /* Now that used has been set, update the streams. */
            total_used = 0;
            for (pi = 0, pp = ETOP_SOURCE(esp, 0); pi < num_sources;
                 ++pi, pp -= 2 ) {
                (void)sbufskip(pp->value.pfile, used[pi]);
                total_used += used[pi];
            }
            if (code == gs_error_Remap_Color)
                return code;
        }
        if (at_eof_count >= num_sources || (at_eof_count && total_used == 0))
            code = 1;
        if (code) {
            int code1;

            esp = zimage_pop_estack(esp);
            code1 = image_cleanup(i_ctx_p);
            return (code < 0 ? code : code1 < 0 ? code1 : o_pop_estack);
        }
    }
}

/* Process data from an image with string data sources. */
/* This may still encounter a RemapColor callback. */
static int
image_string_continue(i_ctx_t *i_ctx_p)
{
    gs_image_enum *penum = r_ptr(esp, gs_image_enum);
    int num_sources = ETOP_NUM_SOURCES(esp)->value.intval;
    gs_const_string sources[GS_IMAGE_MAX_COMPONENTS];
    uint used[GS_IMAGE_MAX_COMPONENTS];

    /* Pass no data initially, to find out how much is retained. */
    memset(sources, 0, sizeof(sources[0]) * num_sources);
    for (;;) {
        int px;
        int code = gs_image_next_planes(penum, sources, used);

        if (code == gs_error_Remap_Color)
            return code;
    stop_now:
        if (code) {		/* Stop now. */
            esp -= NUM_PUSH(num_sources);
            image_cleanup(i_ctx_p);
            return (code < 0 ? code : o_pop_estack);
        }
        for (px = 0; px < num_sources; ++px)
            if (sources[px].size == 0) {
                const ref *psrc = ETOP_SOURCE(esp, px);
                uint size = r_size(psrc);

                if (size == 0) {	    /* empty source */
                    code = 1;
                    goto stop_now;
                }
                sources[px].data = psrc->value.bytes;
                sources[px].size = size;
            }
    }
}

/* Clean up after enumerating an image */
static int
image_cleanup(i_ctx_t *i_ctx_p)
{
    es_ptr ep_top = esp + NUM_PUSH(EBOT_NUM_SOURCES(esp)->value.intval);
    gs_image_enum *penum = r_ptr(ep_top, gs_image_enum);

    return gs_image_cleanup_and_free_enum(penum, igs);
}

/* ------ Initialization procedure ------ */

const op_def zimage_op_defs[] =
{
    {"1.image1", zimage1},
    {"1.imagemask1", zimagemask1},
                /* Internal operators */
    {"1%image_proc_continue", image_proc_continue},
    {"0%image_file_continue", image_file_continue},
    {"0%image_string_continue", image_string_continue},
    op_def_end(0)
};
