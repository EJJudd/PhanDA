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


/* Client interface to Pattern color */

#ifndef gspcolor_INCLUDED
#  define gspcolor_INCLUDED

#include "gsccolor.h"
#include "gsrefct.h"
#include "gsuid.h"
#include "gsgstate.h"
#include "gsmatrix.h"

/* ---------------- Types and structures ---------------- */

/*
 * We originally defined the gs_client_pattern structure before we
 * realized that we would have to accommodate multiple PatternTypes.
 * In version 5.68, we bit the bullet and made an incompatible change
 * to this structure so that multiple PatternTypes could be supported.
 * In order to make this work:
 *
 *	Clients creating instances of any Pattern template structure
 *	(gs_patternN_template_t) must call gs_patternN_init to
 *	initialize all the members, before filling in any of the
 *	members themselves.
 *
 * This is a non-backward-compatible requirement relative to previous
 * versions, but it was unavoidable.
 */

/*
 * Define the abstract pattern template (called "prototype pattern" in Red
 * Book).
 */

typedef struct gs_pattern_type_s gs_pattern_type_t;

#define gs_pattern_template_common\
  const gs_pattern_type_t *type;\
  int PatternType;		/* copied from the type structure */\
  gs_uid uid

typedef struct gs_pattern_template_s {
    gs_pattern_template_common;
} gs_pattern_template_t;

/* The descriptor is public for subclassing. */
extern_st(st_pattern_template);
#define public_st_pattern_template() /* in gspcolor.c */\
  gs_public_st_ptrs1(st_pattern_template, gs_pattern_template_t,\
    "gs_pattern_template_t", pattern_template_enum_ptrs,\
    pattern_template_reloc_ptrs, uid.xvalues)
#define st_pattern_template_max_ptrs 2

typedef void (*gs_pinst_free_proc_t) (gs_memory_t * mem, void *pinst);

/* Definition of Pattern instances. */
#define gs_pattern_instance_common\
    rc_header rc;\
    /* Following are set by makepattern */\
    const gs_pattern_type_t *type;  /* from template */\
    gs_gstate *saved;\
    void *client_data;		/* additional data for rendering */\
    gs_pinst_free_proc_t notify_free;\
    gs_id pattern_id
struct gs_pattern_instance_s {
    gs_pattern_instance_common;
};

/* The following is public for subclassing. */
extern_st(st_pattern_instance);
#define public_st_pattern_instance() /* in gspcolor.c */\
  gs_public_st_ptrs2(st_pattern_instance, gs_pattern_instance_t,\
    "gs_pattern_instance_t", pattern_instance_enum_ptrs,\
    pattern_instance_reloc_ptrs, saved, client_data)
#define st_pattern_instance_max_ptrs 1

/* ---------------- Procedures ---------------- */

/* Set a Pattern color or a Pattern color space. */
int gs_setpattern(gs_gstate *, const gs_client_color *);
int gs_setpatternspace(gs_gstate *);

/*
 * Construct a Pattern color of any PatternType.
 * The gs_memory_t argument for gs_make_pattern may be NULL, meaning use the
 * same allocator as for the gs_gstate argument.  Note that gs_make_pattern
 * uses rc_alloc_struct_1 to allocate pattern instances.
 */
int gs_make_pattern(gs_client_color *, const gs_pattern_template_t *,
                    const gs_matrix *, gs_gstate *, gs_memory_t *);
const gs_pattern_template_t *gs_get_pattern(const gs_client_color *);

/*
 * Adjust the reference count of a pattern. This is intended to support
 * applications (such as PCL) which maintain client colors outside of the
 * graphic state. Since the pattern instance structure is opaque to these
 * applications, they need some way to release or retain the instances as
 * needed.
 */
void gs_pattern_reference(gs_client_color * pcc, int delta);

#endif /* gspcolor_INCLUDED */
