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


/* Private line parameter definitions */

#ifndef gxline_INCLUDED
#  define gxline_INCLUDED

#include "math_.h"
#include "gslparam.h"
#include "gsmatrix.h"

/* Line parameter structures */
/* gx_dash_params are never instantiated by themselves. */
typedef struct gx_dash_params_s {
    float *pattern;
    uint pattern_size;
    float offset;
    bool adapt;
    /* The rest of the parameters are computed from the above */
    float pattern_length;	/* total of all pattern elements */
    bool init_ink_on;		/* true if ink is initially on */
    int init_index;		/* initial index in pattern */
    float init_dist_left;
} gx_dash_params;

#define gx_dash_params_initial\
  NULL, 0, 0.0, 0/*false*/, 0.0, 1/*true*/, 0, 0.0
typedef struct gx_line_params_s {
    float half_width;		/* one-half line width */
    gs_line_cap start_cap;      /* Cap to use on start of line */
    gs_line_cap end_cap;        /* Cap to use on end of line */
    gs_line_cap dash_cap;       /* Cap to use on start/end of dash segment */
    gs_line_join join;
    int curve_join;		/* <0 means use join between segments of */
                                /* flattened curves, >=0 means use this join */
    float miter_limit;
    float miter_check;		/* computed from miter limit, */
                                /* see gx_set_miter_limit and gs_stroke */
    float dot_length;
    bool dot_length_absolute;	/* if true, dot_length is 1/72" units */
    gs_matrix dot_orientation;	/* dot length is aligned with (1,0); */
                                /* must be xxyy or xyyx */
    gx_dash_params dash;
} gx_line_params;

#define gx_set_line_width(plp, wid)\
  ((plp)->half_width = fabs(wid) / 2)
#define gx_current_line_width(plp)\
  ((plp)->half_width * 2)
int gx_set_miter_limit(gx_line_params *, double);

#define gx_current_miter_limit(plp)\
  ((plp)->miter_limit)
int gx_set_dash(gx_dash_params *, const float *, uint, double, gs_memory_t *);

#define gx_set_dash_adapt(pdp, adpt) ((pdp)->adapt = (adpt))
int gx_set_dot_length(gx_line_params *, double, bool);

/* See gsline.c for the computation of miter_check. */
#define gx_line_params_initial\
 0.0, gs_cap_butt, gs_cap_butt, gs_cap_butt, gs_join_miter, -1,\
 10.0, (float)0.20305866, 0.0, 0/*false*/,\
 { identity_matrix_body }, { gx_dash_params_initial }

#endif /* gxline_INCLUDED */
