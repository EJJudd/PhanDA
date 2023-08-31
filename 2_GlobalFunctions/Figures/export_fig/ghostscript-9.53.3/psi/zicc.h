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


/* Definitions for setcolorspace */

#ifndef zicc_INCLUDED
#  define zicc_INCLUDED

#include "iref.h"

int seticc(i_ctx_t * i_ctx_p, int ncomps, ref *ICCdict, float *range_buff);
int seticc_lab(i_ctx_t * i_ctx_p, float *white, float *black, float *range_buff);
int seticc_cal(i_ctx_t * i_ctx_p, float *white, float *black, float *gamma,
               float *matrix, int num_colorants,ulong dictkey);
#endif /* zicc_INCLUDED */
