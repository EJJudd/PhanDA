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


/* Interpreter language level interface */

#ifndef ilevel_INCLUDED
#  define ilevel_INCLUDED

#include "imemory.h"

/* The current interpreter language level */
#define LANGUAGE_LEVEL (i_ctx_p->language_level)
#define LL2_ENABLED (LANGUAGE_LEVEL >= 2)
#define LL3_ENABLED (LANGUAGE_LEVEL >= 3)
#define level2_enabled LL2_ENABLED	/* backward compatibility */

#endif /* ilevel_INCLUDED */
