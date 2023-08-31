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



#ifndef idisp_INCLUDED
#  define idisp_INCLUDED

#include "imain.h"

#ifndef display_callback_DEFINED
# define display_callback_DEFINED
typedef struct display_callback_s display_callback;
#endif

/* Called from imain.c to reopen the device after initialisation if the.
 * device requires this. This gives it a chance to refetch any callbacks. */
int reopen_device_if_required(gs_main_instance *minst);

#endif /* idisp_INCLUDED */
