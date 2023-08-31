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


/* $Id: zcie.h  $ */
/* Definitions for CIEBased color spaces */

#ifndef zcie_INCLUDED
#  define zcie_INCLUDED

#include "iref.h"

int cieaspace(i_ctx_t *i_ctx_p, ref *CIEdict, uint64_t dictkey);
int cieabcspace(i_ctx_t *i_ctx_p, ref *CIEDict, uint64_t dictkey);
int ciedefspace(i_ctx_t *i_ctx_p, ref *CIEDict, uint64_t dictkey);
int ciedefgspace(i_ctx_t *i_ctx_p, ref *CIEDict, uint64_t dictkey);

#endif /* zcie_INCLUDED */
