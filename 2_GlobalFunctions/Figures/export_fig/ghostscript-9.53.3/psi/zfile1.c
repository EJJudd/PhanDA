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


/* Special file operators */

#include "memory_.h"
#include "string_.h"
#include "ghost.h"
#include "gp.h"
#include "ierrors.h"
#include "oper.h"
#include "ialloc.h"
#include "opdef.h"
#include "opcheck.h"
#include "store.h"
#include "gpmisc.h"

/* <string> <string> <bool> .file_name_combine <string> true */
/* <string> <string> <bool> .file_name_combine <string> <string> false */
static int
zfile_name_combine(i_ctx_t *i_ctx_p)
{
    uint plen, flen, blen, blen0;
    const byte *prefix, *fname;
    byte *buffer;
    os_ptr op = osp;
    bool no_sibling;

    check_type(op[ 0], t_boolean);
    check_type(op[-1], t_string);
    check_type(op[-2], t_string);
    plen = r_size(op - 2);
    flen = r_size(op - 1);
    blen = blen0 = plen + flen + 2; /* Inserts separator and ending zero byte. */
    buffer = ialloc_string(blen, "zfile_name_combine");
    if (buffer == 0)
        return_error(gs_error_VMerror);
    prefix = op[-2].value.const_bytes;
    fname =  op[-1].value.const_bytes;
    no_sibling = op[0].value.boolval;
    if (gp_file_name_combine((const char *)prefix, plen,
                             (const char *)fname, flen, no_sibling,
                             (char *)buffer, &blen) != gp_combine_success) {
        make_bool(op, false);
    } else {
        buffer = iresize_string(buffer, blen0, blen, "zfile_name_combine");
        if (buffer == 0)
            return_error(gs_error_VMerror);
        make_string(op - 2, a_all | icurrent_space, blen, buffer);
        make_bool(op - 1, true);
        pop(1);
    }
    return 0;
}

/* This is compiled conditionally to let PS library to know
 * whether it works with the new gp_combine_file_name.
 */

/* <string> .file_name_is_absolute <bool> */
static int
zfile_name_is_absolute(i_ctx_t *i_ctx_p)
{   os_ptr op = osp;

    check_type(op[0], t_string);
    make_bool(op, gp_file_name_is_absolute((const char *)op->value.const_bytes,
                                        r_size(op)));
    return 0;
}

static int
push_string(i_ctx_t *i_ctx_p, const char *v)
{   os_ptr op = osp;
    int len = strlen(v);

    push(1);
    make_const_string(op, avm_foreign | a_readonly,
                      len, (const byte *)v);
    return 0;
}

/* - .file_name_separator <string> */
static int
zfile_name_separator(i_ctx_t *i_ctx_p)
{   return push_string(i_ctx_p, gp_file_name_separator());
}

/* - .file_name_directory_separator <string> */
static int
zfile_name_directory_separator(i_ctx_t *i_ctx_p)
{   return push_string(i_ctx_p, gp_file_name_directory_separator());
}

/* - .file_name_current <string> */
static int
zfile_name_current(i_ctx_t *i_ctx_p)
{   return push_string(i_ctx_p, gp_file_name_current());
}

/* - .file_name_parent <string> */
static int
zfile_name_parent(i_ctx_t *i_ctx_p)
{   return push_string(i_ctx_p, gp_file_name_parent());
}

const op_def zfile1_op_defs[] =
{
    {"3.file_name_combine", zfile_name_combine},
    {"1.file_name_is_absolute", zfile_name_is_absolute},
    {"0.file_name_separator", zfile_name_separator},
    {"0.file_name_directory_separator", zfile_name_directory_separator},
    {"0.file_name_current", zfile_name_current},
    {"0.file_name_parent", zfile_name_parent},
    op_def_end(0)
};
