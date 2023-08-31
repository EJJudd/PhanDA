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


/* Token scanner state and interface */
/* Requires gsstruct.h, ostack.h, stream.h */

#ifndef iscan_INCLUDED
#  define iscan_INCLUDED

#include "sa85x.h"
#include "sstring.h"
#include "inamestr.h"
#include "iref.h"

/*
 * Define the state of the scanner.  Before first calling gs_scan_token,
 * the caller must initialize the state by calling scanner_state_init.
 * Most of the state is only used if scanning is suspended because of
 * an interrupt or a callout.
 *
 * Note that the scanner state includes a reference to the stream being
 * scanned.  We do this primarily for Adobe compatibility, so that we can
 * remove the source object from the operand stack after the initial checks.
 *
 * We expose the entire state definition to the caller so that
 * the state can normally be allocated on the stack.
 */
typedef struct scanner_state_s scanner_state;

/*
 * Define a structure for dynamically growable strings.
 * If is_dynamic is true, base/next/limit point to a string on the heap;
 * if is_dynamic is false, base/next/limit point either to the local buffer
 * or (only while control is inside gs_scan_token) into the source stream
 * buffer.
 */
#define max_comment_line 255    /* max size of an externally processable comment */
typedef struct dynamic_area_s {
    byte *base;
    byte *next;
    byte *limit;
    bool is_dynamic;
    byte buf[max_name_string];  /* initial buffer, enough for a valid string */
    gs_memory_t *memory;
} dynamic_area;

#define da_size(pda) ((uint)((pda)->limit - (pda)->base))
typedef dynamic_area *da_ptr;

/* Define state specific to binary tokens and binary object sequences. */
typedef struct scan_binary_state_s {
    int num_format;
    int (*cont)(i_ctx_t *, ref *, scanner_state *);
    ref bin_array;
    uint index;
    uint max_array_index;       /* largest legal index in objects */
    uint min_string_index;      /* smallest legal index in strings */
    uint top_size;
    uint size;
    int token_type;             /* binary token type for error reporting */
    ulong lsize;                /* b.o.s. size ibid. */
} scan_binary_state;

/* Define the scanner state. */
struct scanner_state_s {
    ref s_file;                 /* source file */
    uint s_pstack;              /* stack depth when starting current */
                                /* procedure, after pushing old pstack */
    uint s_pdepth;              /* pstack for very first { encountered, */
                                /* for error recovery */
    int s_options;
    enum {
        scanning_none,
        scanning_binary,
        scanning_comment,
        scanning_name,
        scanning_string
    } s_scan_type;
    dynamic_area s_da;
    union sss_ {                /* scanning sub-state */
        scan_binary_state binary;       /* binary */
        struct sns_ {           /* name */
            int s_name_type;    /* number of /'s preceding a name */
            bool s_try_number;  /* true if should try scanning name */
            /* as number */
        } s_name;
        stream_state st;        /* string */
        stream_A85D_state a85d; /* string */
        stream_AXD_state axd;   /* string */
        stream_PSSD_state pssd; /* string */
    } s_ss;
    /* The following are used only to return information for errors. */
    struct se_ {                /* scanner error */
        ref object;             /* normally t__invalid */
        bool is_name;           /* true if 'string' is name, false if string */
#define SCANNER_MAX_ERROR_STRING 120 /* adhoc, for Adobe-compatible messages */
        char string[SCANNER_MAX_ERROR_STRING+1]; /* normally empty */
    } s_error;
#define SCAN_INIT_ERROR(pstate)\
  (make_t(&(pstate)->s_error.object, t__invalid),\
   (pstate)->s_error.is_name = false,\
   (pstate)->s_error.string[0] = 0)
};

/* Remember the allocator for proper freeing. */
typedef struct scanner_state_dynamic_s {
    scanner_state state;
    gs_memory_t *mem;
} scanner_state_dynamic;

/* The type descriptor is public only for checking. */
extern_st(st_scanner_state_dynamic);
#define public_st_scanner_state_dynamic()       /* in iscan.c */\
  gs_public_st_complex_only(st_scanner_state_dynamic, scanner_state_dynamic, "scanner state",\
    scanner_clear_marks, scanner_enum_ptrs, scanner_reloc_ptrs, 0)

/* Initialize a scanner with a given set of options. */
#define SCAN_FROM_STRING 1      /* true if string is source of data */
                                /* (for Level 1 `\' handling) */
#define SCAN_CHECK_ONLY 2       /* true if just checking for syntax errors */
                                /* and complete statements (no value) */
#define SCAN_PROCESS_COMMENTS 4 /* return scan_Comment for comments */
                                /* (all comments or only non-DSC) */
#define SCAN_PROCESS_DSC_COMMENTS 8  /* return scan_DSC_Comment */
#define SCAN_PDF_RULES 16       /* PDF scanning rules (for names) */
#define SCAN_PDF_INV_NUM 32     /* Adobe ignores invalid numbers */
                                /* This is for compatibility with Adobe */
                                /* Acrobat Reader                       */
#define SCAN_PDF_UNSIGNED 64    /* Scan 2147483648..4294967295 as unsigned numbers */
                                /* This is needed in some contexts for */
                                /* compatibility with Adobe */
#define SCAN_CPSI_MODE 128      /* Flag to indicate CPSI compatible integer parsing */

void gs_scanner_init_options(scanner_state *sstate, const ref *fop,
                             int options);
#define gs_scanner_init(sstate, fop)\
  gs_scanner_init_options(sstate, fop, 0)
void gs_scanner_init_stream_options(scanner_state *sstate, stream *s,
                                    int options);
#define gs_scanner_init_stream(sstate, s)\
  gs_scanner_init_stream_options(sstate, s, 0)

/*
 * Read a token from a stream.  As usual, 0 is a normal return,
 * <0 is an error.  There are also some special return codes:
 */
#define scan_BOS 1              /* binary object sequence */
#define scan_EOF 2              /* end of stream */
#define scan_Refill 3           /* get more input data, then call again */
#define scan_Comment 4          /* comment, non-DSC if processing DSC */
#define scan_DSC_Comment 5      /* DSC comment */
int gs_scan_token(i_ctx_t *i_ctx_p, ref * pref, scanner_state * pstate);

/*
 * Read a token from a string.  Return like gs_scan_token, but also
 * update the string to move past the token (if no error).
 */
int gs_scan_string_token_options(i_ctx_t *i_ctx_p, ref * pstr, ref * pref,
                                 int options);
#define gs_scan_string_token(i_ctx_p, pstr, pref)\
  gs_scan_string_token_options(i_ctx_p, pstr, pref, 0)

/*
 * Return the "error object" to be stored in $error.command instead of
 * --token--, if any, or -1 if no special error object is required.
 */
int gs_scanner_error_object(i_ctx_t *i_ctx_p, const scanner_state *pstate,
                            ref *pseo);

/*
 * Handle a scan_Refill return from gs_scan_token.
 * This may return o_push_estack, 0 (meaning just call gs_scan_token
 * again), or an error code.
 */
int gs_scan_handle_refill(i_ctx_t *i_ctx_p, scanner_state * pstate,
                          bool save, op_proc_t cont);

/*
 * Define the procedure "hook" for parsing DSC comments.  If not NULL,
 * this procedure is called for every DSC comment seen by the scanner.
 */
extern int (*gs_scan_dsc_proc) (const byte *, uint);

/*
 * Define the procedure "hook" for parsing general comments.  If not NULL,
 * this procedure is called for every comment seen by the scanner.
 * If both gs_scan_dsc_proc and gs_scan_comment_proc are set,
 * gs_scan_comment_proc is called only for non-DSC comments.
 */
extern int (*gs_scan_comment_proc) (const byte *, uint);

#endif /* iscan_INCLUDED */
