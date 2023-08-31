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


/* Reference counting definitions */

#ifndef gsrefct_INCLUDED
#  define gsrefct_INCLUDED

#include "memento.h"
#include "std.h"

/*
 * A reference-counted object must include the following header:
 *      rc_header rc;
 * The header need not be the first element of the object.
 *
 * Reference-counted objects have a freeing procedure that gets called when
 * the reference count reaches zero.  In retrospect, we probably should have
 * used finalization for this, but it's too difficult to change now.
 * Because of the interaction between these two features, the freeing
 * procedure for reference-counted objects that do use finalization must
 * free the object itself first, before decrementing the reference counts
 * of referenced objects (which of course requires saving pointers to those
 * objects before freeing the containing object).
 *
 * To add a reference to an object, copy its pointer and call:
 *      rc_increment(pobj);
 *
 * We provide two decrement macros for reference-counted objects:
 *      rc_decrement(pobj);
 *      rc_decrement_only(pobj);
 * Both decrement the reference count, and free the object if that
 * was the last reference. The difference is that rc_decrement()
 * also assigns zero to its argument, while rc_decrement_only()
 * does not. The first must be used if the argument could be traceable
 * by the allocator to avoid examining stale memory when a garbage
 * collector runs. The second must be used if the pointer argument
 * is read-only.
 */
typedef struct rc_header_s rc_header;
struct rc_header_s {
    long ref_count;
    gs_memory_t *memory;
#define rc_free_proc(proc)\
  void proc(gs_memory_t *, void *, client_name_t)
    rc_free_proc((*free));
};

#ifdef DEBUG
void rc_trace_init_free(const void *vp, const rc_header *prc);
void rc_trace_free_struct(const void *vp, const rc_header *prc,
                          client_name_t cname);
void rc_trace_increment(const void *vp, const rc_header *prc);
void rc_trace_adjust(const void *vp, const rc_header *prc, int delta, const char *cname);
#define IF_RC_DEBUG(call) BEGIN if (gs_debug_c('^')) { dlputs(""); call; } END
#else
#define IF_RC_DEBUG(call) DO_NOTHING
#endif

/* ------ Allocate/free ------ */

rc_free_proc(rc_free_struct_only);
/*
 * rc_init[_free] really should be used only to initialize
 * stack-allocated structures; with heap-allocated structures, it's
 * better to use a finalize method so that the garbage collector can
 * clean them up if the refcount fails to reach zero.
 */
#define rc_init_free(vp, mem, rcinit, proc)\
  BEGIN\
    (vp)->rc.ref_count = rcinit;\
    (vp)->rc.memory = mem;\
    (vp)->rc.free = proc;\
    IF_RC_DEBUG(rc_trace_init_free(vp, &(vp)->rc));\
  END
#define rc_init(vp, mem, rcinit)\
  rc_init_free(vp, mem, rcinit, rc_free_struct_only)

#define rc_alloc_struct_n(vp, typ, pstyp, mem, errstat, cname, rcinit)\
  BEGIN\
    if ( ((vp) = gs_alloc_struct(mem, typ, pstyp, cname)) == 0 ) {\
      errstat;\
    } else {\
      rc_init(vp, mem, rcinit);\
    }\
  END
#define rc_alloc_struct_0(vp, typ, pstype, mem, errstat, cname)\
  rc_alloc_struct_n(vp, typ, pstype, mem, errstat, cname, 0)
#define rc_alloc_struct_1(vp, typ, pstype, mem, errstat, cname)\
  rc_alloc_struct_n(vp, typ, pstype, mem, errstat, cname, 1)

#define rc_free_struct(vp, cname)\
  BEGIN\
    IF_RC_DEBUG(rc_trace_free_struct(vp, &(vp)->rc, cname));\
    (vp)->rc.free((vp)->rc.memory, (void *)(vp), cname);\
  END

/* ------ Reference counting ------ */

/* Increment a reference count. */
#define RC_DO_INCREMENT(vp)\
  BEGIN\
    (vp)->rc.ref_count++;\
    IF_RC_DEBUG(rc_trace_increment(vp, &(vp)->rc));\
  END
#define rc_increment(vp)\
  BEGIN\
    if (vp) RC_DO_INCREMENT(vp);\
  END

/* Increment a reference count, allocating the structure if necessary. */
#define rc_allocate_struct(vp, typ, pstype, mem, errstat, cname)\
  BEGIN\
    if (vp)\
      RC_DO_INCREMENT(vp);\
    else\
      rc_alloc_struct_1(vp, typ, pstype, mem, errstat, cname);\
  END

/* Guarantee that a structure is allocated and is not shared. */
#define RC_DO_ADJUST(vp, delta, cname)\
  BEGIN\
    IF_RC_DEBUG(rc_trace_adjust(vp, &(vp)->rc, delta, cname));\
    (vp)->rc.ref_count += (delta);\
  END
#define rc_unshare_struct(vp, typ, pstype, mem, errstat, cname)\
  BEGIN\
    if ( (vp) == 0 || (vp)->rc.ref_count > 1 || (vp)->rc.memory != (mem) ) {\
      typ *new;\
      if ( vp ) RC_DO_ADJUST(vp, -1, cname);\
      rc_alloc_struct_1(new, typ, pstype, mem, errstat, cname);\
      (vp) = new;\
    }\
  END

/* Adjust a reference count either up or down. */
#ifdef DEBUG
#  define rc_check_(vp)\
     BEGIN\
       if (gs_debug_c('?') && (vp)->rc.ref_count < 0)\
         lprintf2(PRI_INTPTR" has ref_count of %ld!\n", (intptr_t)(vp),\
                  (vp)->rc.ref_count);\
     END
#else
#  define rc_check_(vp) DO_NOTHING
#endif
#define rc_adjust_(vp, delta, cname, body)\
  BEGIN\
    if (vp) {\
      RC_DO_ADJUST(vp, delta, cname);\
      if (!(vp)->rc.ref_count) {\
        rc_free_struct(vp, cname);\
        body;\
      } else\
        rc_check_(vp);\
    }\
  END
#define rc_adjust(vp, delta, cname)\
  rc_adjust_(vp, delta, cname, (vp) = 0)
#define rc_adjust_only(vp, delta, cname)\
  rc_adjust_(vp, delta, cname, DO_NOTHING)
#define rc_adjust_const(vp, delta, cname)\
  rc_adjust_only(vp, delta, cname)
#define rc_decrement(vp, cname)\
  rc_adjust(vp, -1, cname)
#define rc_decrement_only(vp, cname)\
  rc_adjust_only(vp, -1, cname)

/*
 * Assign a pointer, adjusting reference counts.  vpfrom might be a local
 * variable with a copy of the last reference to the object, and freeing
 * vpto might decrement the object's reference count and cause it to be
 * freed (incorrectly); for that reason, we do the increment first.
 */
#define rc_assign(vpto, vpfrom, cname)\
  BEGIN\
    if ((vpto) != (vpfrom)) {\
      rc_increment(vpfrom);\
      rc_decrement_only(vpto, cname);\
      (vpto) = (vpfrom);\
    }\
  END
/*
 * Adjust reference counts for assigning a pointer,
 * but don't do the assignment.  We use this before assigning
 * an entire structure containing reference-counted pointers.
 */
#define rc_pre_assign(vpto, vpfrom, cname)\
  BEGIN\
    if ((vpto) != (vpfrom)) {\
      rc_increment(vpfrom);\
      rc_decrement_only(vpto, cname);\
    }\
  END

#endif /* gsrefct_INCLUDED */
