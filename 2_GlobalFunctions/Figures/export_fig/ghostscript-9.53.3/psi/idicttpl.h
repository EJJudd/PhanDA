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


/* A template for packed dictionary search method */

#if defined(found) && defined(deleted) && defined(missing)

/*
 * Define template for searching a packed dictionary.
 *
 * Free variables:
 *      ref_packed kpack - holds the packed key.
 *      uint hash - holds the hash of the name.
 *      dict *pdict - points to the dictionary.
 *      uint size - holds npairs(pdict).
 *
 * Template parameters are :
 *	found   - the found key action.
 *	deleted - the deleted key action.
 *	missing - the missed key action.
 *
 * Note that the template is *not* enclosed in {}, so that we can access
 * the values of kbot and kp after leaving the template.
 */

    const ref_packed *kbot = pdict->keys.value.packed;
    const int start = dict_hash_mod(hash, size) + 1;
    register const ref_packed *kp = kbot + start;
    int wrap = 0;

    again:
    for (; ; kp-- ) {
        if_debug2('D', "[D]probe "PRI_INTPTR": 0x%x\n", (intptr_t)kp, *kp);
        if ( *kp == kpack ) {
            found;
        } else if ( !r_packed_is_name(kp) ) {
            /* Empty, deleted, or wraparound. Figure out which. */
            if ( *kp == packed_key_empty )
                missing;
            if ( kp == kbot ) {
                if (wrap)
                    break;
                else {
                    wrap++;
                    kp += size; /* wrap */
                    goto again; /* skip "kp--". */
                }
            } else {
                deleted;
            }
        }
   }
#else
int dummy;
#endif
