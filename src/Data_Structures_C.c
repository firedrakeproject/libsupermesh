/*
  Copyright (C) 2016 The University of Edinburgh

  The file is part of libsupermesh
    
  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation;
  version 2.1 of the License.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/

/*
 * The following code is derived from femtools/Data_structures_C.c in Fluidity
 * git revision 4e6c1d2b022df3a519cdec120fad28e60d1b08d9 (dated 2015-02-25)
 */
 
/*
 * Fluidity copyright information (note that AUTHORS mentioned in the following
 * has been renamed to fluidity_AUTHORS):
 */

/*
  Copyright (C) 2006 Imperial College London and others.
  
  Please see the AUTHORS file in the main source directory for a full list
  of copyright holders.

  Prof. C Pain
  Applied Modelling and Computation Group
  Department of Earth Science and Engineering
  Imperial College London

  amcgsoftware@imperial.ac.uk
  
  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation,
  version 2.1 of the License.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
  USA
*/

#include "Data_Structures_C.h"

#ifdef LIBSUPERMESH_ENABLE_JUDY

#include <Judy.h>
#include <stdio.h>

#include "libsupermesh_debug_C.h"

/* To understand these, read
   http://judy.sourceforge.net/doc/Judy1_3x.htm and
   http://judy.sourceforge.net/doc/JudyL_3x.htm */

void libsupermesh_integer_set_create(Pvoid_t *i) {
  *i = (Pvoid_t) NULL;
}

void libsupermesh_integer_set_delete(Pvoid_t *i) {
  Word_t mem_freed;
  Pvoid_t ptr = (Pvoid_t) *i;
  J1FA(mem_freed, ptr);
  *i = ptr;
}

void libsupermesh_integer_set_insert(Pvoid_t *i, const int *v, int *changed) {
  Word_t index = *v;
  Pvoid_t ptr = (Pvoid_t) *i;
  J1S(*changed, ptr, index);
  *i = ptr;
}

void libsupermesh_integer_set_length(Pvoid_t *i, int *l) {
  Word_t len;
  Pvoid_t ptr = (Pvoid_t) *i;
  J1C(len, ptr, 0, -1);
  *l = len;
  *i = ptr;
}

void libsupermesh_integer_set_fetch(Pvoid_t *i, const int *idx, int *v) {
  Word_t index = *idx, value;
  int worked;
  Pvoid_t ptr = (Pvoid_t) *i;
  J1BC(worked, ptr, index, value);
  if(!worked) {
    fprintf(stderr, "Failed to fetch integer set element with index '%i'\n", *idx);
    libsupermesh_abort("Failed to fetch integer set element");
  }
  *v = value;
  *i = ptr;
}

void libsupermesh_integer_set_remove(Pvoid_t *i, const int *v, int *status) {
  Word_t index = *v;
  int stat;
  Pvoid_t ptr = (Pvoid_t) *i;
  J1U(stat, ptr, index); 
  *status = stat;
  *i = ptr;
}

void libsupermesh_integer_set_has_value(Pvoid_t *i, const int *v, int *present) {
  Word_t value = *v;
  int wpresent;
  Pvoid_t ptr = (Pvoid_t) *i;
  J1T(wpresent, ptr, value); 
  *present = wpresent;
}

void libsupermesh_integer_hash_table_create(Pvoid_t *i) {
  *i = (Pvoid_t) NULL;
}

void libsupermesh_integer_hash_table_delete(Pvoid_t *i) {
  Word_t mem_freed;
  Pvoid_t ptr = (Pvoid_t) *i;
  JLFA(mem_freed, ptr);
  *i = ptr;
}

void libsupermesh_integer_hash_table_insert(Pvoid_t *i, const int *k, const int *v) {
  Word_t key = *k;
  PWord_t pvalue;
  Pvoid_t ptr = (Pvoid_t) *i;
  JLI(pvalue, ptr, key);
  *pvalue = *v;
  *i = ptr;
}

void libsupermesh_integer_hash_table_length(Pvoid_t *i, int *l) {
  Word_t len;
  Pvoid_t ptr = (Pvoid_t) *i;
  JLC(len, ptr, 0, -1);
  *l = len;
  *i = ptr;
}

void libsupermesh_integer_hash_table_fetch(Pvoid_t *i, const int *k, int *v) {
  Word_t key = *k, value;
  PWord_t pvalue = &value;
  Pvoid_t ptr = (Pvoid_t) *i;
  JLG(pvalue, ptr, key); 
  if(!pvalue) {
    fprintf(stderr, "Failed to fetch integer hash table element with key '%i'\n", *k);
    libsupermesh_abort("Failed to fetch integer hash table element");
  }
  *v = *pvalue;
  *i = ptr;
}

void libsupermesh_integer_hash_table_remove(Pvoid_t *i, const int *k, int *status) {
  Word_t key = *k;
  int stat;
  Pvoid_t ptr = (Pvoid_t) *i;
  JLD(stat, ptr, key); 
  *status = stat;
  *i = ptr;
}

void libsupermesh_integer_hash_table_has_key(Pvoid_t *i, const int *k, int *present) {
  Word_t key = *k, value;
  PWord_t pvalue = &value;
  Pvoid_t ptr = (Pvoid_t) *i;
  JLG(pvalue, ptr, key); 
  *present = (pvalue != NULL);
  *i = ptr;
}

void libsupermesh_integer_hash_table_fetch_pair(Pvoid_t *i, const int *idx, int *k, int *v) {
  Word_t nth = *idx; /* what Judy calls nth is what I am calling idx */
  Word_t index; /* what Judy calls index is what I am calling k */
  PWord_t pvalue;
  Pvoid_t ptr = (Pvoid_t) *i;
  JLBC(pvalue, ptr, nth, index); 
  if(!pvalue) {
    fprintf(stderr, "Failed to fetch integer hash table element with index '%i'\n", *idx);
    libsupermesh_abort("Failed to fetch integer hash table element");
  }
  *k = index;
  *v = *pvalue;
  *i = ptr;
}

#endif
