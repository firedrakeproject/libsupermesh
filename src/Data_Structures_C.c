#include "libsupermesh_configuration.h"

#ifdef ENABLE_JUDY

#include "Judy.h"
#include "stdio.h"
#include "assert.h"

/* To understand these, read
   http://judy.sourceforge.net/doc/Judy1_3x.htm and
   http://judy.sourceforge.net/doc/JudyL_3x.htm */

void libsupermesh_integer_set_create(Pvoid_t* i)
{
  *i = (Pvoid_t) NULL;
}

void libsupermesh_integer_set_delete(Pvoid_t* i)
{
  Word_t mem_freed;
  Pvoid_t ptr = (Pvoid_t) *i;
  J1FA(mem_freed, ptr);
  *i = ptr;
}

void libsupermesh_integer_set_insert(Pvoid_t* i, int* v, int* changed)
{
  Word_t index = *v;
  Pvoid_t ptr = (Pvoid_t) *i;
  J1S(*changed, ptr, index);
  *i = ptr;
}

void libsupermesh_integer_set_length(Pvoid_t* i, int* l)
{
  Word_t len;
  Pvoid_t ptr = (Pvoid_t) *i;
  J1C(len, ptr, 0, -1);
  *l = len;
  *i = ptr;
}

void libsupermesh_integer_set_fetch(Pvoid_t* i, int* idx, int* v)
{
  Word_t index = *idx, value;
  int worked;
  Pvoid_t ptr = (Pvoid_t) *i;
  J1BC(worked, ptr, index, value); 
  assert(worked == 1);
  *v = value;
  *i = ptr;
}

void libsupermesh_integer_set_remove(Pvoid_t* i, int* v, int* status)
{
  Word_t index = *v;
  int stat;
  Pvoid_t ptr = (Pvoid_t) *i;
  J1U(stat, ptr, index); 
  *status = stat;
  *i = ptr;
}

void libsupermesh_integer_set_has_value(Pvoid_t* i, int* v, int* present)
{
  Word_t value = *v;
  int wpresent;
  Pvoid_t ptr = (Pvoid_t) *i;
  J1T(wpresent, ptr, value); 
  *present = wpresent;
}

void libsupermesh_integer_hash_table_create(Pvoid_t* i)
{
  assert(sizeof(int*) == sizeof(Pvoid_t));
  *i = (Pvoid_t) NULL;
}

void libsupermesh_integer_hash_table_delete(Pvoid_t* i)
{
  Word_t mem_freed;
  Pvoid_t ptr = (Pvoid_t) *i;
  JLFA(mem_freed, ptr);
  *i = ptr;
}

void libsupermesh_integer_hash_table_insert(Pvoid_t* i, int* k, int* v)
{
  Word_t key = *k;
  PWord_t pvalue;
  Pvoid_t ptr = (Pvoid_t) *i;
  JLI(pvalue, ptr, key);
  *pvalue = *v;
  *i = ptr;
}

void libsupermesh_integer_hash_table_length(Pvoid_t* i, int* l)
{
  Word_t len;
  Pvoid_t ptr = (Pvoid_t) *i;
  JLC(len, ptr, 0, -1);
  *l = len;
  *i = ptr;
}

void libsupermesh_integer_hash_table_fetch(Pvoid_t* i, int* k, int* v)
{
  Word_t key = *k, value;
  PWord_t pvalue = &value;
  Pvoid_t ptr = (Pvoid_t) *i;
  JLG(pvalue, ptr, key); 
  if (pvalue == NULL)
  {
    fprintf(stderr, "Error: hash table has no key %d\n", *k);
    assert(pvalue != NULL);
  }
  *v = *pvalue;
  *i = ptr;
}

void libsupermesh_integer_hash_table_remove(Pvoid_t* i, int* k, int* status)
{
  Word_t key = *k;
  int stat;
  Pvoid_t ptr = (Pvoid_t) *i;
  JLD(stat, ptr, key); 
  *status = stat;
  *i = ptr;
}

void libsupermesh_integer_hash_table_has_key(Pvoid_t* i, int* k, int* present)
{
  Word_t key = *k, value;
  PWord_t pvalue = &value;
  Pvoid_t ptr = (Pvoid_t) *i;
  JLG(pvalue, ptr, key); 
  *present = (pvalue != NULL);
  *i = ptr;
}

void libsupermesh_integer_hash_table_fetch_pair(Pvoid_t* i, int* idx, int* k, int* v)
{
  Word_t nth = *idx; /* what Judy calls nth is what I am calling idx */
  Word_t index; /* what Judy calls index is what I am calling k */
  PWord_t pvalue;
  Pvoid_t ptr = (Pvoid_t) *i;
  JLBC(pvalue, ptr, nth, index); 
  assert(pvalue != NULL);
  *k = index;
  *v = *pvalue;
  *i = ptr;
}

#endif