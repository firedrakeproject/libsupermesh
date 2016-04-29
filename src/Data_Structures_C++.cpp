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
 
// Fluidity copyright information (note that AUTHORS mentioned in the following
// has been renamed to fluidity_AUTHORS):

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

#include "libsupermesh_configuration.h"

#ifndef LIBSUPERMESH_ENABLE_JUDY

#include <cstdlib>
#include <iostream>
#include <map>
#include <set>

#include "libsupermesh_debug_C.h"

using namespace std;

namespace libsupermesh {

class integer_set {
  public:
    inline integer_set(void) {
      this->v_array = NULL;
      
      return;
    }
    
    inline ~integer_set(void) {
      if(this->v_array) {
        free(this->v_array);
      }
      
      return;
    }
    
    inline void insert(const int &v, int &changed) {
      if(this->v_array) {
        free(this->v_array);
        this->v_array = NULL;
      }
      pair<set<int>::const_iterator, bool> position = this->v_set.insert(v);
      changed = position.second ? 1 : 0;
      
      return;
    }
    
    inline void length(int &l) const {
      l = this->v_set.size();
      
      return;
    };
    
    inline void fetch(const int &idx, int &v) {
      if(!this->v_array) {
        this->v_array = (int*)malloc(this->v_set.size() * sizeof(int));
        if(!this->v_array) {
          libsupermesh_abort("malloc failure");
        }
        set<int>::size_type i = 0;
        for(set<int>::const_iterator iter = this->v_set.begin();iter != this->v_set.end();iter++) {
          this->v_array[i] = *iter;
          i++;
        }
      }
      v = this->v_array[idx - 1];
      
      return;
    };
    
    inline void remove(const int &v, int &status) {
      if(this->v_array) {
        free(this->v_array);
        this->v_array = NULL;
      }
      set<int>::size_type count = this->v_set.erase(v);
      status = (count > 0) ? 1 : 0;
      
      return;
    }
    
    inline void has_value(const int &v, int &present) const {
      set<int>::const_iterator iter = this->v_set.find(v);
      present = (iter == this->v_set.end()) ? 0 : 1;
      
      return;
    }
    
  private:
    set<int> v_set;
    int *v_array;
};

class integer_hash_table {
  public:
    inline integer_hash_table(void) {
      this->k_array = NULL;
      
      return;
    }
    
    inline ~integer_hash_table(void) {
      if(this->k_array) {
        free(this->k_array);
      }
    }
    
    inline void insert(const int &k, const int &v) {
      if(this->k_array) {
        free(this->k_array);
        this->k_array = NULL;
      }
      this->v_map[k] = v;
      
      return;
    }
    
    inline void length(int &l) const {
      l = this->v_map.size();
      
      return;
    }
    
    inline void fetch(const int &k, int &v) {
      v = this->v_map[k];
      
      return;
    }
    
    inline void remove(const int &k, int &status) {
      if(this->k_array) {
        free(this->k_array);
        this->k_array = NULL;
      }
      map<int, int>::size_type count = this->v_map.erase(k);
      status = (count > 0) ? 1 : 0;
      
      return;
    }
    
    inline void has_key(const int &k, int &present) const {
      map<int, int>::const_iterator iter = this->v_map.find(k);
      present = (iter == this->v_map.end()) ? 0 : 1;
      
      return;
    }
    
    inline void fetch_pair(const int &idx, int &k, int &v) {
      if(!this->k_array) {
        this->k_array = (int*)malloc(this->v_map.size() * sizeof(int));
        if(!this->k_array) {
          libsupermesh_abort("malloc failure");
        }
        map<int, int>::size_type i = 0;
        for(map<int, int>::const_iterator iter = this->v_map.begin();iter != this->v_map.end();iter++) {
          this->k_array[i] = iter->first;
          i++;
        }
      }
      k = this->k_array[idx - 1];
      v = this->v_map[k];
      
      return;
    }
    
  private:
    map<int, int> v_map;
    int *k_array;
};

}

extern "C" {
  void libsupermesh_integer_set_create(void **i) {
    (*i) = (void*)(new libsupermesh::integer_set());
    if(!((libsupermesh::integer_set*)(*i))) {
      libsupermesh_abort("new failure");
    }
    
    return;
  }
  
  void libsupermesh_integer_set_delete(void **i) {
    delete ((libsupermesh::integer_set*)(*i));
    
    return;
  }
  
  void libsupermesh_integer_set_insert(void **i, const int *v, int *changed) {
    ((libsupermesh::integer_set*)(*i))->insert(*v, *changed);
    
    return;
  }
  
  void libsupermesh_integer_set_length(const void **i, int *l) {
    ((libsupermesh::integer_set*)(*i))->length(*l);
    
    return;
  }
  
  void libsupermesh_integer_set_fetch(void **i, const int *idx, int *v) {
    ((libsupermesh::integer_set*)(*i))->fetch(*idx, *v);
    
    return;
  }
  
  void libsupermesh_integer_set_remove(void **i, const int *v, int *status) {
    ((libsupermesh::integer_set*)(*i))->remove(*v, *status);
    
    return;
  }

  void libsupermesh_integer_set_has_value(const void **i, const int *v, int *present) {
    ((libsupermesh::integer_set*)(*i))->has_value(*v, *present);
    
    return;
  }
  
  void libsupermesh_integer_hash_table_create(void **i) {
    (*i) = (void*)(new libsupermesh::integer_hash_table());
    if(!((libsupermesh::integer_hash_table*)(*i))) {
      libsupermesh_abort("new failure");
    }
    
    return;
  }
  
  void libsupermesh_integer_hash_table_delete(void **i) {
    delete ((libsupermesh::integer_hash_table*)(*i));
    
    return;
  }
  
  void libsupermesh_integer_hash_table_insert(void **i, const int *k, const int *v) {
    ((libsupermesh::integer_hash_table*)(*i))->insert(*k, *v);
    
    return;
  }
  
  void libsupermesh_integer_hash_table_length(const void **i, int *l) {    
    ((libsupermesh::integer_hash_table*)(*i))->length(*l);
    
    return;
  }
  
  void libsupermesh_integer_hash_table_fetch(void **i, const int *k, int *v) {
    ((libsupermesh::integer_hash_table*)(*i))->fetch(*k, *v);
    
    return;
  }
  
  void libsupermesh_integer_hash_table_remove(void **i, const int *k, int *status) {
    ((libsupermesh::integer_hash_table*)(*i))->remove(*k, *status);
    
    return;
  }
  
  void libsupermesh_integer_hash_table_has_key(const void **i, const int *k, int *present) {
    ((libsupermesh::integer_hash_table*)(*i))->has_key(*k, *present);

    return;
  }

  void libsupermesh_integer_hash_table_fetch_pair(void **i, const int *idx, int *k, int *v) {
    ((libsupermesh::integer_hash_table*)(*i))->fetch_pair(*idx, *k, *v);

    return;
  }
}

#endif
