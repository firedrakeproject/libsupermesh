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

#include "Data_Structures_C.h"

#ifndef LIBSUPERMESH_ENABLE_JUDY

#include <iostream>
#include <map>
#include <set>

#include "libsupermesh_debug_C.h"

using namespace std;

namespace libsupermesh {

class integer_set {
  public:
    inline integer_set(void) {
      this->value_array = NULL;
    }
    
    inline ~integer_set(void) {
      if(this->value_array) {
        delete[] this->value_array;
      }
    }
    
    inline void insert(const int &value, int &changed) {
      pair<set<int>::const_iterator, bool> position = this->value_set.insert(value);
      changed = position.second ? 1 : 0;
      if(changed && this->value_array) {
        delete[] this->value_array;
        this->value_array = NULL;
      }
    }
    
    inline void size(int &size) const {
      size = this->value_set.size();
    };
    
    inline void fetch(const int &index, int &value) {
      if(!this->value_array) {
        this->value_array = new int[this->value_set.size()];
        set<int>::size_type i = 0;
        for(set<int>::const_iterator iter = this->value_set.begin();iter != this->value_set.end();iter++) {
          this->value_array[i++] = *iter;
        }
      }
      if(index < 0 || index > this->value_set.size()) {
        std::cerr << "Failed to fetch integer set element with index " << index << endl;
        libsupermesh_abort("Failed to fetch integer set element");
      }
      value = this->value_array[index - 1];
    };
    
    inline void remove(const int &value) {
      if(this->value_array) {
        delete[] this->value_array;
        this->value_array = NULL;
      }
      set<int>::size_type count = this->value_set.erase(value);
      if(count == 0) {
        cerr << "Failed to remove integer set element with value " << value << endl;
        libsupermesh_abort("Failed to remove integer set element");
      }
    }
    
    inline void has_value(const int &value, int &present) const {
      set<int>::const_iterator iter = this->value_set.find(value);
      present = (iter == this->value_set.end()) ? 0 : 1;
    }
    
  private:
    set<int> value_set;
    int *value_array;
};

class integer_hash_table {
  public:
    inline integer_hash_table(void) {
      this->key_array = NULL;
    }
    
    inline ~integer_hash_table(void) {
      if(this->key_array) {
        delete[] this->key_array;
      }
    }
    
    inline void insert(const int &key, const int &value) {
      if(this->key_array) {
        delete[] this->key_array;
        this->key_array = NULL;
      }
      this->value_map[key] = value;
    }
    
    inline void size(int &size) const {
      size = this->value_map.size();
    }
    
    inline void fetch(const int &key, int &value) {
      if(this->value_map.count(key) == 0) {
        cerr << "Failed to fetch integer hash table element with key " << key << endl;
        libsupermesh_abort("Failed to fetch integer hash table element");
      }
      value = this->value_map[key];
    }
    
    inline void remove(const int &key) {
      if(this->key_array) {
        delete[] this->key_array;
        this->key_array = NULL;
      }
      map<int, int>::size_type count = this->value_map.erase(key);
      if(count == 0) {
        cerr << "Failed to remove integer hash table element with key " << key << endl;
        libsupermesh_abort("Failed to remove integer hash table element");
      }
    }
    
    inline void has_key(const int &key, int &present) const {
      map<int, int>::const_iterator iter = this->value_map.find(key);
      present = (iter == this->value_map.end()) ? 0 : 1;
    }
    
    inline void fetch_pair(const int &index, int &key, int &value) {
      if(!this->key_array) {
        this->key_array = new int[this->value_map.size()];
        map<int, int>::size_type i = 0;
        for(map<int, int>::const_iterator iter = this->value_map.begin();iter != this->value_map.end();iter++) {
          this->key_array[i++] = iter->first;
        }
      }
      if(index < 1 || index > this->value_map.size()) {
        cerr << "Failed to fetch integer hash table element with index " << index << endl;
        libsupermesh_abort("Failed to fetch integer hash table element");
      }
      key = this->key_array[index - 1];
      value = this->value_map[key];
    }
    
  private:
    map<int, int> value_map;
    int *key_array;
};

}

extern "C" {
  void libsupermesh_integer_set_new(void **i) {
    *i = static_cast<void*>(new libsupermesh::integer_set());
  }
  
  void libsupermesh_integer_set_delete(void **i) {
    delete (static_cast<libsupermesh::integer_set*>(*i));
  }
  
  void libsupermesh_integer_set_insert(void **i, int value, int *changed) {
    (static_cast<libsupermesh::integer_set*>(*i))->insert(value, *changed);
  }
  
  void libsupermesh_integer_set_size(void **i, int *size) {
    (static_cast<libsupermesh::integer_set*>(*i))->size(*size);
  }
  
  void libsupermesh_integer_set_fetch(void **i, int index, int *value) {
    (static_cast<libsupermesh::integer_set*>(*i))->fetch(index, *value);
  }
  
  void libsupermesh_integer_set_remove(void **i, int value) {
    (static_cast<libsupermesh::integer_set*>(*i))->remove(value);
  }

  void libsupermesh_integer_set_has_value(void **i, int value, int *present) {
    (static_cast<libsupermesh::integer_set*>(*i))->has_value(value, *present);
  }
  
  void libsupermesh_integer_hash_table_new(void **i) {
    *i = static_cast<void*>(new libsupermesh::integer_hash_table());
  }
  
  void libsupermesh_integer_hash_table_delete(void **i) {
    delete (static_cast<libsupermesh::integer_hash_table*>(*i));
  }
  
  void libsupermesh_integer_hash_table_insert(void **i, int key, int value) {
    (static_cast<libsupermesh::integer_hash_table*>(*i))->insert(key, value);
  }
  
  void libsupermesh_integer_hash_table_size(void **i, int *size) {    
    (static_cast<libsupermesh::integer_hash_table*>(*i))->size(*size);
  }
  
  void libsupermesh_integer_hash_table_fetch(void **i, int key, int *value) {
    (static_cast<libsupermesh::integer_hash_table*>(*i))->fetch(key, *value);
  }
  
  void libsupermesh_integer_hash_table_remove(void **i, int key) {
    (static_cast<libsupermesh::integer_hash_table*>(*i))->remove(key);
  }
  
  void libsupermesh_integer_hash_table_has_key(void **i, int key, int *present) {
    (static_cast<libsupermesh::integer_hash_table*>(*i))->has_key(key, *present);
  }

  void libsupermesh_integer_hash_table_fetch_pair(void **i, int index, int *key, int *value) {
    (static_cast<libsupermesh::integer_hash_table*>(*i))->fetch_pair(index, *key, *value);
  }
}

#endif
