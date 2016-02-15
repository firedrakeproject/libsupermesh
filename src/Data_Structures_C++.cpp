#include "libsupermesh_configuration.h"

#ifndef ENABLE_JUDY

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <map>
#include <set>

using namespace std;

namespace libsupermesh {

class integer_set {
  public:
    inline integer_set() {
      this->v_array = NULL;
      
      return;
    }
    
    inline ~integer_set() {
      if(this->v_array) {
        free(this->v_array);
      }
      
      return;
    }
    
    inline void insert(const int& v, int& c) {
      if(this->v_array) {
        free(this->v_array);
        this->v_array = NULL;
      }
      pair<set<int>::iterator, bool> position = this->v_set.insert(v);
      c = position.second ? 1 : 0;
      
      return;
    }
    
    inline void length(int& l) {
      l = this->v_set.size();
      
      return;
    };
    
    inline void fetch(const int& idx, int& val) {
      if(!this->v_array) {
        this->v_array = (int*)malloc(this->v_set.size() * sizeof(int));
        if(!this->v_array) {
          cerr << "malloc failure" << endl;
          exit(1);
        }
        size_t i = 0;
        for(set<int>::iterator iter = this->v_set.begin();iter != this->v_set.end();iter++) {
          this->v_array[i] = *iter;
          i++;
        }
      }
      val = this->v_array[idx - 1];
      
      return;
    };
    
    inline void remove(const int& idx, int& status) {
      if(this->v_array) {
        free(this->v_array);
        this->v_array = NULL;
      }
      set<int>::size_type count = this->v_set.erase(idx);
      status = (count > 0) ? 1 : 0;
      
      return;
    }
    
    inline void has_value(const int& val, int& present) {
      set<int>::iterator iter = this->v_set.find(val);
      present = (iter == this->v_set.end()) ? 0 : 1;
      
      return;
    }
    
  private:
    set<int> v_set;
    int *v_array;
};

class integer_hash_table {
  public:
    inline integer_hash_table() {
      this->k_array = NULL;
      
      return;
    }
    
    inline ~integer_hash_table() {
      if(this->k_array) {
        free(this->k_array);
      }
    }
    
    inline void insert(const int& k, const int& v) {
      if(this->k_array) {
        free(this->k_array);
        this->k_array = NULL;
      }
      this->v_map[k] = v;
      
      return;
    }
    
    inline void length(int& l) {
      l = this->v_map.size();
      
      return;
    }
    
    inline void fetch(const int& k, int& v) {
      v = this->v_map[k];
      
      return;
    }
    
    inline void remove(const int& k, int& status) {
      if(this->k_array) {
        free(this->k_array);
        this->k_array = NULL;
      }
      map<int, int>::size_type count = this->v_map.erase(k);
      status = (count > 0) ? 1 : 0;
      
      return;
    }
    
    inline void has_key(const int& k, int& present) {
      map<int, int>::iterator iter = this->v_map.find(k);
      present = (iter == this->v_map.end()) ? 0 : 1;
      
      return;
    }
    
    inline void fetch_pair(const int& idx, int& k, int& v) {
      if(!this->k_array) {
        this->k_array = (int*)malloc(this->v_map.size() * sizeof(int));
        if(!this->k_array) {
          cerr << "malloc failure" << endl;
          exit(1);
        }
        size_t i = 0;
        for(map<int, int>::iterator iter = this->v_map.begin();iter != this->v_map.end();iter++) {
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
  void libsupermesh_integer_set_create_c(void **i) {
    (*i) = (void*)(new libsupermesh::integer_set());
    if(!(*i)) {
      cerr << "new failure" << endl;
      exit(1);
    }
    
    return;
  }
  
  void libsupermesh_integer_set_delete_c(void **i) {
    delete ((libsupermesh::integer_set*)(*i));
    
    return;
  }
  
  void libsupermesh_integer_set_insert_c(void **i, int *v, int *c) {
    ((libsupermesh::integer_set*)(*i))->insert(*v, *c);
    
    return;
  }
  
  void libsupermesh_integer_set_length_c(void **i, int *l) {
    ((libsupermesh::integer_set*)(*i))->length(*l);
    
    return;
  }
  
  void libsupermesh_integer_set_fetch_c(void **i, int *idx, int *val) {
    ((libsupermesh::integer_set*)(*i))->fetch(*idx, *val);
    
    return;
  }
  
  void libsupermesh_integer_set_remove_c(void **i, int *idx, int *status) {
    ((libsupermesh::integer_set*)(*i))->remove(*idx, *status);
    
    return;
  }

  void libsupermesh_integer_set_has_value_c(void **i, int *val, int *present) {
    ((libsupermesh::integer_set*)(*i))->has_value(*val, *present);
    
    return;
  }
  
  void libsupermesh_integer_hash_table_create_c(void **i) {
    (*i) = (void*)(new libsupermesh::integer_hash_table);
    if(!(*i)) {
      cerr << "new failure" << endl;
      exit(1);
    }
    
    return;
  }
  
  void libsupermesh_integer_hash_table_delete_c(void **i) {
    delete ((libsupermesh::integer_hash_table*)(*i));
    
    return;
  }
  
  void libsupermesh_integer_hash_table_insert_c(void **i, int *k, int *v) {
    ((libsupermesh::integer_hash_table*)(*i))->insert(*k, *v);
    
    return;
  }
  
  void libsupermesh_integer_hash_table_length_c(void **i, int *l) {    
    ((libsupermesh::integer_hash_table*)(*i))->length(*l);
    
    return;
  }
  
  void libsupermesh_integer_hash_table_fetch_c(void **i, int *k, int *v) {
    ((libsupermesh::integer_hash_table*)(*i))->fetch(*k, *v);
    
    return;
  }
  
  void libsupermesh_integer_hash_table_remove_c(void **i, int *k, int *status) {
    ((libsupermesh::integer_hash_table*)(*i))->remove(*k, *status);
    
    return;
  }
  
  void libsupermesh_integer_hash_table_has_key_c(void **i, int *k, int *present) {
    ((libsupermesh::integer_hash_table*)(*i))->has_key(*k, *present);

    return;
  }

  void libsupermesh_integer_hash_table_fetch_pair_c(void **i, int *idx, int *key, int *val) {
    ((libsupermesh::integer_hash_table*)(*i))->fetch_pair(*idx, *key, *val);

    return;
  }
}

#endif