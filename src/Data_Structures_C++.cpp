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
    
    inline void insert(const int& v, int& changed) {
      if(this->v_array) {
        free(this->v_array);
        this->v_array = NULL;
      }
      pair<set<int>::const_iterator, bool> position = this->v_set.insert(v);
      changed = position.second ? 1 : 0;
      
      return;
    }
    
    inline void length(int& l) const {
      l = this->v_set.size();
      
      return;
    };
    
    inline void fetch(const int& idx, int& v) {
      if(!this->v_array) {
        this->v_array = (int*)malloc(this->v_set.size() * sizeof(int));
        if(!this->v_array) {
          cerr << "malloc failure" << endl;
          exit(1);
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
    
    inline void remove(const int& v, int& status) {
      if(this->v_array) {
        free(this->v_array);
        this->v_array = NULL;
      }
      set<int>::size_type count = this->v_set.erase(v);
      status = (count > 0) ? 1 : 0;
      
      return;
    }
    
    inline void has_value(const int& v, int& present) const {
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
    
    inline void insert(const int& k, const int& v) {
      if(this->k_array) {
        free(this->k_array);
        this->k_array = NULL;
      }
      this->v_map[k] = v;
      
      return;
    }
    
    inline void length(int& l) const {
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
    
    inline void has_key(const int& k, int& present) const {
      map<int, int>::const_iterator iter = this->v_map.find(k);
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
      cerr << "new failure" << endl;
      exit(1);
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
      cerr << "new failure" << endl;
      exit(1);
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