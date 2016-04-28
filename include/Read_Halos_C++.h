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
 * The following code is derived from include/Halos_IO.h and include/Tokenize.h
 * in Fluidity git revision 4e6c1d2b022df3a519cdec120fad28e60d1b08d9
 * (dated 2015-02-25)
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

#ifndef LIBSUPERMESH_READ_HALOS_CPP_H
#define LIBSUPERMESH_READ_HALOS_CPP_H

#include "libsupermesh_configuration.h"
#include "tinyxml.h"

#include <cassert>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace libsupermesh {  
  void ReadHalos(const std::string &filename, int &process, int &nprocs,
    std::map<int, int> &npnodes,
    std::map<int, std::vector<std::vector<int> > > &send,
    std::map<int, std::vector<std::vector<int> > > &recv);
  
  struct HaloData {
      int process, nprocs;
      std::map<int, int> npnodes;
      std::map<int, std::vector<std::vector<int> > > send, recv;
  };
  
  void Tokenize(const std::string &str, std::vector<std::string> &tokens, const std::string &delimiters = " ");
}

extern "C" {
  void libsupermesh_read_halo(void **data, const char *basename,
    const int *basename_len, const int *process, const int *nprocs);
  void libsupermesh_halo_sizes(const void **data, const int *level, 
    const int *nprocs, int *nsends, int *nreceives);
  void libsupermesh_halo_data(const void **data, const int *level,
    const int *nprocs, const int *nsends, const int *nreceives, int *npnodes,
    int *send, int *recv);
  void libsupermesh_deallocate_halo(void **data);
}

#endif
