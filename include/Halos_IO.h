/*  Copyright (C) 2006 Imperial College London and others.
    
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

#pragma once

#include "libsupermesh_configuration.h"
#include "Tokenize.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <deque>
#include <iostream>
#include <map>
#include <set>
#include <string.h>
#include <vector>

#include "tinyxml.h"

namespace libsupermesh {

  enum HaloReadError{
    HALO_READ_SUCCESS = 0,
    HALO_READ_FILE_NOT_FOUND = -1,
    HALO_READ_FILE_INVALID = -2,
  };

  //* Read halo information
  /** Read from a halo file.
    * \param filename Halo file name
    * \param process The process number
    * \param nprocs The number of processes
    * \param npnodes Number of private nodes, by tag
    * \param send Sends, by tag and process
    * \param recv Receives, by tag and process
    * \return 0 on success, non-zero on failure
    */
  HaloReadError ReadHalos(const std::string& filename, int& process, int& nprocs, std::map<int, int>& npnodes, std::map<int, std::vector<std::vector<int> > >& send, std::map<int, std::vector<std::vector<int> > >& recv);
  
  struct HaloData{
      int process, nprocs;
      std::map<int, int> npnodes;
      std::map<int, std::vector<std::vector<int> > > send, recv;
  };
  
  extern HaloData* readHaloData;
}

extern "C" {
  void libsupermesh_halo_reader_reset();
  int libsupermesh_halo_reader_set_input(char* filename, int* filename_len, int* process, int* nprocs);
  void libsupermesh_halo_reader_query_output(int* level, int* nprocs, int* nsends, int* nreceives);
  void libsupermesh_halo_reader_get_output(int* level, int* nprocs, int* nsends, int* nreceives,
    int* npnodes, int* send, int* recv);
}
