/*
 * The following code is derived from include/Halos_IO.h and include/Tokenize.h
 * in Fluidity git revision 4e6c1d2b022df3a519cdec120fad28e60d1b08d9 (dated
 * 2015-02-25)
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
*/

#pragma once

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
