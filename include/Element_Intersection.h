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

#include "confdefs.h"
#include "MeshDataStream.h"

#include <strings.h>
#include <spatialindex/SpatialIndex.h>

#include <cassert>
#include <vector>

namespace LibSupermesh
{
  // StorageManager parameters
  const int capacity = 10;
  const int writeThrough = false;

  // R-Tree parameters
  const LibSupermesh_SpatialIndex::RTree::RTreeVariant variant = LibSupermesh_SpatialIndex::RTree::RV_RSTAR;
  // Minimum fraction (of maximum) of entries in any node (index or leaf)
  const double fillFactor = 0.7;
  // Node index capacity in the rtree
  const unsigned long indexCapacity = 10;
  // Node leaf capacity in the rtree
  const unsigned long leafCapacity = 10;
  
  // Customised version of PyListVisitor class in
  // wrapper.cc in Rtree 0.4.1
  class ElementListVisitor : public LibSupermesh_SpatialIndex::IVisitor, public std::vector< int >
  {
    public:
      inline ElementListVisitor() {}      
      inline virtual ~ElementListVisitor() {}      
      inline virtual void visitNode(const LibSupermesh_SpatialIndex::INode& node) {}      
      inline virtual void visitData(const LibSupermesh_SpatialIndex::IData& data) {push_back(data.getIdentifier());}
      inline virtual void visitData(std::vector< const LibSupermesh_SpatialIndex::IData* >& vector) {}
  };

  // Interface to spatialindex to calculate element intersection lists between
  // meshes using bulk storage
  // Uses code from gispatialindex.{cc,h} in Rtree 0.4.1
  class ElementIntersectionFinder
  {
    public:
      ElementIntersectionFinder();
      ~ElementIntersectionFinder();

      void Reset();
      void SetInput(const double*& positions, const int& nnodes, const int& dim,
                    const int*& enlist, const int& nelements, const int& loc);
      void SetTestElement(const double*& positions, const int& dim, const int& loc);
      void QueryOutput(int& nelms) const;
      void GetOutput(int& id, const int& index) const;
    protected:
      void Initialise();
      void Free();
    
      int dim, loc;
      LibSupermesh_SpatialIndex::IStorageManager* storageManager;
      LibSupermesh_SpatialIndex::StorageManager::IBuffer* storage;
      LibSupermesh_SpatialIndex::ISpatialIndex* rTree;
      ElementListVisitor visitor;
  };
}

extern LibSupermesh::ElementIntersectionFinder elementIntersectionFinder_LibSuperMesh;

extern "C"
{
  void libsupermesh_cintersection_finder_reset();
  void libsupermesh_cintersection_finder_set_input(const double* positions, const int* enlist, const int* dim, const int* loc, const int* nnodes, const int* nelements);  
  void libsupermesh_cintersection_finder_find(const double* positions, const int* dim, const int* loc);  
  void libsupermesh_cintersection_finder_query_output(int* nelms);
  void libsupermesh_cintersection_finder_get_output(int* id, const int* index);
}
