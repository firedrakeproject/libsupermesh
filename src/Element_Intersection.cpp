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

#include "Element_Intersection.h"

using namespace LibSupermesh;
using namespace libsupermesh::SpatialIndex;

MeshDataStream::MeshDataStream(const double*& positions, const int& nnodes, const int& dim,
                               const int*& enlist, const int& nelements, const int& loc)
{
  assert(positions);
  assert(nnodes >= 0);
  assert(dim >= 0);
  assert(enlist);
  assert(nelements >= 0);
#ifdef DDEBUG
  if(nelements > 0)
  {
    assert(loc >= 1);
  }
  for(int i = 0;i < nelements*loc;i++)
  {
    assert(enlist[i] <= nnodes);
  }
#endif

  this->positions = positions;
  this->nnodes = nnodes;
  this->dim = dim;
  this->enlist = enlist;
  this->nelements = nelements;
  this->loc = loc;

  index = 0;
}

MeshDataStream::~MeshDataStream()
{
}

IData* MeshDataStream::getNext()
{
  if(index >= nelements)
  {
    return NULL;
  }

  int node;

  double high[dim], low[dim];

  node = enlist[loc * index] - 1;

  for(int i = 0;i < dim;i++)
  {
    high[i] = positions[dim * node + i];
    low[i] = positions[dim * node + i];
  }
  for(int i = 1; i < loc; i++)
  {
    node = enlist[loc * index + i] - 1;
    for(int j = 0;j < dim;j++)
    {
      if(positions[dim * node + j] > high[j])
      {
        high[j] = positions[dim * node + j];
      }
      else if(positions[dim * node + j] < low[j])
      {
        low[j] = positions[dim * node + j];
      }
    }
  }
      
  libsupermesh::SpatialIndex::Region region = libsupermesh::SpatialIndex::Region(low, high, dim);
  IData* data = new RTree::Data(0, 0, region, ++index);

  return data;
}

bool MeshDataStream::hasNext()
{
  return index < nelements;
}

uint32_t MeshDataStream::size()
{
  return nelements;
}

void MeshDataStream::rewind()
{
  index = 0;
  
  return;
}

ElementIntersectionFinder::ElementIntersectionFinder()
{
  Initialise();  
  
  return;
}

ElementIntersectionFinder::~ElementIntersectionFinder()
{
  Free();
  
  return;
}

void ElementIntersectionFinder::Reset()
{
  Free();
  Initialise();
  
  return;
}

void ElementIntersectionFinder::SetInput(const double*& positions, const int& nnodes, const int& dim,
                                         const int*& enlist, const int& nelements, const int& loc)
{
  
  assert(positions);
  assert(enlist);
  assert(nnodes >= 0);
  assert(dim >= 0);
  assert(nelements >= 0);
  assert(loc >= 0);

  Reset();
  this->dim = dim;
  this->loc = loc;
  
  MeshDataStream stream(positions, nnodes, dim,
                        enlist, nelements, loc);  
  // As in regressiontest/rtree/RTreeBulkLoad.cc in spatialindex 1.2.0
  id_type id = 1;
  rTree = RTree::createAndBulkLoadNewRTree(RTree::BLM_STR, stream, *storageManager, fillFactor, indexCapacity, leafCapacity, dim, variant, id);
  
  return;
}

void ElementIntersectionFinder::SetTestElement(const double*& positions, const int& dim, const int& loc)
{
  assert(positions);
  assert(dim == this->dim);
  
  visitor.clear();
  
  double high[dim], low[dim];
  for(int i = 0;i < dim;i++)
  {
    high[i] = positions[i];
    low[i] = positions[i];
  }
  for(int i = 0;i < loc;i++)
  {
    for(int j = 0;j < dim;j++)
    {
      if(positions[dim * i + j] > high[j])
      {
        high[j] = positions[dim * i + j];
      }
      else if(positions[dim * i + j] < low[j])
      {
        low[j] = positions[dim * i + j];
      }
    }
  }
  
  libsupermesh::SpatialIndex::Region* region = new libsupermesh::SpatialIndex::Region(low, high, dim);
  rTree->intersectsWithQuery(*region, visitor);
  
  delete region;
  
  return;
}

void ElementIntersectionFinder::QueryOutput(int& nelms) const
{
  nelms = visitor.size();
  
  return;
}

void ElementIntersectionFinder::GetOutput(int& id, const int& index) const
{
  assert(index > 0);
  assert(index <= (int) visitor.size());
  
  id = visitor[index-1];

  return;
}

void ElementIntersectionFinder::Initialise()
{
  storageManager = StorageManager::createNewMemoryStorageManager();
  storage = StorageManager::createNewRandomEvictionsBuffer(*storageManager, capacity, writeThrough);
  rTree = NULL;
  
  dim = 0;
  loc = 0;
}

void ElementIntersectionFinder::Free()
{
  if(rTree)
  {
    delete rTree;
    rTree = NULL;
  }
  delete storage;
  delete storageManager;

  visitor.clear();
  
  return;
}

ElementIntersectionFinder elementIntersectionFinder_LibSuperMesh;

extern "C"
{
  void libsupermesh_cintersection_finder_reset()
  {
    elementIntersectionFinder_LibSuperMesh.Reset();
    
    return;
  }

  void libsupermesh_cintersection_finder_set_input(const double* positions, const int* enlist, const int* dim, const int* loc, const int* nnodes, const int* nelements)
  {
    assert(*dim >= 0);
    assert(*loc >= 0);
    assert(*nnodes >= 0);
    assert(*nelements >= 0);
    
    elementIntersectionFinder_LibSuperMesh.SetInput(positions, *nnodes, *dim, enlist, *nelements, *loc);
    
    return;
  }

  void libsupermesh_cintersection_finder_find(const double* positions, const int* dim, const int* loc)
  {
    assert(*dim >= 0);
    assert(*loc >= 0);
    
    elementIntersectionFinder_LibSuperMesh.SetTestElement(positions, *dim, *loc);
    
    return;
  }

  void libsupermesh_cintersection_finder_query_output(int* nelms)
  {
    elementIntersectionFinder_LibSuperMesh.QueryOutput(*nelms);
    
    return;
  }

  void libsupermesh_cintersection_finder_get_output(int* id, const int* index)
  {
    elementIntersectionFinder_LibSuperMesh.GetOutput(*id, *index);
    
    return;
  }
}
