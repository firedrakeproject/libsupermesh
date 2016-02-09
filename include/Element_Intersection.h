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
#include "Wm4Intersector.h"
#include "Wm4IntrTriangle2Triangle2.h"
#include "Wm4Triangle2.h"
#include "Wm4IntrQuad2Quad2.h"
#include "Wm4Quad2.h"
#include "Wm4IntrTetrahedron3Tetrahedron3.h"
#include "Wm4Tetrahedron3.h"
#include "Wm4Vector3.h"

#include <strings.h>
#include <spatialindex/SpatialIndex.h>

#include <cassert>
#include <iostream>
#include <vector>

#define GEOM_REAL double

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

  class ElementIntersector
  {
    public:
      virtual ~ElementIntersector();
      
      virtual unsigned int GetDim() const = 0;
      
      virtual void SetInput(double*& positionsA, double*& positionsB, const int& dim, const int& loc);
      virtual void Intersect() = 0;
      virtual void QueryOutput(int& nnodes, int& nelms) const = 0; 
      virtual void GetOutput(double*& positions, int*& enlist) const = 0;
    protected:
      ElementIntersector();
      
      double* positionsA;
      double* positionsB;
      int loc;
      int dim;
  };
  
  class ElementIntersector1D : public ElementIntersector
  {
    public:
      ElementIntersector1D();
      virtual ~ElementIntersector1D();
      
      inline virtual unsigned int GetDim() const
      {
        return 1;
      }
      
      inline virtual void SetInput(double*& positionsA, double*& positionsB, const int& dim, const int& loc)
      {
        assert(dim == 1);
        assert(loc == 2);
        ElementIntersector::SetInput(positionsA, positionsB, dim, loc);
        
        return;
      }
      
      virtual void Intersect();
      virtual void QueryOutput(int& nnodes, int& nelms) const; 
      virtual void GetOutput(double*& positions, int*& enlist) const;
    protected:
      double positionsC[2];
      bool intersection;
  };

  class ElementIntersector2D : public ElementIntersector
  {
    public:
      ElementIntersector2D();
      virtual ~ElementIntersector2D();
      
      inline virtual unsigned int GetDim() const
      {
        return 2;
      }
      
      inline virtual void SetInput(double*& positionsA, double*& positionsB, const int& dim, const int& loc)
      {
        assert(dim == 2);
        assert(loc == 3 || loc == 4);
        ElementIntersector::SetInput(positionsA, positionsB, dim, loc);
      
        return;
      }
      virtual void Intersect();
      virtual void QueryOutput(int& nnodes, int& nelms) const; 
      virtual void GetOutput(double*& positions, int*& enlist) const;
      
      typedef LibSupermesh_Wm4::IntrTriangle2Triangle2<GEOM_REAL> IntrTriangle2Triangle2;
      typedef LibSupermesh_Wm4::Triangle2<GEOM_REAL> Triangle2;
      typedef LibSupermesh_Wm4::Vector2<GEOM_REAL> Vector2;
      typedef LibSupermesh_Wm4::IntrQuad2Quad2<GEOM_REAL> IntrQuad2Quad2;
      typedef LibSupermesh_Wm4::Quad2<GEOM_REAL> Quad2;
      typedef LibSupermesh_Wm4::Intersector<GEOM_REAL, Vector2> Intersector2d;

    protected:
       Intersector2d* intersection;
  };

  class ElementIntersector3D : public ElementIntersector
  {
    public:
      inline ElementIntersector3D() {return;}
      inline virtual ~ElementIntersector3D() {return;}

      inline virtual unsigned int GetDim() const
      {
        return 3;
      }
      
      inline virtual void SetInput(double*& positionsA, double*& positionsB, const int& dim, const int& loc)
      {
        assert(dim == 3);
      
        ElementIntersector::SetInput(positionsA, positionsB, dim, loc);
        
        return;
      }

      virtual void Intersect() = 0;
      virtual void QueryOutput(int& nnodes, int& nelms) const = 0; 
      virtual void GetOutput(double*& positions, int*& enlist) const = 0;
  };

  class WmElementIntersector3D : public ElementIntersector3D
  {
    public:
      WmElementIntersector3D();
      virtual ~WmElementIntersector3D();
      
      inline virtual void SetInput(double*& positionsA, double*& positionsB, const int& dim, const int& loc)
      {
        assert(loc == 4);
        ElementIntersector3D::SetInput(positionsA, positionsB, dim, loc);
        
        return;
      }

      virtual void Intersect();
      virtual void QueryOutput(int& nnodes, int& nelms) const; 
      virtual void GetOutput(double*& positions, int*& enlist) const;
      
      typedef LibSupermesh_Wm4::IntrTetrahedron3Tetrahedron3<GEOM_REAL> IntrTetrahedron3Tetrahedron3;
      typedef LibSupermesh_Wm4::Tetrahedron3<GEOM_REAL> Tetrahedron3;
      typedef LibSupermesh_Wm4::Vector3<GEOM_REAL> Vector3;
    protected:
      IntrTetrahedron3Tetrahedron3* intersection;
      std::vector<GEOM_REAL>* volumes;
      int nodes, elements;
  };
}

extern LibSupermesh::ElementIntersector* elementIntersector_LibSuperMesh;

extern LibSupermesh::ElementIntersectionFinder elementIntersectionFinder_LibSuperMesh;

extern "C"
{
#define cLibSuperMeshIntersectorGetDimension libsupermesh_cintersector_get_dimension
  int cLibSuperMeshIntersectorGetDimension();

#define cLibSuperMeshIntersectorSetDimension libsupermesh_cintersector_set_dimension
  void cLibSuperMeshIntersectorSetDimension(const int* dim);

#define cLibSuperMeshIntersectorSetInput libsupermesh_cintersector_set_input
  void cLibSuperMeshIntersectorSetInput(double* positionsA, double* positionsB, const int* dim, const int* loc);

#define cLibSuperMeshIntersectorDrive libsupermesh_cintersector_drive
  void cLibSuperMeshIntersectorDrive();
  
#define cLibSuperMeshIntersectorQuery libsupermesh_cintersector_query
  void cLibSuperMeshIntersectorQuery(int* nnodes, int* nelms);
  
#define cLibSuperMeshIntersectorGetOutput libsupermesh_cintersector_get_output
  void cLibSuperMeshIntersectorGetOutput(const int* nnodes, const int* nelms, const int* dim, const int* loc, double* positions, int* enlist);

  void libsupermesh_cintersection_finder_reset();
  void libsupermesh_cintersection_finder_set_input(const double* positions, const int* enlist, const int* dim, const int* loc, const int* nnodes, const int* nelements);  
  void libsupermesh_cintersection_finder_find(const double* positions, const int* dim, const int* loc);  
  void libsupermesh_cintersection_finder_query_output(int* nelms);
  void libsupermesh_cintersection_finder_get_output(int* id, const int* index);
}
