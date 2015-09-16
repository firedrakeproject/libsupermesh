/******************************************************************************
 * Project:  libsidx - A C API wrapper around libspatialindex
 * Purpose:  C++ object declarations to implement the wrapper.
 * Author:   Howard Butler, hobu.inc@gmail.com
 ******************************************************************************
 * Copyright (c) 2009, Howard Butler
 *
 * All rights reserved.
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
******************************************************************************/

#pragma once

class Index
{

public:
    Index(const LibSupermesh_Tools::PropertySet& poProperties);
    Index(const LibSupermesh_Tools::PropertySet& poProperties, int (*readNext)(LibSupermesh_SpatialIndex::id_type *id, double **pMin, double **pMax, uint32_t *nDimension, const uint8_t **pData, uint32_t *nDataLength));
    ~Index();

    const LibSupermesh_Tools::PropertySet& GetProperties() { return m_properties; }

    bool insertFeature(uint64_t id, double *min, double *max);
    
    RTIndexType GetIndexType();
    void SetIndexType(RTIndexType v);

    RTStorageType GetIndexStorage();
    void SetIndexStorage(RTStorageType v);
    
    RTIndexVariant GetIndexVariant();
    void SetIndexVariant(RTStorageType v);
    
    LibSupermesh_SpatialIndex::ISpatialIndex& index() {return *m_rtree;}
    LibSupermesh_SpatialIndex::StorageManager::IBuffer& buffer() {return *m_buffer;}

private:

    void Initialize();
    LibSupermesh_SpatialIndex::IStorageManager* m_storage;
    LibSupermesh_SpatialIndex::StorageManager::IBuffer* m_buffer;
    LibSupermesh_SpatialIndex::ISpatialIndex* m_rtree;
    
    LibSupermesh_Tools::PropertySet m_properties;


    void Setup();
    LibSupermesh_SpatialIndex::IStorageManager* CreateStorage();
    LibSupermesh_SpatialIndex::StorageManager::IBuffer* CreateIndexBuffer(LibSupermesh_SpatialIndex::IStorageManager& storage);
    LibSupermesh_SpatialIndex::ISpatialIndex* CreateIndex();
};
