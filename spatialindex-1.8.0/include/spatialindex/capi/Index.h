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
    Index(const libsupermesh::Tools::PropertySet& poProperties);
    Index(const libsupermesh::Tools::PropertySet& poProperties, int (*readNext)(libsupermesh::SpatialIndex::id_type *id, double **pMin, double **pMax, uint32_t *nDimension, const uint8_t **pData, uint32_t *nDataLength));
    ~Index();

    const libsupermesh::Tools::PropertySet& GetProperties() { return m_properties; }

    bool insertFeature(uint64_t id, double *min, double *max);
    
    RTIndexType GetIndexType();
    void SetIndexType(RTIndexType v);

    RTStorageType GetIndexStorage();
    void SetIndexStorage(RTStorageType v);
    
    RTIndexVariant GetIndexVariant();
    void SetIndexVariant(RTStorageType v);
    
    libsupermesh::SpatialIndex::ISpatialIndex& index() {return *m_rtree;}
    libsupermesh::SpatialIndex::StorageManager::IBuffer& buffer() {return *m_buffer;}

private:

    void Initialize();
    libsupermesh::SpatialIndex::IStorageManager* m_storage;
    libsupermesh::SpatialIndex::StorageManager::IBuffer* m_buffer;
    libsupermesh::SpatialIndex::ISpatialIndex* m_rtree;
    
    libsupermesh::Tools::PropertySet m_properties;


    void Setup();
    libsupermesh::SpatialIndex::IStorageManager* CreateStorage();
    libsupermesh::SpatialIndex::StorageManager::IBuffer* CreateIndexBuffer(libsupermesh::SpatialIndex::IStorageManager& storage);
    libsupermesh::SpatialIndex::ISpatialIndex* CreateIndex();
};
