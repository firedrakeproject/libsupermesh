/*
  For libsupermesh copyright information see COPYING in the libsupermesh root
  directory

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
 * The following code is derived from include/fdebug.h and debug/C++_Debug.cpp
 * in Fluidity git revision 4e6c1d2b022df3a519cdec120fad28e60d1b08d9 (dated
 * 2015-02-25)
 */
 
/*
 * Fluidity copyright information (note that AUTHORS mentioned in the following
 * has been renamed to Fluidity_AUTHORS):
 */

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

#ifndef LIBSUPERMESH_LIBSUPERMESH_DEBUG_C_H
#define LIBSUPERMESH_LIBSUPERMESH_DEBUG_C_H

#include "libsupermesh_configuration.h"

#ifdef __cplusplus
extern "C" {
#endif
  void libsupermesh_print_backtrace(int max_size);
  void libsupermesh_abort_pinpoint(const char *error, const char *file, int line_number);
#ifdef __cplusplus
}
#endif

#define libsupermesh_abort(X) libsupermesh_abort_pinpoint(X, __FILE__, __LINE__)

#ifdef assert
#undef assert
#endif
#ifdef NDEBUG
#define assert(X)
#else
#define assert(X) if(!(X)) libsupermesh_abort("Failed assertion " #X)
#endif

#endif
