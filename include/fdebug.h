!    Copyright (C) 2006 Imperial College London and others.
!    
!    Please see the AUTHORS file in the main source directory for a full list
!    of copyright holders.
!
!    Prof. C Pain
!    Applied Modelling and Computation Group
!    Department of Earth Science and Engineering
!    Imperial College London
!
!    amcgsoftware@imperial.ac.uk
!    
!    This library is free software; you can redistribute it and/or
!    modify it under the terms of the GNU Lesser General Public
!    License as published by the Free Software Foundation,
!    version 2.1 of the License.
!
!    This library is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!    Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public
!    License along with this library; if not, write to the Free Software
!    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
!    USA

#pragma once

#include "confdefs.h"

#ifndef __FILE__
#define __FILE__ "unknown"
#endif

#ifndef __LINE__
#define __LINE__ "unknown"
#endif

#define ewrite(priority, format) if(priority <= current_debug_level) write(debug_unit(priority), format)

#define FLAbort(X) call FLAbort_pinpoint(X, __FILE__, __LINE__)
#define FLExit(X) call FLExit_pinpoint(X, __FILE__, __LINE__)

#ifdef NDEBUG
#define assert(X)
#else
#define assert(X) if(.not. (X)) FLAbort("Failed assertion " // "X")
#endif
