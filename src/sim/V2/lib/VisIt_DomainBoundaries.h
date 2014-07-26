/*****************************************************************************
*
* Copyright (c) 2000 - 2014, Lawrence Livermore National Security, LLC
* Produced at the Lawrence Livermore National Laboratory
* LLNL-CODE-442911
* All rights reserved.
*
* This file is  part of VisIt. For  details, see https://visit.llnl.gov/.  The
* full copyright notice is contained in the file COPYRIGHT located at the root
* of the VisIt distribution or at http://www.llnl.gov/visit/copyright.html.
*
* Redistribution  and  use  in  source  and  binary  forms,  with  or  without
* modification, are permitted provided that the following conditions are met:
*
*  - Redistributions of  source code must  retain the above  copyright notice,
*    this list of conditions and the disclaimer below.
*  - Redistributions in binary form must reproduce the above copyright notice,
*    this  list of  conditions  and  the  disclaimer (as noted below)  in  the
*    documentation and/or other materials provided with the distribution.
*  - Neither the name of  the LLNS/LLNL nor the names of  its contributors may
*    be used to endorse or promote products derived from this software without
*    specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT  HOLDERS AND CONTRIBUTORS "AS IS"
* AND ANY EXPRESS OR  IMPLIED WARRANTIES, INCLUDING,  BUT NOT  LIMITED TO, THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND  FITNESS FOR A PARTICULAR  PURPOSE
* ARE  DISCLAIMED. IN  NO EVENT  SHALL LAWRENCE  LIVERMORE NATIONAL  SECURITY,
* LLC, THE  U.S.  DEPARTMENT OF  ENERGY  OR  CONTRIBUTORS BE  LIABLE  FOR  ANY
* DIRECT,  INDIRECT,   INCIDENTAL,   SPECIAL,   EXEMPLARY,  OR   CONSEQUENTIAL
* DAMAGES (INCLUDING, BUT NOT  LIMITED TO, PROCUREMENT OF  SUBSTITUTE GOODS OR
* SERVICES; LOSS OF  USE, DATA, OR PROFITS; OR  BUSINESS INTERRUPTION) HOWEVER
* CAUSED  AND  ON  ANY  THEORY  OF  LIABILITY,  WHETHER  IN  CONTRACT,  STRICT
* LIABILITY, OR TORT  (INCLUDING NEGLIGENCE OR OTHERWISE)  ARISING IN ANY  WAY
* OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
* DAMAGE.
*
*****************************************************************************/

#ifndef VISIT_DOMAINBOUNDARIES_H
#define VISIT_DOMAINBOUNDARIES_H

#ifdef __cplusplus
extern "C"
{
#endif

int VisIt_DomainBoundaries_alloc(visit_handle*);
int VisIt_DomainBoundaries_free(visit_handle);
/* Pass 0=rectilinear, 1=curvilinear for the type. */
int VisIt_DomainBoundaries_set_type(visit_handle, int);
int VisIt_DomainBoundaries_set_numDomains(visit_handle, int);

/* Set extents for structured mesh. */
int VisIt_DomainBoundaries_set_rectIndices(visit_handle, int dom, const int e[6]);

/* Set extents for an AMR patch. */
int VisIt_DomainBoundaries_set_amrIndices(visit_handle, int patch, int level, const int e[6]);

/* Functions to add a custom number of neighbors for a domain. */
int VisIt_DomainBoundaries_set_extents(visit_handle, int dom, const int e[6]);
int VisIt_DomainBoundaries_add_neighbor(visit_handle, int dom, int d, int mi, 
                                        const int orientation[3], const int extents[6]);
int VisIt_DomainBoundaries_finish(visit_handle, int dom);

#ifdef __cplusplus
}
#endif

#endif
