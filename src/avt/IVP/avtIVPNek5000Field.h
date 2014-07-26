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

// ************************************************************************* //
//                            avtIVPNek5000Field.h                             //
// ************************************************************************* //

#ifndef AVT_IVP_NEK5000_FIELD_H
#define AVT_IVP_NEK5000_FIELD_H

#include "visit-config.h"

#include "avtIVPVTKField.h"

#include <vtkDataSet.h>
#include <vtkPointData.h>

#include <ivp_exports.h>

#include <map>
#include <vector> 

struct findpts_local_data;

// ****************************************************************************
//  Class:  avtIVPNek5000Field
//
//  Purpose:
//    A wrapper class to allow the use of vtkDataSets as IVP fields for 
//    streamline integration. Uses vtkInterpolatedVelocityField on top of 
//    the supplied vtkDataSet. 
//
//  Creationist: Allen Sanderson
//  Creation:    May 1, 2013
//
// ****************************************************************************

class IVP_API avtIVPNek5000Field: public avtIVPVTKField
{
 public:
  avtIVPNek5000Field( vtkDataSet* ds, avtCellLocator* loc ); 
  avtIVPNek5000Field( float *elementsPtr, int nelements, int dim, int planes );

  ~avtIVPNek5000Field();

  Result operator()( const double &t, const avtVector &v, avtVector &retV ) const;

  avtVector ConvertToCartesian(const avtVector& pt) const;
  avtVector ConvertToCylindrical(const avtVector& pt) const;

 protected:

  // Storage for the Nek5000 mesh and the vector.
  double *nek_pts[3], *nek_vec[3];

  // Variables needed for Nek5000 find points
  mutable struct findpts_local_data *nek_fld;
};

#endif
