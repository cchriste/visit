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
//                    avtSphericalCompactnessFactorQuery.h                   //
// ************************************************************************* //

#ifndef AVT_SPHERICAL_COMPACTNESS_FACTOR_QUERY_H
#define AVT_SPHERICAL_COMPACTNESS_FACTOR_QUERY_H

#include <avtTwoPassDatasetQuery.h>

#include <query_exports.h>

class vtkDataSet;
class vtkCell;

class avtBinaryMultiplyFilter;
class avtVMetricVolume;
class avtRevolvedVolume;


// ****************************************************************************
//  Class: avtSphericalCompactnessFactorQuery
//
//  Purpose:
//      Calculates the spherical compactness factor of a data set.  This will
//      calculate the volume of the data set.  It will then calculate what the
//      radius of that volume would be if it was a perfect sphere.  Then it
//      will calculate what percentage of the data set is within that sphere.
//      (Note that this last step is the hardest, since it requires determining
//      where the center of the idealized sphere should be placed.)
//
//  Programmer: Hank Childs
//  Creation:   July 14, 2005
//
//  Modifications:
//    Cyrus Harrison, Wed Jul 16 15:52:57 PDT 2014
//    Added support for user selected center.
//
// ****************************************************************************

class QUERY_API avtSphericalCompactnessFactorQuery 
    : public avtTwoPassDatasetQuery
{
  public:
                                    avtSphericalCompactnessFactorQuery();
    virtual                        ~avtSphericalCompactnessFactorQuery();

    virtual const char             *GetType(void)
                         {return "avtSphericalCompactnessFactorQuery";};
    virtual const char             *GetDescription(void)
                         {return "Calculating Spherical Compactness Factor";};

    virtual void              SetInputParams(const MapNode &);

  protected:

    bool                            overrideCentroid;

    double                          centroid[3];
    double                          sphere_center[3];
    double                          radius;
    double                          total_volume;
    double                          volume_inside;
    bool                            is2D;
    avtRevolvedVolume              *rev_volume;
    avtVMetricVolume               *volume;

    virtual void                    Execute1(vtkDataSet *, const int);
    virtual void                    Execute2(vtkDataSet *, const int);
    virtual void                    PreExecute(void);
    virtual void                    MidExecute(void);
    virtual void                    PostExecute(void);
    virtual avtDataObject_p         ApplyFilters(avtDataObject_p);
};


#endif


