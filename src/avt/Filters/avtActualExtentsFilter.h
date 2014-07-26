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
//                          avtActualExtentsFilter.h                         //
// ************************************************************************* //

#ifndef AVT_ACTUAL_EXTENTS_FILTER_H
#define AVT_ACTUAL_EXTENTS_FILTER_H

#include <filters_exports.h>

#include <avtDatasetToDatasetFilter.h>


// ****************************************************************************
//  Class: avtActualExtentsFilter
//
//  Purpose:
//    Calculates the actual extents, both spatial and data.  Stores them
//    in the output's info. 
//
//  Programmer: Kathleen Bonnell 
//  Creation:   October 2, 2001 
//
//  Modifications:
//    Jeremy Meredith, Thu Feb 15 11:44:28 EST 2007
//    Added support for rectilinear grids with an inherent transform.
//
//    Hank Childs, Tue Aug 26 14:37:27 PDT 2008
//    Implement ModifyContract to prevent base class from declaring that it
//    can only work on floats.
//
//    Hank Childs, Thu Aug 26 13:47:30 PDT 2010
//    Renamed to avtActualExtentsFilter.
//
// ****************************************************************************

class AVTFILTERS_API avtActualExtentsFilter : public avtDatasetToDatasetFilter
{
  public:
                          avtActualExtentsFilter(){};
    virtual              ~avtActualExtentsFilter(){}; 

    virtual const char   *GetType(void) {return "avtActualExtentsFilter";};
    virtual const char   *GetDescription(void) 
                              { return "Calculating Actual Extents."; };

  protected:
    avtContract_p         lastContract;

    virtual void          Execute(void);
    virtual void          UpdateExtents(void);
    virtual bool          FilterUnderstandsTransformedRectMesh();
    virtual avtContract_p ModifyContract(avtContract_p);
};


#endif


