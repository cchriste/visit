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
//                  avtMeshLogFilter.h                                       //
// ************************************************************************* //

#ifndef AVT_MESH_LOG_FILTER_H
#define AVT_MESH_LOG_FILTER_H

#include <filters_exports.h>

#include <avtDataTreeIterator.h>

#include <enumtypes.h>


// ****************************************************************************
//  Class: avtMeshLogFilter
//
//  Purpose:
//    Scales the mesh's coordinates by base 10 logarithm.
//
//  Programmer: Kathleen Bonnell 
//  Creation:   March 6, 2007 
//
//  Modifications:
//    Eric Brugger, Mon Jul 21 14:26:01 PDT 2014
//    Modified the class to work with avtDataRepresentation.
//
// ****************************************************************************

class AVTFILTERS_API avtMeshLogFilter : public avtDataTreeIterator
{
  public:
                          avtMeshLogFilter();
    virtual              ~avtMeshLogFilter(){}; 

    virtual const char   *GetType(void) {return "avtMeshLogFilter";};
    virtual const char   *GetDescription(void) 
                              { return "Calculating base 10 logarithm."; };

    void                  SetXScaleMode(ScaleMode val) { xScaleMode = val; }
    void                  SetYScaleMode(ScaleMode val) { yScaleMode = val; }

    void                  SetUseInvLogX(bool val) { useInvLogX = val; }
    void                  SetUseInvLogY(bool val) { useInvLogY = val; }

  protected:
    ScaleMode             xScaleMode;
    ScaleMode             yScaleMode;
    bool                  useInvLogX;
    bool                  useInvLogY;

    virtual void          PostExecute(void);
    virtual avtDataRepresentation *ExecuteData(avtDataRepresentation *);
    virtual void          UpdateDataObjectInfo(void);

  private:
};


#endif


