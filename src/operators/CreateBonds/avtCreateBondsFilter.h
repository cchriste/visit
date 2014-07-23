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
//  File: avtCreateBondsFilter.h
// ************************************************************************* //

#ifndef AVT_CreateBonds_FILTER_H
#define AVT_CreateBonds_FILTER_H

#include <avtPluginDataTreeIterator.h>

#include <CreateBondsAttributes.h>

class vtkPoints;
class vtkDataSet;
class vtkPolyData;


// ****************************************************************************
//  Class: avtCreateBondsFilter
//
//  Purpose:
//      A plugin operator for CreateBonds.
//
//  Programmer: Jeremy Meredith
//  Creation:   August 29, 2006
//
//  Modifications:
//    Jeremy Meredith, Mon Feb 11 16:39:51 EST 2008
//    The manual bonding matches now supports wildcards, alleviating the need
//    for a "simple" mode.  (The default for the manual mode is actually
//    exactly what the simple mode was going to be anyway.)
//
//    Jeremy Meredith, Tue Jan 26 16:30:40 EST 2010
//    Separated bond distance check to a more manual version.  This allows
//    for some optimizations, particularly in the new "periodic bond check"
//    feature.
//
//    Jeremy Meredith, Wed Jan 27 16:47:10 EST 2010
//    Added a new, fast bond creation routine.  Kept the old one around
//    as a fallback in case bonding distance paramemters make the fast
//    routine a poor choice, or in case we need to add some feature which
//    be too hard to add in the more complex fast routine.
//
//    Eric Brugger, Wed Jul 23 11:24:54 PDT 2014
//    Modified the class to work with avtDataRepresentation.
//
// ****************************************************************************

class avtCreateBondsFilter : public avtPluginDataTreeIterator
{
  public:
                         avtCreateBondsFilter();
    virtual             ~avtCreateBondsFilter();

    static avtFilter    *Create();

    virtual const char  *GetType(void)  { return "avtCreateBondsFilter"; };
    virtual const char  *GetDescription(void)
                             { return "CreateBonds"; };

    virtual void         SetAtts(const AttributeGroup*);
    virtual bool         Equivalent(const AttributeGroup*);

  protected:
    CreateBondsAttributes   atts;

    virtual avtDataRepresentation *ExecuteData(avtDataRepresentation *);
    vtkDataSet           *ExecuteData_Fast(vtkPolyData*, float maxBondDist,
                                          float,float,float,float,float,float);
    vtkDataSet           *ExecuteData_Slow(vtkPolyData*);
    virtual avtContract_p
                          ModifyContract(avtContract_p spec);
    virtual void          UpdateDataObjectInfo(void);

    bool AtomBondDistances(int eA, int eB, double &dmin, double &dmax);
};


#endif
