/*****************************************************************************
*
* Copyright (c) 2000 - 2016, Lawrence Livermore National Security, LLC
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
//  File: avtCoordSystemConvert.h
// ************************************************************************* //

#ifndef AVT_CoordConvert_FILTER_H
#define AVT_CoordConvert_FILTER_H

#include <filters_exports.h>

#include <avtDataTreeIterator.h>

// ****************************************************************************
//  Class: avtCoordSystemConvert
//
//  Purpose:
//      A plugin operator for CoordConvert.
//
//  Programmer: childs -- generated by xml2info
//  Creation:   Fri Jun 27 16:41:32 PST 2003
//
//  Modifications:
//
//    Hank Childs, Wed Jun  8 14:50:52 PDT 2005
//    Added UpdateDataObjectInfo.
//
//    Jeremy Meredith, Fri Aug  7 15:33:18 EDT 2009
//    Added selectable vector transform method.
//    Made the coord system enum be within the class namespace.
//
//    Eric Brugger, Mon Jul 21 10:19:29 PDT 2014
//    Modified the class to work with avtDataRepresentation.
//
// ****************************************************************************

class AVTFILTERS_API avtCoordSystemConvert : public avtDataTreeIterator
{
  public:
    enum CoordSystem
    {
        CARTESIAN = 0,                /* 0 */
        CYLINDRICAL,                  /* 1 */
        SPHERICAL                     /* 2 */
    };

    enum VectorTransformMethod
    {
        None,
        AsPoint,
        AsDisplacement,
        AsDirection
    };

                         avtCoordSystemConvert();
    virtual             ~avtCoordSystemConvert();

    virtual const char  *GetType(void)  { return "avtCoordSystemConvert"; };
    virtual const char  *GetDescription(void)
                             { return "Converting Coordinate System"; };

    void                 SetInputCoordSys(CoordSystem a) { inputSys = a; };
    void                 SetOutputCoordSys(CoordSystem a) { outputSys = a; };
    void                 SetContinuousPhi(bool a) { continuousPhi = a; };
    void                 SetVectorTransformMethod(VectorTransformMethod m)
                                               { vectorTransformMethod = m; }

  protected:
    CoordSystem           inputSys;
    CoordSystem           outputSys;
    bool                  continuousPhi;

    VectorTransformMethod vectorTransformMethod;

    virtual avtDataRepresentation *ExecuteData(avtDataRepresentation *);
    virtual void          PostExecute(void);

    virtual void          UpdateDataObjectInfo(void);
};
#endif
