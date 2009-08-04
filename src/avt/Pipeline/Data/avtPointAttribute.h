/*****************************************************************************
*
* Copyright (c) 2000 - 2009, Lawrence Livermore National Security, LLC
* Produced at the Lawrence Livermore National Laboratory
* LLNL-CODE-400124
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
//                            avtPointAttribute.h                            //
// ************************************************************************* //

#ifndef AVT_POINT_ATTRIBUTE_H
#define AVT_POINT_ATTRIBUTE_H

#include <pipeline_exports.h>

#include <string>

#include <avtDynamicAttribute.h>

#include <Point.h>


// ****************************************************************************
//  Class: avtPointAttribute
//
//  Purpose:
//      Contains the knowledge of how to convert a standard point attribute
//      into an actual point.
//
//  Programmer: Hank Childs
//  Creation:   March 15, 2002
//
//  Modifications:
//
//    Hank Childs, Thu Feb  5 17:11:06 PST 2004
//    Moved inlined constructor and destructor definitions to .C files
//    because certain compilers have problems with them.
//
// ****************************************************************************

class PIPELINE_API avtPointAttribute : public avtDynamicAttribute
{
  public:
                             avtPointAttribute(Point &);
    virtual                 ~avtPointAttribute();

    void                     GetPoint(double [3]);

    virtual avtContract_p
                             ModifyContract(avtContract_p);
    virtual void             ResolveAttributes(avtDataObject_p);

  protected:
    Point                   &attributes;
    double                   point[3];
    std::string              origVar;
    bool                     resolvedAttributes;

    void                     ResolvePoint(void);
    void                     ResolveMin(avtDataObject_p &);
    void                     ResolveMax(avtDataObject_p &);
    void                     ResolveZone(avtDataObject_p &);
    void                     ResolveNode(avtDataObject_p &);
};


#endif


