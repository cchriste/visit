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
//                            avtTruecolorFilter.h                           //
// ************************************************************************* //

#ifndef AVT_Truecolor_FILTER_H
#define AVT_Truecolor_FILTER_H

#include <avtDataTreeIterator.h>


// ****************************************************************************
//  Class: avtTruecolorFilter
//
//  Purpose:
//      This operator is the implied operator associated with an Truecolor plot.
//
//  Programmer: Chris Wojtan
//  Creation:   Monday June 15, 2004
//
//  Modifications:
//
//    Chris Wojtan Mon Jun 21 15:45 PDT 2004
//    Added "variable_name" member variable
//    Added SetVarName member function
//
//    Hank Childs, Fri May 20 14:55:06 PDT 2005
//    Remove UpdateDataObjectInfo (it was empty).
//
//    Eric Brugger, Tue Aug 19 11:46:04 PDT 2014
//    Modified the class to work with avtDataRepresentation.
//
// ****************************************************************************

class avtTruecolorFilter : public avtDataTreeIterator
{
  public:
                              avtTruecolorFilter();
    virtual                  ~avtTruecolorFilter();

    virtual const char       *GetType(void)   { return "avtTruecolorFilter"; };
    virtual const char       *GetDescription(void)
                                  { return "Performing Truecolor"; };
    void                      SetVarName(const char*name)
                                  {variable_name = name;}

  protected:
    virtual avtDataRepresentation *ExecuteData(avtDataRepresentation *);
    const char               *variable_name;
};


#endif


