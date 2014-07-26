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
//                         avtShiftCenteringFilter.h                         //
// ************************************************************************* //

#ifndef AVT_SHIFT_CENTERING_FILTER_H
#define AVT_SHIFT_CENTERING_FILTER_H

#include <filters_exports.h>

#include <avtDataTreeIterator.h>


// ****************************************************************************
//  Class: avtShiftCenteringFilter
//
//  Purpose:
//    A filter which creates node-centered data from point-centered data or
//    vice-versa depending upon the desired centering. 
//
//  Programmer: Kathleen Bonnell 
//  Creation:   April 19, 2001 
//
//  Modifications:
//
//    Jeremy Meredith, Fri Jun 29 15:11:44 PDT 2001
//    Added description.
//
//    Hank Childs, Fri Feb  8 18:33:39 PST 2002
//    Removed all references to PC atts.
//
//    Hank Childs, Wed Feb 27 13:03:32 PST 2002
//    Added UpdateDataObjectInfo.
//
//    Hank Childs, Thu Feb  5 17:11:06 PST 2004
//    Moved inlined destructor definition to .C file because certain compilers
//    have problems with them.
//
//    Hank Childs, Wed Aug 11 09:47:46 PDT 2004
//    Added ModifyContract.
//
//    Sean Ahern, Wed Sep 10 13:13:18 EDT 2008
//    For ease of reading code, I forced the argument to be avtCentering, not
//    integer.
//
//    Eric Brugger, Tue Jul 22 08:00:59 PDT 2014
//    Modified the class to work with avtDataRepresentation.
//
// ****************************************************************************

class AVTFILTERS_API avtShiftCenteringFilter : public avtDataTreeIterator
{
  public:
                            avtShiftCenteringFilter(avtCentering);
    virtual                ~avtShiftCenteringFilter();

    virtual const char     *GetType(void) { return "avtShiftCenteringFilter"; };
    virtual const char     *GetDescription(void) 
                                  { return "Re-centering data"; };

  protected:
    avtCentering            centeringTarget;

    virtual avtDataRepresentation *ExecuteData(avtDataRepresentation *);
    virtual void            UpdateDataObjectInfo(void);
    virtual avtContract_p
                            ModifyContract(avtContract_p);
};


#endif


