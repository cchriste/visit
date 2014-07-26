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

#ifndef VIEW_STACK_H
#define VIEW_STACK_H
#include <viewer_exports.h>
#include <avtViewCurve.h>
#include <avtView2D.h>
#include <avtView3D.h>
#include <avtViewAxisArray.h>

#define VSTACK_SIZE 15

// ****************************************************************************
// Class: ViewStack
//
// Purpose:
//   Contains the stacks that allow us to stack different view types.
//
// Notes:      
//
// Programmer: Brad Whitlock
// Creation:   Tue Mar 7 17:24:46 PST 2006
//
// Modifications:
//    Jeremy Meredith, Thu Jan 31 14:56:06 EST 2008
//    Added new axis array window mode.
//   
// ****************************************************************************

class VIEWER_API ViewStack
{
public:
    ViewStack();
    ViewStack(bool);
    ViewStack(const ViewStack &);
    ~ViewStack();

    void Clear();

    bool PopViewCurve(avtViewCurve &);
    bool PopView2D(avtView2D &);
    bool PopView3D(avtView3D &);
    bool PopViewAxisArray(avtViewAxisArray &);

    void PushViewCurve(const avtViewCurve &);
    void PushView2D(const avtView2D &);
    void PushView3D(const avtView3D &);
    void PushViewAxisArray(const avtViewAxisArray &);

    bool HasViewCurves() const;
    bool HasView2Ds() const;
    bool HasView3Ds() const;
    bool HasViewAxisArrays() const;

    // Assignment operator.
    void operator = (const ViewStack &);
private:
    bool         preventPopFirst;

    avtViewCurve     viewCurveStack[VSTACK_SIZE];
    int              viewCurveStackTop;
    avtView2D        view2DStack[VSTACK_SIZE];
    int              view2DStackTop;
    avtView3D        view3DStack[VSTACK_SIZE];
    int              view3DStackTop;
    avtViewAxisArray viewAxisArrayStack[VSTACK_SIZE];
    int              viewAxisArrayStackTop;
};

#endif
