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
//                         VisWindowToolProxy.h                              //
// ************************************************************************* //

#ifndef VIS_WINDOW_TOOL_PROXY_H
#define VIS_WINDOW_TOOL_PROXY_H
#include <viswindow_exports.h>
#include <VisWindowProtectionProxy.h>

// ****************************************************************************
// Class: VisWindowToolProxy
//
// Purpose:
//   This is a proxy specific to VisWindow tools.  They may use it
//   to access the VisWindow.  This object is not a friend to VisWindow,
//   but it may access it through its base class' methods and friend 
//   status.
//
// Programmer: Brad Whitlock
// Creation:   Tue Oct 2 16:52:41 PST 2001
//
// Modifications:
//
//   Hank Childs, Wed Sep  4 08:40:47 PDT 2002
//   Added RecalculateRenderOrder.
//
//   Kathleen Bonnell, Fri Jun  6 15:31:37 PDT 2003  
//   Added GetFullFrameMode and GetScaleFactorAndType. 
//
//   Kathleen Bonnell, Wed Dec  3 16:49:14 PST 2003
//   Added TransparenciesExist.
//
//   Mark Blair, Wed Aug 30 14:09:00 PDT 2006
//   Added GetSize, GetWindowSize, GetPlotListIndex.
//
//   Mark Blair, Wed Oct 25 15:12:55 PDT 2006
//   Added GetPlotInfoAtts.
//
//   Brad Whitlock, Wed Jan  7 14:51:13 PST 2009
//   Removed GetPlotListIndex, GetPlotInfoAtts.
//
//   Brad Whitlock, Fri Apr  3 10:54:57 PDT 2009
//   I added GetBoundingBoxMode.
//
//   Jeremy Meredith, Tue Feb  2 11:12:10 EST 2010
//   Added GetToolUpdateMode.
//
//   Jeremy Meredith, Wed May 19 14:15:58 EDT 2010
//   Support 3D axis scaling (3D equivalent of full-frame mode).
//
// ****************************************************************************

class VISWINDOW_API VisWindowToolProxy : public VisWindowProtectionProxy
{
  public:
                        VisWindowToolProxy(VisWindow *vw) 
                             : VisWindowProtectionProxy(vw) {;};

    void                DisableUpdates()
                                  { ProxiedDisableUpdates(); };

    void                EnableUpdates()
                                  { ProxiedEnableUpdates(); };

    void                EndBoundingBox()
                                  { ProxiedEndBoundingBox(); };

    void                GetForegroundColor(double *fg)
                                  { ProxiedGetForegroundColor(fg); };

    vtkRenderer        *GetBackground()
                                  { return ProxiedGetBackground(); };

    bool                GetBoundingBoxMode() const { return ProxiedGetBoundingBoxMode(); }

    void                GetBounds(double bounds[6])
                                  { ProxiedGetBounds(bounds); };

    vtkRenderer        *GetCanvas()
                                  { return ProxiedGetCanvas(); };

    vtkRenderer        *GetForeground()
                                  { return ProxiedGetForeground(); };

    WINDOW_MODE         GetMode() { return ProxiedGetMode(); };

    void                GetSize(int &width, int &height) const
                                  { ProxiedGetSize(width, height); };

    void                GetWindowSize(int &width, int &height) const
                                  { ProxiedGetWindowSize(width, height); };

    bool                HasPlots(){ return ProxiedHasPlots(); };
    
    void                Render()  { ProxiedRender(); };

    void                RecalculateRenderOrder()
                                  { ProxiedRecalculateRenderOrder(); };

    void                StartBoundingBox()
                                  { ProxiedStartBoundingBox(); };

    bool                UpdatesEnabled()
                                  { return ProxiedUpdatesEnabled(); };

    bool                GetFullFrameMode()
                                  { return ProxiedGetFullFrameMode(); };

    void                GetScaleFactorAndType(double &scale, int &type)
                            { ProxiedGetScaleFactorAndType(scale, type); };

    bool                Get3DAxisScalingFactors(double s[3]) 
                                  { return ProxiedGet3DAxisScalingFactors(s); }

    bool                TransparenciesExist()
                                  { return ProxiedTransparenciesExist(); };

    TOOLUPDATE_MODE     GetToolUpdateMode() const { return ProxiedGetToolUpdateMode(); }
};


#endif


