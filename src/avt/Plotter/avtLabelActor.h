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
//                              avtLabelActor.h                              //
// ************************************************************************* //


#ifndef AVT_LABEL_ACTOR_H
#define AVT_LABEL_ACTOR_H
#include <plotter_exports.h>
#include <visitstream.h>
#include <ref_ptr.h>

class vtkFollower;
class vtkPolyDataMapper;
class vtkRenderer;


// ****************************************************************************
//  Class:  avtLabelActor
//
//  Purpose:  Responsible for creating label actors used in decorations. 
//
//  Programmer:  Kathleen Bonnell 
//  Creation:    July 12, 2002 
//
//  Modifications:
//    Kathleen Bonnell, Fri Jul 19 08:39:04 PDT 2002
//    Add ComputeScaleFactor.
//
//    Eric Brugger, Tue Dec  9 16:16:49 PST 2008
//    Added the ability to display a marker instead of a text string.
//
//    Eric Brugger, Mon Mar  9 17:26:13 PDT 2009
//    Added an overloaded version of SetForegroundColor that allows setting
//    rgba instead of just rgb.
//
//    Eric Brugger, Thu Feb 19 13:24:04 PST 2013
//    Added the ability to set a scale factor and the line width. 
//
// ****************************************************************************

class PLOTTER_API avtLabelActor
{
  public:
                       avtLabelActor();
    virtual           ~avtLabelActor();

    void               Add(vtkRenderer *ren);
    void               Remove();
    void               Hide();
    void               UnHide();

    void               SetAttachmentPoint(const double newPos[3]);
    const double *     GetAttachmentPoint() { return attach; };
    void               SetScale(double);
    void               SetScaleFactor(double);
    void               SetLineWidth(int);
    void               SetDesignator(const char *l);
    void               SetMarker(const int index);
    void               SetForegroundColor(double fgr, double fgg, double fgb);
    void               SetForegroundColor(double fgr, double fgg, double fgb,
                           double fga);
    void               SetForegroundColor(double fg[3]);
    void               Shift(const double vec[3]);
    double             ComputeScaleFactor();

  protected:
    double             attach[3];
    double             scaleFactor;
    vtkFollower       *labelActor;

    vtkRenderer       *renderer; 

  private:
};

typedef ref_ptr<avtLabelActor> avtLabelActor_p;

#endif
