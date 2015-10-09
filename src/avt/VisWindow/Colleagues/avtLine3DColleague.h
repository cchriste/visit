/*****************************************************************************
*
* Copyright (c) 2000 - 2015, Lawrence Livermore National Security, LLC
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

#ifndef AVT_LINE3D_COLLEAGUE_H
#define AVT_LINE3D_COLLEAGUE_H

#include <ColorAttribute.h>
#include <avtAnnotationColleague.h>
#include <viswindow_exports.h>

class vtkActor;
class vtkConeSource;
class vtkLineSource;
class vtkPolyDataMapper;
class vtkTubeFilter;

// ****************************************************************************
// Class: avtLine3DColleague
//
// Purpose:
//   This colleague is a 3D line that can be shown in the vis window.
//
// Notes:      
//
// Programmer: Kathleen Biagas 
// Creation:   July 13, 2015 
//
// Modifications:
//   Kathleen Biagas, Tue Jul 14 16:35:47 PDT 2015
//   Add support for arrows and tube.
//
// ****************************************************************************

class VISWINDOW_API avtLine3DColleague : public avtAnnotationColleague
{
public:
    avtLine3DColleague(VisWindowColleagueProxy &);
    virtual ~avtLine3DColleague();

    virtual void AddToRenderer();
    virtual void RemoveFromRenderer();
    virtual void Hide();

    virtual std::string TypeName() const { return "Line3D"; }

    // Methods to set and get the annotation's properties.
    virtual void SetOptions(const AnnotationObject &annot);
    virtual void GetOptions(AnnotationObject &annot);

    // Methods that are called in response to vis window events.
    virtual void SetForegroundColor(double r, double g, double b);
    virtual void HasPlots(void);
    virtual void NoPlots(void);

protected:
    vtkLineSource       *lineSource;
    vtkConeSource       *arrow1Source;
    vtkConeSource       *arrow2Source;
    vtkPolyDataMapper   *lineMapper;
    vtkPolyDataMapper   *arrow1Mapper;
    vtkPolyDataMapper   *arrow2Mapper;
    vtkActor            *lineActor;
    vtkActor            *arrow1Actor;
    vtkActor            *arrow2Actor;
    vtkTubeFilter       *tubeFilter;

    bool                 addedToRenderer;
    bool                 useForegroundForLineColor;
    bool                 useArrow1;
    bool                 useArrow2;
    bool                 arrow1Added;
    bool                 arrow2Added;
    int                  lineType;
    ColorAttribute       lineColor;

    bool ShouldBeAddedToRenderer() const;

};


#endif


