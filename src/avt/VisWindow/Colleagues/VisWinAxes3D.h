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
//                            VisWinAxes3D.h                                 //
// ************************************************************************* //

#ifndef VIS_WIN_AXES3D_H
#define VIS_WIN_AXES3D_H
#include <viswindow_exports.h>


#include <VisWinColleague.h>


class vtkVisItCubeAxesActor;
class vtkPolyDataMapper;
class vtkOutlineSource;

// ****************************************************************************
//  Class: VisWinAxes3D
//
//  Purpose:
//      This is a concrete colleague for the mediator VisWindow.  It places
//      a 3D axes around the bounding box for the plots.
//
//  Programmer: Kathleen Bonnell 
//  Creation:   June 28, 2001
//
//  Modifications:
//
//    Kathleen Bonnell, Fri Aug  3 14:55:59 PDT 2001
//    Changed from using a 2d cube axes actor to using a 3d version.
//
//    Kathleen Bonnell, Tue Oct 30 10:30:10 PST 2001 
//    Moved AdjustValues, AdjustRange and related members to the
//    more appropriate location of vtkVisItCubeAxesActor
//
//    Kathleen Bonnell, Wed May 28 15:52:32 PDT 2003   
//    Added method 'ReAddToWindow'. 
//    
//    Kathleen Bonnell, Tue Dec 16 11:34:33 PST 2003 
//    Added method 'SetLabelScaling'.
//
//    Brad Whitlock, Thu Jul 28 10:10:40 PDT 2005
//    Added methods to set the units and title.
//
//    Brad Whitlock, Tue Mar 25 16:27:40 PDT 2008
//    Added methods for line width, font, font size.
// 
//    Eric Brugger, Wed Oct 15 13:05:33 PDT 2008
//    Added SetAutoSetTicks, SetMajorTickMinimum, SetMajorTickMaximum,
//    SetMajorTickSpacing and SetMinorTickSpacing.
//
//    Jeremy Meredith, Wed May  5 14:32:23 EDT 2010
//    Added support for title visibility separate from label visibility.
//
//    Jeremy Meredith, Wed May 19 14:15:58 EDT 2010
//    Account for 3D axis scaling (3D equivalent of full-frame mode).
//
//    Hank Childs, Mon May 23 13:26:09 PDT 2011
//    Add method for overriding bounding box location.
//
// ****************************************************************************

class VISWINDOW_API VisWinAxes3D : public VisWinColleague
{
  public:
                              VisWinAxes3D(VisWindowColleagueProxy &);
    virtual                  ~VisWinAxes3D();

    virtual void              SetForegroundColor(double, double, double);
    virtual void              UpdateView(void);
    virtual void              UpdatePlotList(std::vector<avtActor_p> &);

    virtual void              Start3DMode(void);
    virtual void              Stop3DMode(void);

    virtual void              HasPlots(void);
    virtual void              NoPlots(void);
    virtual void              ReAddToWindow(void);

    void                      SetBBoxLocation(bool, const double *);
    void                      SetBounds(double [6], double scales[3]);
    void                      SetXTickVisibility(int, int);
    void                      SetYTickVisibility(int, int);
    void                      SetZTickVisibility(int, int);
    void                      SetXLabelVisibility(int);
    void                      SetYLabelVisibility(int);
    void                      SetZLabelVisibility(int);
    void                      SetXTitleVisibility(int);
    void                      SetYTitleVisibility(int);
    void                      SetZTitleVisibility(int);
    void                      SetXGridVisibility(int);
    void                      SetYGridVisibility(int);
    void                      SetZGridVisibility(int);
    void                      SetVisibility(int);
    void                      SetAutoSetTicks(int);
    void                      SetMajorTickMinimum(double, double, double);
    void                      SetMajorTickMaximum(double, double, double);
    void                      SetMajorTickSpacing(double, double, double);
    void                      SetMinorTickSpacing(double, double, double);
    void                      SetBBoxVisibility(int);
    void                      SetFlyMode(int);
    void                      SetTickLocation(int);
    void                      SetLabelScaling(bool, int, int, int);

    void                      SetXTitle(const std::string &, bool);
    void                      SetXUnits(const std::string &, bool);
    void                      SetYTitle(const std::string &, bool);
    void                      SetYUnits(const std::string &, bool);
    void                      SetZTitle(const std::string &, bool);
    void                      SetZUnits(const std::string &, bool);

    void                      SetLineWidth(int);
    void                      SetTitleTextAttributes(
                                  const VisWinTextAttributes &xAxis, 
                                  const VisWinTextAttributes &yAxis,
                                  const VisWinTextAttributes &zAxis);
    void                      SetLabelTextAttributes(
                                  const VisWinTextAttributes &xAxis, 
                                  const VisWinTextAttributes &yAxis,
                                  const VisWinTextAttributes &zAxis);

    void                      Set3DAxisScalingFactors(bool scale,
                                                       const double s[3]);
  protected:
    void UpdateTitleTextAttributes(double fr, double fg, double fb);
    void UpdateLabelTextAttributes(double fr, double fg, double fb);

    vtkVisItCubeAxesActor    *axes;
    vtkOutlineSource         *axesBoxSource;
    vtkPolyDataMapper        *axesBoxMapper;
    vtkActor                 *axesBox;

    bool                      addedAxes3D;
    double                    currentBounds[6];
    double                    currentScaleFactors[3];
    bool                      visibility;

    bool                      boundsOverridden;
    double                    overrideBounds[6];

    std::string               userXTitle;
    std::string               userYTitle;
    std::string               userZTitle;
    std::string               userXUnits;
    std::string               userYUnits;
    std::string               userZUnits;
    bool                      userXTitleFlag;
    bool                      userYTitleFlag;
    bool                      userZTitleFlag;
    bool                      userXUnitsFlag;
    bool                      userYUnitsFlag;
    bool                      userZUnitsFlag;

    VisWinTextAttributes      titleTextAttributes[3];
    VisWinTextAttributes      labelTextAttributes[3];

    void                      AddAxes3DToWindow(void);
    void                      RemoveAxes3DFromWindow(void);
    bool                      ShouldAddAxes3D(void);
};


#endif


