/*****************************************************************************
*
* Copyright (c) 2000 - 2010, Lawrence Livermore National Security, LLC
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
//                                avtCurvePlot.h                             //
// ************************************************************************* //

#ifndef AVT_CURVE_PLOT_H
#define AVT_CURVE_PLOT_H


#include <avtLegend.h>
#include <avtPlot.h>
#include <avtCurveRenderer.h>

#include <CurveAttributes.h>

class     avtCurveFilter;
class     avtCurveLegend;
class     avtLabeledCurveMapper;
class     avtUserDefinedMapper;
class     avtWarpFilter;
class     vtkProperty;


// ****************************************************************************
//  Class:  avtCurvePlot
//
//  Purpose:
//      A concrete type of avtPlot, this is the Curve plot.
//
//  Programmer: kbonnell -- generated by xml2info
//  Creation:   Sat Apr 20 13:01:58 PST 2002
//
//  Modifications:
//    Kathleen Bonnell, Fri Jul 12 16:53:11 PDT 2002 
//    Added labeled curve mapper.
//
//    Kathleen Bonnell, Tue Oct 22 08:33:26 PDT 2002
//    Added ApplyRenderingTransformation. 
//    
//    Kathleen Bonnell, Thu Oct 27 15:12:13 PDT 2005 
//    Added Legend. 
//    
//    Kathleen Bonnell, Wed Jul 12 08:30:04 PDT 2006 
//    Added warp filter. 
//    
//    Kathleen Bonnell, Tue Apr  3 17:17:33 PDT 2007 
//    Added CanDoCurveViewScaling. 
//    
//    Kathleen Bonnell, Tue Mar  3 13:37:13 PST 2009
//    Removed CanDo2DViewScaling (moved into Viewer PluginInfo)
//
// ****************************************************************************


class avtCurvePlot : public avtLineDataPlot
{
  public:
                                avtCurvePlot();
    virtual                    ~avtCurvePlot();

    virtual const char         *GetName(void) { return "CurvePlot"; };
    virtual void                ReleaseData(void);

    static avtPlot             *Create();

    virtual void                SetAtts(const AttributeGroup*);
    void                        SetLineWidth(int);
    void                        SetLineStyle(int);

  protected:
    CurveAttributes                atts;

    avtCurveLegend                *curveLegend;
    avtLegend_p                    curveLegendRefPtr;

    avtCurveRenderer_p              renderer;
    avtUserDefinedMapper           *mapper;
    avtLabeledCurveMapper          *decoMapper;
    vtkProperty                    *property;

    avtCurveFilter                 *CurveFilter;
    avtWarpFilter                  *WarpFilter;

    virtual avtMapper          *GetMapper(void);
    virtual avtDataObject_p     ApplyOperators(avtDataObject_p);
    virtual avtDataObject_p     ApplyRenderingTransformation(avtDataObject_p);
    virtual void                CustomizeBehavior(void);

    virtual avtLegend_p         GetLegend(void) { return curveLegendRefPtr; };
    virtual avtDecorationsMapper *GetDecorationsMapper(void); 
};

#endif


