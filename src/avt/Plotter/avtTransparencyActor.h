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
//                          avtTransparencyActor.h                           //
// ************************************************************************* //

#ifndef AVT_TRANSPARENCY_ACTOR_H
#define AVT_TRANSPARENCY_ACTOR_H

#include <plotter_exports.h>

#include <vector>
#include <map>
#include <cstring>

class     vtkActor;
class     vtkAppendPolyData;
class     vtkAxisDepthSort;
class     vtkCamera;
class     vtkDataSet;
class     vtkDataSetMapper;
class     vtkDepthSortPolyData;
class     vtkMatrix4x4;
class     vtkPolyData;
class     vtkPolyDataMapper;
class     vtkRenderer;
class     vtkParallelImageSpaceRedistributor;

class     ColorAttribute;


// ****************************************************************************
//  Class: avtTransparencyActor
//
//  Purpose:
//      Rendering transparent polygons correctly requires sorting them in front
//      to back order.  Since the AVT model is to keep groups of polygons
//      separated, this sorting would be impossible.  This class collects the
//      inputs from all of the plots and keeps them together so they can be
//      sorted.
//
//  Programmer: Hank Childs
//  Creation:   July 3, 2002
//
//  Modifications:
//
//    Hank Childs, Thu Jul 11 16:02:45 PDT 2002
//    Added the concept of visibility.  This is "different" from useActor
//    because the VisWindow has both concepts and they don't always agree.  At
//    this point it was easier to adopt its concept than re-write that code.
//
//    Hank Childs, Sun Jul 14 15:49:58 PDT 2002
//    Make use of one VTK module to do all the sorting.
//
//    Jeremy Meredith, Fri Jul 26 14:28:09 PDT 2002
//    Made perfect sorting a permanent mode instead of a one-frame mode.
//
//    Kathleen Bonnell, Wed Jul 16 16:39:02 PDT 2003
//    Added ScaleByVector method. 
//
//    Kathleen Bonnell, Wed Dec  3 16:42:38 PST 2003 
//    Added TransparenciesExist method.
//
//    Chris Wojtan, Fri Jun 25 15:15 PDT 2004
//    Added is2Dimensional bool and a function to set and retrieve its value
//
//    Hank Childs, Wed Sep  8 17:55:34 PDT 2004
//    No longer inline is2Dimensional.
//
//    Chris Wojtan, Wed Jul 7 10:27 PDT 2004
//    Added MyParallelFilter for improved parallel transparency
//
//    Jeremy Meredith, Thu Oct 21 12:09:00 PDT 2004
//    Renamed the parallel filter.
//
//    Tom Fogal, Sun May 24 19:34:25 MDT 2009
//    Added const to method.
//
//    Tom Fogal, Mon May 25 14:03:14 MDT 2009
//    Added caching for transparency calculation.
//
//    Hank Childs, Wed Feb 17 18:19:53 CST 2010
//    Add support for changes in specular lighting.
//
// ****************************************************************************

class PLOTTER_API avtTransparencyActor
{
  public:
                                     avtTransparencyActor();
    virtual                         ~avtTransparencyActor();

    int                              AddInput(std::vector<vtkDataSet *> &,
                                              std::vector<vtkDataSetMapper *>&,
                                              std::vector<vtkActor *> &);
    void                             ReplaceInput(int,
                                              std::vector<vtkDataSet *> &,
                                              std::vector<vtkDataSetMapper *>&,
                                              std::vector<vtkActor *> &);

    void                             TurnOffInput(int);
    void                             TurnOnInput(int);
    void                             RemoveInput(int);
    void                             InputWasModified(int, double=-1.0);
    void                             SetVisibility(int, bool);
    void                             VisibilityOff(void);
    void                             VisibilityOn(void);

    void                             PrepareForRender(vtkCamera *);
    bool                             UsePerfectSort(bool);

    bool                             TransparenciesExist(void);
    bool                             TransparenciesMightExist(void) const;
    void                             InvalidateTransparencyCache() {
                                         cachedTransparencies = false;
                                     }

    void                             AddToRenderer(vtkRenderer *);
    void                             RemoveFromRenderer(vtkRenderer *);
    void                             ScaleByVector(const double vec[3]);
    void                             SetSpecularProperties(bool flag,
                                           double coeff, double power,
                                           const ColorAttribute &color);

    void                             SuspendRendering() { renderingSuspended = true;  }
    void                             ResumeRendering()  { renderingSuspended = false; }

    bool                             GetIs2Dimensional() const
                                                    { return is2Dimensional; };
    void                             SetIs2Dimensional(bool val);

  protected:
    std::vector<std::vector <vtkDataSet *> >         datasets;
    std::vector<std::vector <vtkDataSetMapper *> >   mappers;
    std::vector<std::vector <vtkActor *> >           actors;

    std::map<int,double>                             inputsOpacities;

    std::vector<std::vector <vtkPolyData *> >        preparedDataset;

    std::vector<bool>                                useActor;
    std::vector<bool>                                visibility;
    std::vector<bool>                                lastExecutionActorList;
    bool                                             inputModified;

    vtkAppendPolyData                               *appender;
    vtkParallelImageSpaceRedistributor              *parallelFilter;
    vtkActor                                        *myActor;
    vtkPolyDataMapper                               *myMapper;

    vtkAxisDepthSort                                *axisSort;
    vtkDepthSortPolyData                            *perfectSort;
    bool                                             usePerfectSort;
    bool                                             is2Dimensional;
    vtkMatrix4x4                                    *lastCamera;

    bool                                             renderingSuspended;
    bool                                             transparenciesExist;
    bool                                             cachedTransparencies;

    void                                             SetUpActor(void);
    void                                             PrepareDataset(size_t, size_t);
    void                                             DetermineTransparencies();
};
#endif
