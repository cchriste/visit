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
//                        avtLabeledCurveMapper.h                            //
// ************************************************************************* //

#ifndef AVT_LABELEDCURVE_MAPPER_H
#define AVT_LABELEDCURVE_MAPPER_H
#include <plotter_exports.h>


#include <avtDecorationsMapper.h>

class     vtkMaskPoints;
class     vtkDataSet;


// ****************************************************************************
//  Class: avtLabeledCurveMapper
//
//  Purpose:
//      A mapper for labeled curves.  This extends the functionality of a 
//      decorations mapper by mapping a label-glyph onto a curve dataset.
//
//  Programmer: Kathleen Bonnell 
//  Creation:   July 12, 2002 
//
//  Modifications:
//    Kathleen Bonnell, Tue Jul 23 15:01:55 PDT 2002 
//    Added member to keep track of label visibility.
// 
//    Kathleen Bonnell, Thu Feb 17 09:09:41 PST 2011
//    Add ability to set multiple labels, keep track of which input an actor
//    is associated with.
//
// ****************************************************************************

class PLOTTER_API  avtLabeledCurveMapper : public avtDecorationsMapper
{
  public:
                               avtLabeledCurveMapper();
    virtual                   ~avtLabeledCurveMapper();

    void                       SetLabelColor(double [3]);
    void                       SetLabelColor(double, double, double);
    void                       SetLabel(std::string &);
    void                       SetLabels(std::vector<std::string> &);
    void                       SetScale(double);
    void                       SetLabelVisibility(bool);

  protected:
    std::string                label; 
    std::vector<std::string>   labels; 
    std::vector<int>           actorsInputNum; 
    bool                       labelVis;
    double                     labelColor[3];
    double                     scale;

    vtkMaskPoints            **filter;
    int                        nFilters;

    virtual void               CustomizeMappers(void);

    virtual void               SetDatasetInput(vtkDataSet *, int);
    virtual void               SetUpFilters(int);
};


#endif


