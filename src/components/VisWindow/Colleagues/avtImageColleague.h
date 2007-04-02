/*****************************************************************************
*
* Copyright (c) 2000 - 2006, The Regents of the University of California
* Produced at the Lawrence Livermore National Laboratory
* All rights reserved.
*
* This file is part of VisIt. For details, see http://www.llnl.gov/visit/. The
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
*    documentation and/or materials provided with the distribution.
*  - Neither the name of the UC/LLNL nor  the names of its contributors may be
*    used to  endorse or  promote products derived from  this software without
*    specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT  HOLDERS AND CONTRIBUTORS "AS IS"
* AND ANY EXPRESS OR  IMPLIED WARRANTIES, INCLUDING,  BUT NOT  LIMITED TO, THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND  FITNESS FOR A PARTICULAR  PURPOSE
* ARE  DISCLAIMED.  IN  NO  EVENT  SHALL  THE  REGENTS  OF  THE  UNIVERSITY OF
* CALIFORNIA, THE U.S.  DEPARTMENT  OF  ENERGY OR CONTRIBUTORS BE  LIABLE  FOR
* ANY  DIRECT,  INDIRECT,  INCIDENTAL,  SPECIAL,  EXEMPLARY,  OR CONSEQUENTIAL
* DAMAGES (INCLUDING, BUT NOT  LIMITED TO, PROCUREMENT OF  SUBSTITUTE GOODS OR
* SERVICES; LOSS OF  USE, DATA, OR PROFITS; OR  BUSINESS INTERRUPTION) HOWEVER
* CAUSED  AND  ON  ANY  THEORY  OF  LIABILITY,  WHETHER  IN  CONTRACT,  STRICT
* LIABILITY, OR TORT  (INCLUDING NEGLIGENCE OR OTHERWISE)  ARISING IN ANY  WAY
* OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
* DAMAGE.
*
*****************************************************************************/

#ifndef AVT_IMAGE_COLLEAGUE_H
#define AVT_IMAGE_COLLEAGUE_H

#include <string>

#include <ColorAttribute.h>
#include <avtAnnotationColleague.h>
#include <viswindow_exports.h>

class vtkActor2D;
class vtkImageData;
class vtkImageMapper;
class vtkImageReader2;
class vtkImageResample;

// ****************************************************************************
// Class: avtImageColleague
//
// Purpose:
//   This colleague is a image that can be shown in the vis window.
//
// Notes:      
//
// Programmer: John C. Anderson
// Creation:   Thu Jul 15 08:04:46 PDT 2004
//
// Modifications:
//   
// ****************************************************************************

class VISWINDOW_API avtImageColleague : public avtAnnotationColleague
{
public:
    avtImageColleague(VisWindowColleagueProxy &);
    virtual ~avtImageColleague();

    virtual void AddToRenderer();
    virtual void RemoveFromRenderer();
    virtual void Hide();

    // Methods to set and get the annotation's properties.
    virtual void SetOptions(const AnnotationObject &annot);
    virtual void GetOptions(AnnotationObject &annot);

    // Methods that are called in response to vis window events.
    virtual void HasPlots(void);
    virtual void NoPlots(void);

protected:
    void CreateActorAndMapper();
    bool UpdateImage(std::string);

    vtkActor2D                 *actor;
    vtkImageMapper             *mapper;
    vtkImageResample           *resample;

    std::string                 currentImage;
    vtkImageData               *iData;

    int                         width, height;

    bool                        useOpacityColor;
    ColorAttribute              opacityColor;

    bool                        maintainAspectRatio;

    bool                        addedToRenderer;

    bool ShouldBeAddedToRenderer() const;
};


#endif


