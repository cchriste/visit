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

// ************************************************************************* //
//                       VisWinRenderingWithoutWindow.h                      //
// ************************************************************************* //

#ifndef VIS_WIN_RENDERING_WITHOUT_WINDOW_H
#define VIS_WIN_RENDERING_WITHOUT_WINDOW_H
#include <viswindow_exports.h>
#include <VisWinRendering.h>


class     vtkRenderWindow;


// ****************************************************************************
//  Class: VisWinRenderingWithoutWindow
//
//  Purpose:
//      A derived type of VisWinRendering that assumes that there will be no
//      window on the screen.
//
//  Programmer: Hank Childs
//  Creation:   January 29, 2002
//
//  Modifications:
//      Sean Ahern, Mon May 20 13:35:16 PDT 2002
//      Added empty functions for Raise/Lower.
//
//      Kathleen Bonnell, Tue Feb 11 11:28:03 PST 2003  
//      Removed member 'iren'.  Made GetRenderWindowInteractor return NULL. 
//
// ****************************************************************************

class VISWINDOW_API VisWinRenderingWithoutWindow : public VisWinRendering
{
  public:
                                       VisWinRenderingWithoutWindow(
                                                    VisWindowColleagueProxy &);
    virtual                           ~VisWinRenderingWithoutWindow();

    virtual void                       Iconify(void) {;};
    virtual void                       DeIconify(void) {;};
    virtual void                       Show(void) {;};
    virtual void                       Hide(void) {;};
    virtual void                       Raise(void) {;};
    virtual void                       Lower(void) {;};

    virtual void                       SetResizeEvent(void(*callback)(void *), void *){;};
    virtual void                       SetCloseCallback(void(*callback)(void *), void *)
                                                                           {;};

  protected:
    vtkRenderWindow                   *renWin;

    virtual vtkRenderWindow           *GetRenderWindow(void);
    virtual vtkRenderWindowInteractor *GetRenderWindowInteractor(void)
                                           { return NULL; };

    virtual void                       RealizeRenderWindow(void) {;};
};


#endif


