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

#ifndef QVIS_GAUSSIAN_OPACITY_BAR_H
#define QVIS_GAUSSIAN_OPACITY_BAR_H
#include <gui_exports.h>

#include <QvisAbstractOpacityBar.h>

class QPixmap;

// ****************************************************************************
//  Class:  QvisGaussianOpacityBar
//
//  Purpose:
//    Gaussian-max implementation of QvisAbstractOpacityBar
//
//  Programmer:  Jeremy Meredith
//  Creation:    January 30, 2001
//
//  Modifications:
//    Jeremy Meredith, Mon Feb  7 10:37:27 PST 2005
//    Removed mouseReleased because it was already in the base class.
//
//    Brad Whitlock, Wed Jun  4 10:32:52 PDT 2008
//    Qt 4.
//
//    Brad Whitlock, Thu Dec 18 14:07:00 PST 2008
//    I changed the drawOpacities method.
//
// ****************************************************************************

class GUI_API QvisGaussianOpacityBar : public QvisAbstractOpacityBar
{
    Q_OBJECT
  public:
                  QvisGaussianOpacityBar(QWidget *parent=NULL);
                 ~QvisGaussianOpacityBar();
    float        *getRawOpacities(int);
    int           getNumberOfGaussians();
    void          getGaussian(int, float*,float*,float*,float*,float*);
    void          setAllGaussians(int, float*);

  protected:
    virtual void  mouseMoveEvent(QMouseEvent*);
    virtual void  mousePressEvent(QMouseEvent*);
    virtual void  mouseReleaseEvent(QMouseEvent*);
    virtual void  drawOpacities();
    void          drawControlPoints();

  private:
    enum Mode     {modeNone, modeX, modeH, modeW, modeWR, modeWL, modeB};
    // encapsulation of gaussian parameters
    class Gaussian
    {
      public:
        float x;
        float h;
        float w;
        float bx;
        float by;
      public:
        Gaussian(float x_,float h_,float w_,float bx_,float by_) : x(x_),h(h_),w(w_),bx(bx_),by(by_) {};
        Gaussian() {};
        ~Gaussian() {};
    };

    // the list of gaussians
    int         ngaussian;
    Gaussian    gaussian[200];

    // the current interaction mode and the current gaussian
    Mode        currentMode;
    int         currentGaussian;

    // GUI interaction variables
    bool        mousedown;
    int         lastx;
    int         lasty;

    // helper functions
    bool findGaussianControlPoint(int,int, int*,Mode*);
    void removeGaussian(int);
    void addGaussian(float,float,float,float,float);
};

#endif
