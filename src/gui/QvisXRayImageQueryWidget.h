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

#ifndef QVIS_XRAYIMAGE_QUERY_WIDGET_H
#define QVIS_XRAYIMAGE_QUERY_WIDGET_H
#include <gui_exports.h>
#include <QWidget>
#include <vectortypes.h>

// Forward declarations.
class QCheckBox;
class QComboBox;
class QLineEdit;
class MapNode;

// ****************************************************************************
// Class: QvisXRayImageQueryWidget
//
// Purpose:
//   This widget provides options for performing a XRayImage query.
//
// Notes:      
//
// Programmer: Kathleen Biagas
// Creation:   June 17, 2011 
//
// Modifications:
//   Kathleen Biagas, Wed Oct 17 12:12:10 PDT 2012
//   Added upVector.
//
//   Eric Brugger, Fri May 22 15:50:50 PDT 2015
//   I updated the window to use the new view description and support the
//   recently added background intensity parameter.
//
//   Eric Brugger, Wed May 27 17:27:31 PDT 2015
//   I added an option to family output files.
//
//   Eric Brugger, Thu Jun  4 17:23:58 PDT 2015
//   I added an option to enable outputting the ray bounds to a vtk file.
//
// ****************************************************************************

class GUI_API QvisXRayImageQueryWidget : public QWidget
{
    Q_OBJECT
public:
    QvisXRayImageQueryWidget(QWidget *parent = 0, Qt::WindowFlags f = 0);
    virtual ~QvisXRayImageQueryWidget();

    bool GetQueryParameters(MapNode &params);


private:
    bool             GetDoubleValues(int whichWidget, doubleVector &pt);
    bool             GetDoubleValues(int whichWidget, int n, double *pt);
    bool             GetIntValues(int whichWidget, int *pt);

    QComboBox       *imageFormat;
    QCheckBox       *divideFlag;
    QLineEdit       *backgroundIntensities;
    QLineEdit       *normal;
    QLineEdit       *focus;
    QLineEdit       *viewUp;
    QLineEdit       *viewAngle;
    QLineEdit       *parallelScale;
    QLineEdit       *nearPlane;
    QLineEdit       *farPlane;
    QLineEdit       *imagePan;
    QLineEdit       *imageZoom;
    QCheckBox       *perspective;
    QCheckBox       *family;
    QCheckBox       *outputRayBounds;
    QLineEdit       *imageSize;
};

#endif
