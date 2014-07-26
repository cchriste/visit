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

#ifndef QVIS_CLIP_WINDOW_H
#define QVIS_CLIP_WINDOW_H
#include <QvisOperatorWindow.h>
#include <QGroupBox>

class QCheckBox;
class QComboBox;
class QButtonGroup;
class QGrid;
class QLineEdit;
class QLabel;
class QGroupBox;
class QVBoxLayout;
class QPlaneGroup;
class ClipAttributes;


// ****************************************************************************
// Class: QPlaneGroup
//
// Purpose: 
//   Widget that encapsualtes the options for a single clipping plane.
//
// Programmer: Cyrus Harrison
// Creation:   Thu Aug 21 13:45:03 PDT 2008
//
// Modifications:
//
// ****************************************************************************
class QPlaneGroup: public QGroupBox
{
Q_OBJECT
public:
    QPlaneGroup(const QString &title,QWidget *parent=0);
    
    virtual ~QPlaneGroup();
    
    void SetOrigin(double val[3]);
    void SetNormal(double val[3]);
    bool GetOrigin(double val[3]);
    bool GetNormal(double val[3]);

signals:
    void OriginChanged();
    void NormalChanged();
    
private:
    QLineEdit *origin;
    QLineEdit *normal;
    
};


// ****************************************************************************
// Class: QvisClipWindow
//
// Purpose:
//   This class is a postable window that watches clip operator
//   attributes and always represents their current state.
//
// Notes:      
//
// Programmer: Kathleen Bonnell 
// Creation:   May 7, 2001 
//
// Modifications:
//   Brad Whitlock, Fri Apr 12 13:04:35 PST 2002
//   Made it inherit from QvisOperatorWindow.
//
//   Kathleen Bonnell, Mon Dec  6 14:35:14 PST 2004 
//   Made plane*Status be checkable QVGroupBox, instead of QButtonGroup.
//   Renamed plane*StatusCliced slots to plane*StatusToggled.
//
//   Brad Whitlock, Tue Dec 21 09:13:49 PDT 2004
//   Added Qt version-specific coding so we can still use versions older than
//   3.2.
//
//   Cyrus Harrison, Wed Mar  5 10:25:39 PST 2008
//   Removed tabWidget and slot for tabWidgetChanged 
//   (to Match Sean's changes in QvisClipWindow.C)
//
//   Cyrus Harrison, Thu Aug 21 09:48:43 PDT 2008
//   Qt4 Port.
//
// ****************************************************************************

class QvisClipWindow : public QvisOperatorWindow
{
    Q_OBJECT
public:
    QvisClipWindow(const int type,
                    ClipAttributes *subj,
                    const QString &caption = QString::null,
                    const QString &shortName = QString::null,
                    QvisNotepadArea *notepad = 0);
    virtual ~QvisClipWindow();
protected:
    virtual void CreateWindowContents();
    void UpdateWindow(bool doAll);
    virtual void GetCurrentValues(int which_widget);
private slots:
    void processPlane1Origin();
    void processPlane2Origin();
    void processPlane3Origin();
    void processPlane1Normal();
    void processPlane2Normal();
    void processPlane3Normal();
    void processCenterText();
    void processRadiusText();
    void planeInverseToggled(bool);
    void planeToolControlledClipPlaneChanged(int val);
    void qualityChanged(int);
    void sliceTypeChanged(int);
    void sphereInverseToggled(bool);
    void plane1StatusToggled(bool);
    void plane2StatusToggled(bool);
    void plane3StatusToggled(bool);

private:
    QButtonGroup *qualityGroup;
    QButtonGroup *typeGroup;
    
    QPlaneGroup  *plane1Group;
    QPlaneGroup  *plane2Group;
    QPlaneGroup  *plane3Group;
    
    QLineEdit    *centerLineEdit;
    QLineEdit    *radiusLineEdit;
    QCheckBox    *planeInverse;
    QButtonGroup *planeToolControlledClipPlane;
    QCheckBox    *sphereInverse;
    QWidget      *planeWidgets;
    QWidget      *sphereWidgets;

    ClipAttributes *atts;
};
#endif
