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

#ifndef QVISMESHMANAGEMENTWINDOW_H
#define QVISMESHMANAGEMENTWINDOW_H
#include <gui_exports.h>
#include <QvisPostableWindowObserver.h>
#include <AttributeSubject.h>

// Forward declarations.
class MeshManagementAttributes;
class QButtonGroup;
class QRadioButton;
class QCheckBox;
class QGroupBox;
class QLabel;
class QLineEdit;
class QTabWidget;
class QVBox;

// ****************************************************************************
// Class: QvisMeshManagementWindow
//
// Purpose: Creates window for mesh management controls 
//
// Programmer: Mark C. Miller 
// Creation:   November 5, 2005
//
// Modifications:
//
//    Mark C. Miller, Sun Dec  3 12:20:11 PST 2006
//    Added makeDefault and reset slots
//
//    Mark C. Miller, Wed Dec 19 11:32:58 PST 2007
//    Made Qt objects and visual controls input to mmatts a little more
//    user-friendly. However, mmatts themselves were not changed.
//
//    Brad Whitlock, Wed Apr  9 11:34:22 PDT 2008
//    QString for caption, shortName.
//
//   Cyrus Harrison, Wed Jul  2 11:16:25 PDT 2008
//   Initial Qt4 Port.
//
//   Cyrus Harrison, Thu Dec 18 09:29:37 PST 2008
//   Removed tabSelected slot b/c it was an empty method.
//
//   Jeremy Meredith, Fri Feb 26 14:13:08 EST 2010
//   Added a new "multi-pass" discretization algorithm
//
// ****************************************************************************

class GUI_API QvisMeshManagementWindow : public QvisPostableWindowObserver
{
    Q_OBJECT
public:
    QvisMeshManagementWindow(MeshManagementAttributes *subj,
                        const QString &caption = QString::null,
                        const QString &shortName = QString::null,
                        QvisNotepadArea *notepad = 0);
    virtual ~QvisMeshManagementWindow();
    virtual void CreateWindowContents();
public slots:
    virtual void apply();
    virtual void makeDefault();
    virtual void reset();
protected:
    void UpdateWindow(bool doAll);
    void Apply(bool ignore = false);
    void GetCurrentValues(const QWidget *widget = 0);
private slots:
    void processSmallestZoneText();
    void processSmallestZoneText(const QString &);
    void processFlatEnoughText();
    void processFlatEnoughText(const QString &);
    void renderCSGDirectChanged(bool);
    void discretizeBoundaryOnlyChanged(bool);
    void discretizationModeChanged(int);
private:
    MeshManagementAttributes *mmAtts;

    QWidget          *pageCSGGroup;
    QCheckBox        *renderCSGDirect;
    QCheckBox        *discretizeBoundaryOnly;
    QLabel           *discretizeModeLabel;
    QButtonGroup     *discretizationMode;
    QRadioButton     *discretizeMultiPass;
    QRadioButton     *discretizeUniform;
    QRadioButton     *discretizeAdaptive;
    QLabel           *smallestZoneLabel;
    QLineEdit        *smallestZoneLineEdit;
    QLabel           *flatEnoughLabel;
    QLineEdit        *flatEnoughLineEdit;

    QTabWidget       *tabs;

};
#endif
