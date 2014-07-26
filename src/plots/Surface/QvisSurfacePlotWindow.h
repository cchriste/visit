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

#ifndef QVIS_SURFACE_PLOT_WINDOW_H
#define QVIS_SURFACE_PLOT_WINDOW_H

#include <QvisPostableWindowObserver.h>
#include <AttributeSubject.h>

// Forward declarations.
class SurfaceAttributes;
class QButtonGroup;
class QCheckBox;
class QComboBox;
class QGroupBox;
class QLabel;
class QLineEdit;
class QRadioButton;
class QvisColorButton;
class QvisColorManagerWidget;
class QvisColorTableWidget;
class QvisLineStyleWidget;
class QvisLineWidthWidget;
class QvisOpacitySlider;

// ****************************************************************************
// Class: QvisSurfacePlotWindow
//
// Purpose: 
//   This class is an observer window that watches material plot attributes
//   and always represents their current state.
//
// Notes:
//
// Programmer: Kathleen Bonnell 
// Creation:   March 06, 2001 
//
// Modifications:
//   Eric Brugger, Fri Mar 16 14:42:45 PST 2001
//   I added a plot type to the constructor for use with the viewer
//   proxy.
//
//   Brad Whitlock, Sat Jun 16 18:21:34 PST 2001
//   I added color table stuff.
//
//   Kathleen Bonnell, Thu Oct 11 12:45:30 PDT 2001
//   Changed limits (min/max) related members to reflect distinction between
//   limits used for scaling and limits used for coloring.   
//   Added limitsSelect
//
//   Kathleen Bonnell, Thu Mar 28 14:03:19 PST 2002 
//   Revert back to using one min/max for both scaling and coloring purposes. 
//   
//   Allen Sanderson, Sun Mar  7 12:49:56 PST 2010
//   Change layout of window for 2.0 interface changes.
//
//   Kathleen Bonnell, Mon Jan 17 17:54:48 MST 2011
//   Change colorTableButton to colorTableWidget to gain invert toggle.
//
// ****************************************************************************

class QvisSurfacePlotWindow : public QvisPostableWindowObserver
{
    Q_OBJECT
public:
    QvisSurfacePlotWindow(const int type, SurfaceAttributes *surfaceAtts_,
                          const QString &caption = QString::null,
                          const QString &shortName = QString::null,
                          QvisNotepadArea *notepad = 0);
    virtual ~QvisSurfacePlotWindow();
    virtual void CreateWindowContents();
public slots:
    virtual void apply();
    virtual void makeDefault();
    virtual void reset();
protected:
    void UpdateWindow(bool doAll);
    void GetCurrentValues(int which_widget);
    void Apply(bool ignore = false);

private slots:
    void lineStyleChanged(int newStyle);
    void lineWidthChanged(int newWidth);
    void legendToggled(bool val);
    void lightingToggled(bool val);
    void scaleClicked(int button);
    void surfaceToggled(bool val);
    void wireframeToggled(bool val);
    void surfaceColorChanged(const QColor &color);
    void wireframeColorChanged(const QColor &color);
    void colorModeChanged(int ); 
    void minToggled(bool);
    void processMinLimitText();
    void maxToggled(bool);
    void processMaxLimitText();
    void processSkewText();
    void colorTableClicked(bool useDefault, const QString &ctName);
    void invertColorTableToggled(bool val);
    void limitsSelectChanged(int);
private:
    int                     plotType;
    SurfaceAttributes      *surfaceAtts;

    // Surface controls
    QGroupBox              *surfaceGroup;
    QButtonGroup           *colorModeButtons;
    QvisColorButton        *surfaceColor;
    QvisColorTableWidget   *colorTableWidget;

    // Wireframe controls
    QGroupBox              *wireframeGroup;
    QvisLineStyleWidget    *lineStyle;
    QvisLineWidthWidget    *lineWidth;
    QvisColorButton        *wireframeColor;

    // Scale controls
    QButtonGroup           *scalingButtons;
    QLineEdit              *skewLineEdit;

    // Limits controls
    QComboBox              *limitsSelect;
    QCheckBox              *minToggle;
    QLineEdit              *minLineEdit;
    QCheckBox              *maxToggle;
    QLineEdit              *maxLineEdit;

    QCheckBox              *legendToggle;
    QCheckBox              *lightingToggle;
};

#endif
