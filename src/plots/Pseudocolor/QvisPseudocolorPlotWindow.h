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

#ifndef QVIS_PSEUDOCOLOR_WINDOW_H
#define QVIS_PSEUDOCOLOR_WINDOW_H
#include <QvisPostableWindowObserver.h>

class QGroupBox;
class QComboBox;
class QLineEdit;
class QCheckBox;
class QButtonGroup;
class QLabel;
class QvisOpacitySlider;
class QvisColorTableWidget;
class QvisPointControl;
class QvisLineStyleWidget;
class QvisLineWidthWidget;
class QvisVariableButton;
class QSpinBox;

class Subject;
class PseudocolorAttributes;

class QvisCollapsibleLayout;

// ****************************************************************************
// Class: QvisPseudocolorPlotWindow
//
// Purpose:
//   This class is a postable window that watches pseudocolot plot
//   attributes and always represents their current state.
//
// Notes:      
//
// Programmer: Brad Whitlock
// Creation:   Fri Jul 28 17:16:38 PST 2000
//
// Modifications:
//   Kathleen Bonnell, Thu Dec 14 17:04:25 PST 2000
//   Added support for setting opacity.
//
//   Eric Brugger, Wed Mar 14 06:59:25 PST 2001
//   I added a plot type to the constructor for use with the viewer
//   proxy.
//
//   Brad Whitlock, Sat Jun 16 15:17:12 PST 2001
//   I added a color table button.
//
//   Kathleen Bonnell, Thu Oct  4 16:28:16 PDT 2001 
//   I added a limits combo box.
//
//   Brad Whitlock, Fri Feb 15 10:27:55 PDT 2002
//   Removed a method.
//
//   Jeremy Meredith, Tue Dec 10 10:22:40 PST 2002
//   Added smoothing level.
//
//   Jeremy Meredith, Fri Dec 20 11:36:03 PST 2002
//   Added scaling of point variables by a scalar field.
//
//   Kathleen Bonnell, Fri Nov 12 11:25:23 PST 2004 
//   Replace individual point-size related widgets and associated slots
//   with QvisPointControl 
//
//   Brad Whitlock, Wed Jul 20 14:23:58 PST 2005
//   Added a new slot to handle a new signal from QvisPointControl.
//
//   Jeremy Meredith, Wed Nov 26 11:28:24 EST 2008
//   Added line style/width controls.
//
//   Jeremy Meredith, Fri Feb 20 15:14:29 EST 2009
//   Added support for using per-color alpha values from a color table
//   (instead of just a single global opacity for the whole plot).
//   There's a new toggle for this, and it overrides the whole-plot opacity.
//
//   Allen Sanderson, Sun Mar  7 12:49:56 PST 2010
//   Change layout of window for 2.0 interface changes.
//
//   Kathleen Bonnell, Mon Jan 17 18:02:39 MST 2011
//   Change colorTableButton to colorTableWidget to gain invert toggle.
//
// ****************************************************************************

class QvisPseudocolorPlotWindow : public QvisPostableWindowObserver
{
    Q_OBJECT
public:
    QvisPseudocolorPlotWindow(const int type, PseudocolorAttributes *_pcAtts,
                              const QString &caption = QString::null,
                              const QString &shortName = QString::null,
                              QvisNotepadArea *notepad = 0);
    virtual ~QvisPseudocolorPlotWindow();
    virtual void CreateWindowContents();
public slots:
    virtual void apply();
    virtual void makeDefault();
    virtual void reset();
protected:

    void CreateGeometryTab(QWidget *);
    void CreateDataTab(QWidget *);
    void CreateExtrasTab(QWidget *);

    void UpdateWindow(bool doAll);
    void GetCurrentValues(int which_widget);
    void Apply(bool ignore = false);
private slots:

    void scaleClicked(int scale);
    void processSkewText();

    void limitsSelectChanged(int);
    void minToggled(bool on);
    void maxToggled(bool on);
    void processMaxLimitText();
    void processMinLimitText();

    void centeringClicked(int button);

    void colorTableClicked(bool useDefault, const QString &ctName);
    void invertColorTableToggled(bool val);

    void opacityTypeChanged(int val);
    void opacityVariableChanged(const QString &var);
    void opacityChanged(int opacity, const void*);

    void opacityMinToggled(bool);
    void opacityMaxToggled(bool);
    void processOpacityVarMin();
    void processOpacityVarMax();

    void pointTypeChanged(int index);
    void pointSizeChanged(double d);
    void pointSizePixelsChanged(int size);
    void pointSizeVarToggled(bool on);
    void pointSizeVarChanged(const QString &);

    void lineTypeChanged(int newType);
    void lineStyleChanged(int newStyle);
    void lineWidthChanged(int newWidth);

    void tubeDisplayDensityChanged(int val);
    void tubeRadiusSizeTypeChanged(int v);
    void tubeRadiusProcessText();
    void tubeRadiusVaryChanged(bool val);
    void tubeRadiusVaryVariableChanged(const QString &var);
    void tubeRadiusVaryFactorProcessText();

    void endPointTypeChanged(int newType);
    void endPointStyleChanged(int newStyle);
    void endPointRadiusSizeTypeChanged(int v);
    void endPointRadiusProcessText();
    void endPointRatioProcessText();

    void smoothingLevelChanged(int index);
    void renderSurfacesChanged(bool);
    void renderWireframeChanged(bool);
    void renderPointsChanged(bool);

    void legendToggled(bool on);
    void lightingToggled(bool on);

private:
    int                   plotType;
    PseudocolorAttributes *pcAtts;

    QvisCollapsibleLayout *propertyLayout;

    QButtonGroup          *scalingButtons;
    QLineEdit             *skewLineEdit;

    QComboBox             *limitsSelect;
    QCheckBox             *minToggle;
    QCheckBox             *maxToggle;
    QLineEdit             *maxLineEdit;
    QLineEdit             *minLineEdit;

    QButtonGroup          *centeringButtons;

    QvisColorTableWidget  *colorTableWidget;

    QComboBox *opacityType;
    QLabel    *opacityVarLabel;
    QvisVariableButton *opacityVar;
    QvisOpacitySlider *opacitySlider;
    QGroupBox *opacityMinMaxGroup;
    QCheckBox *opacityMinToggle;
    QCheckBox *opacityMaxToggle;
    QLineEdit *opacityVarMin;
    QLineEdit *opacityVarMax;

    // QButtonGroup          *opacityButtons;
    // QLabel                *opacitySliderLabel;
    // QvisOpacitySlider     *opacitySlider;

    QvisPointControl      *pointControl;

    QComboBox             *lineType;

    QLabel                *lineStyleLabel;
    QvisLineStyleWidget   *lineStyle;
    QLabel                *lineWidthLabel;
    QvisLineWidthWidget   *lineWidth;


    QLabel             *tubeDisplayDensityLabel;
    QSpinBox           *tubeDisplayDensity;
    QLabel             *tubeRadiusLabel;
    QLineEdit          *tubeRadius;
    QComboBox          *tubeRadiusSizeType;

    QCheckBox          *tubeRadiusVary;
    QLabel             *tubeRadiusVaryVariableLabel;
    QvisVariableButton *tubeRadiusVaryVariable;
    QLabel             *tubeRadiusVaryFactorLabel;
    QLineEdit          *tubeRadiusVaryFactorEdit;
  
    // QLineEdit *ribbonWidth;
    // QComboBox *ribbonSizeType;

    QComboBox *endPointType;
    QLabel    *endPointStyleLabel;
    QComboBox *endPointStyle;
    QLabel    *endPointRadiusLabel;
    QLineEdit *endPointRadius;
    QComboBox *endPointRadiusSizeType;
    QLabel    *endPointRatioLabel;
    QLineEdit *endPointRatio;

    QButtonGroup          *smoothingLevelButtons;
    QCheckBox             *renderSurfaces;
    QCheckBox             *renderWireframe;
    QCheckBox             *renderPoints;

    QCheckBox             *legendToggle;
    QCheckBox             *lightingToggle;
};
#endif
