/*****************************************************************************
*
* Copyright (c) 2000 - 2011, Lawrence Livermore National Security, LLC
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

#ifndef QVISPOINCAREPLOTWINDOW_H
#define QVISPOINCAREPLOTWINDOW_H

#include <QvisPostableWindowObserver.h>
#include <AttributeSubject.h>

class PoincareAttributes;
class QTabWidget;
class QLabel;
class QCheckBox;
class QLineEdit;
class QSpinBox;
class QVBox;
class QButtonGroup;
class QRadioButton;
class QComboBox;
class QGroupBox;
class QvisColorTableButton;
class QvisOpacitySlider;
class QvisColorButton;
class QvisLineWidthWidget;
class QvisVariableButton;
class QvisPointControl;
class QvisLineStyleWidget;

// ****************************************************************************
// Class: QvisPoincarePlotWindow
//
// Purpose:
//    Defines QvisPoincarePlotWindow class.
//
// Notes:      Autogenerated by xml2window.
//
// Programmer: xml2window
// Creation:   omitted
//
// Modifications:
//   
//   Allen Sanderson, Sun Mar  7 12:49:56 PST 2010
//   Change layout of window for 2.0 interface changes.
//
//   Dave Pugmire, Thu Jul  8 09:03:20 EDT 2010
//   Add force node centering option.
//
// ****************************************************************************

class QvisPoincarePlotWindow : public QvisPostableWindowObserver
{
    Q_OBJECT
  public:
    QvisPoincarePlotWindow(const int type,
                         PoincareAttributes *subj,
                         const QString &caption = QString::null,
                         const QString &shortName = QString::null,
                         QvisNotepadArea *notepad = 0);
    virtual ~QvisPoincarePlotWindow();
    virtual void CreateWindowContents();
  public slots:
    virtual void apply();
    virtual void makeDefault();
    virtual void reset();
  protected:
    void UpdateWindow(bool doAll);
    void GetCurrentValues(int which_widget);
    void UpdateIntegrationAttributes();
    void UpdateAlgorithmAttributes();
    void UpdateMeshTypeAttributes();
    void Apply(bool ignore = false);
  private slots:
    void sourceTypeChanged(int val);
    void pointSourceProcessText();
    void lineStartProcessText();
    void lineEndProcessText();
    void pointDensityChanged(int val);
    void minPuncturesChanged(int val);
    void maxPuncturesChanged(int val);
    void puncturePlaneChanged(int val);
    void integrationTypeChanged(int val);
    void maxStepLengthProcessText();
    void relTolProcessText();
    void absTolProcessText();
    void coordinateButtonGroupChanged(int val);
    void analysisChanged(int val);
    void maximumToroidalWindingChanged(int val);
    void overrideToroidalWindingChanged(int val);
    void overridePoloidalWindingChanged(int val);
    void windingPairConfidenceProcessText();
    void rationalTemplateSeedParmProcessText();
    void adjustPlaneChanged(int val);
    void overlapsChanged(int val);
    void meshTypeChanged(int val);
    void numberPlanesChanged(int val);
    void singlePlaneProcessText();
    void minProcessText();
    void maxProcessText();
    void minFlagChanged(bool val);
    void maxFlagChanged(bool val);
    void colorTypeChanged(int val);
    void singleColorChanged(const QColor &color);
    void colorTableNameChanged(bool useDefault, const QString &ctName);
    void setOpaacityClicked(int on);
    void changedOpacity(int opacity, const void *);
    void dataValueChanged(int val);
    void showOPointsChanged(bool val);
    void OPointMaxIterationsChanged(int val);
    void showIslandsChanged(bool val);
    void showChaoticChanged(bool val);
    void showLinesChanged(bool val);
    void lineWidthChanged(int val);
    void lineStyleChanged(int val);
    void showPointsChanged(bool val);
    void pointSizeChanged(double val);
    void pointSizePixelsChanged(int val);
    void pointTypeChanged(int val);
    void verboseFlagChanged(bool val);
    void show1DPlotsChanged(bool val);
    void legendToggled(bool val);
    void lightingToggled(bool val);
    void streamlineAlgorithmChanged(int val);
    void maxSLCountChanged(int val);
    void maxDomainCacheChanged(int val);
    void workGroupSizeChanged(int val);
    void forceNodalChanged(bool);
  private:
    QTabWidget      *propertyTabs;

    int plotType;
    QSpinBox  *minPunctures;
    QSpinBox  *maxPunctures;
    QLabel       *puncturePlaneLabel;
    QWidget      *puncturePlane;
    QButtonGroup *puncturePlaneButtonGroup;
    QWidget   *sourceType;
    QComboBox *sourceTypeCombo;
    QLineEdit *pointSource;
    QLineEdit *lineStart;
    QLineEdit *lineEnd;
    QSpinBox  *pointDensity;
    QWidget   *integrationType;
    QComboBox *integrationTypeCombo;
    QLineEdit *maxStepLength;
    QLineEdit *relTol;
    QLineEdit *absTol;
    QButtonGroup *coordinateButtonGroup;
    QWidget      *analysis;
    QButtonGroup *analysisButtonGroup;
    QSpinBox *maximumToroidalWinding;
    QSpinBox *overrideToroidalWinding;
    QSpinBox *overridePoloidalWinding;
    QLineEdit *windingPairConfidence;
    QLineEdit *rationalTemplateSeedParm;
    QWidget      *overlaps;
    QButtonGroup *overlapsButtonGroup;
    QWidget      *meshType;
    QComboBox *meshTypeCombo;
    QSpinBox  *numberPlanes;
    QLineEdit *singlePlane;
    QSpinBox *adjustPlane;
    QLineEdit *min;
    QLineEdit *max;
    QCheckBox *minFlag;
    QCheckBox *maxFlag;
    QWidget      *colorType;
    QButtonGroup *colorTypeButtonGroup;
    QvisColorButton *singleColor;
    QvisColorTableButton *colorTableName;
    QLabel               *opacityButtonsLabel;
    QButtonGroup         *opacityButtons;
    QRadioButton         *opacityButtonSetExplicit;
    QRadioButton         *opacityButtonColorTable;
    QLabel               *opacitySliderLabel;
    QvisOpacitySlider    *opacitySlider;
    QWidget      *dataValue;
    QComboBox *dataValueCombo;
    QCheckBox *showOPoints;
    QLabel *OPointMaxIterationsLabel;
    QSpinBox *OPointMaxIterations;
    QCheckBox *showChaotic;
    QCheckBox *showIslands;
    QLabel *lineWidthLabel, *lineStyleLabel;
    QCheckBox *showPoints;
    QvisPointControl *pointControl;
    QCheckBox *showLines;
    QvisLineWidthWidget *lineWidth;
    QvisLineStyleWidget *lineStyle;

    QCheckBox *show1DPlots;
    QCheckBox *verboseFlag;
    QCheckBox *legendToggle;
    QCheckBox *lightingToggle;
    QLabel *minPuncturesLabel;
    QLabel *maxPuncturesLabel;
    QLabel *sourceTypeLabel;
    QLabel *pointSourceLabel;
    QLabel *lineStartLabel;
    QLabel *lineEndLabel;
    QLabel *pointDensityLabel;
    QLabel *integrationTypeLabel;
    QLabel *maxStepLengthLabel;
    QLabel *relTolLabel;
    QLabel *absTolLabel;
    QLabel    *forceNodalLabel;
    QCheckBox *forceNodal;
    QLabel *analysisLabel;
    QLabel *maximumToroidalWindingLabel;
    QLabel *overrideToroidalWindingLabel;
    QLabel *overridePoloidalWindingLabel;
    QLabel *windingPairConfidenceLabel;
    QLabel *rationalTemplateSeedParmLabel;
    QLabel *adjustPlaneLabel;
    QLabel *overlapsLabel;
    QLabel *meshTypeLabel;
    QLabel *numberPlanesLabel;
    QLabel *singlePlaneLabel;
    QLabel *minLabel;
    QLabel *maxLabel;
    QLabel *colorTypeLabel;
    QLabel *singleColorLabel;
    QLabel *colorTableNameLabel;
    QLabel *dataValueLabel;

    QLabel    *slAlgoLabel;
    QComboBox *slAlgo;
    QLabel    *maxSLCountLabel;
    QSpinBox  *maxSLCount;
    QLabel    *maxDomainCacheLabel;
    QSpinBox  *maxDomainCache;
    QLabel    *workGroupSizeLabel;
    QSpinBox  *workGroupSize;

    QWidget *firstTab;
    QWidget *secondTab;
    QWidget *thirdTab;
    QWidget *fourthTab;
    QGroupBox *sourceGroup;

    PoincareAttributes *atts;
};
#endif
