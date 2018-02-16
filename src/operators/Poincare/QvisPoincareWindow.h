/*****************************************************************************
*
* Copyright (c) 2000 - 2018, Lawrence Livermore National Security, LLC
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

#ifndef QVIS_POINCARE_WINDOW_H
#define QVIS_POINCARE_WINDOW_H

#include <QvisOperatorWindow.h>
#include <AttributeSubject.h>

// Forward declarations
class PoincareAttributes;

class QLabel;
class QCheckBox;
class QLineEdit;
class QSpinBox;
class QButtonGroup;
class QComboBox;
class QGroupBox;
class QPushButton;
class QListWidget;
class QListWidgetItem;

// ****************************************************************************
// Class: QvisPoincareWindow
//
// Purpose:
//    Defines QvisPoincareWindow class.
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

class QvisPoincareWindow : public QvisOperatorWindow
{
    Q_OBJECT

  public:
    QvisPoincareWindow(const int type,
                         PoincareAttributes *subj,
                         const QString &caption = QString::null,
                         const QString &shortName = QString::null,
                         QvisNotepadArea *notepad = 0);
    virtual ~QvisPoincareWindow();
    virtual void CreateWindowContents();

  protected:
    void CreateIntegrationTab(QWidget *);
    void CreateAnalysisTab(QWidget *);
    void CreateAppearanceTab(QWidget *);
    void CreateAdvancedTab(QWidget *);

    void UpdateWindow(bool doAll);
    void GetCurrentValues(int which_widget);
    void UpdateFieldAttributes();
    void UpdateIntegrationAttributes();
    void UpdateAlgorithmAttributes();
    void UpdateMeshTypeAttributes();

  private slots:
    // Fieldlines
    void sourceTypeChanged(int val);
    void pointSourceProcessText();
    void pointListProcessText();
    void lineStartProcessText();
    void lineEndProcessText();
    void pointDensityChanged(int val);

    void pointListClicked(QListWidgetItem*);
    void pointListDoubleClicked(QListWidgetItem*);
    void addPoint();
    void deletePoint();
    void deletePoints();
    void readPoints();
    void textChanged(const QString &currentText);

  void fieldTypeChanged(int val);
    void fieldConstantProccessText();
    void velocitySourceProcessText();

    void integrationTypeChanged(int val);
    void maxStepLengthProcessText();
    void limitMaxTimeStepChanged(bool val);
    void maxTimeStepProcessText();
    void relTolProcessText();
    void absTolProcessText();
    void absTolSizeTypeChanged(int);

    void minPuncturesChanged(int val);
    void maxPuncturesChanged(int val);

    void puncturePlotTypeChanged(int val);
    void maxStepsProcessText();
    void maxTimeProcessText();
    void limitMaxTimeChanged(bool val);
    void puncturePeriodToleranceProcessText();

    void puncturePlaneChanged(int val);
//    void coordinateButtonGroupChanged(int val);

    // Analysis
    void analysisChanged(int val);

    void maximumToroidalWindingChanged(int val);
    void overrideToroidalWindingChanged(int val);
    void overridePoloidalWindingChanged(int val);
    void windingPairConfidenceProcessText();
    void rationalSurfaceFactorProcessText();

    void showRationalSurfacesChanged(bool val);
    void rationalSurfaceMaxIterationsChanged(int val);

    void showOPointsChanged(bool val);
    void OPointMaxIterationsChanged(int val);

    void performOLineAnalysisChanged(bool val);
    void OLineToroidalWindingChanged(int val);
    void OLineAxisFileNameProcessText();
    void OLineAxisFileDialogButtonClicked();

    void showChaoticChanged(bool val);
    void showIslandsChanged(bool val);
    void summaryFlagChanged(bool val);
    void verboseFlagChanged(bool val);
    void show1DPlotsChanged(bool val);

    // Appearance
    void dataValueChanged(int val);
    void meshTypeChanged(int val);
    void showLinesChanged(bool val);
    void showPointsChanged(bool val);
    void numberPlanesChanged(int val);
    void singlePlaneProcessText();
    void overlapsChanged(int val);

    // Advanced
    void parallelAlgorithmChanged(int val);
    void maxSLCountChanged(int val);
    void maxDomainCacheChanged(int val);
    void workGroupSizeChanged(int val);

    void icButtonGroupChanged(int val);
    void pathlineOverrideStartingTimeFlagChanged(bool val);
    void pathlineOverrideStartingTimeProcessText();
    void pathlinePeriodProcessText();
    void pathlineCMFEButtonGroupChanged(int val);

    void issueWarningForMaxStepsChanged(bool);
    void issueWarningForStepsizeChanged(bool);
    void issueWarningForStiffnessChanged(bool);
    void issueWarningForCriticalPointsChanged(bool);
    void criticalPointThresholdProcessText();

//    void forceNodalChanged(bool);

  private:
    int plotType;

    // Fieldlines
    QWidget   *sourceType;
    QLabel    *sourceTypeLabel;
    QComboBox *sourceTypeCombo;
    QLineEdit *velocitySource;
    QLabel    *velocitySourceLabel;
    QLineEdit *pointSource;
    QLabel    *pointSourceLabel;
    QListWidget *pointList;
    QPushButton *pointListDelPoint, *pointListDelAllPoints, *pointListAddPoint, *pointListReadPoints;
    QLineEdit *lineStart;
    QLabel    *lineStartLabel;
    QLineEdit *lineEnd;
    QLabel    *lineEndLabel;
    QSpinBox  *pointDensity;
    QLabel    *pointDensityLabel;

    QComboBox *fieldType;
    QLabel    *fieldConstantLabel;
    QLineEdit *fieldConstant;
//    QCheckBox *forceNodal;

    QComboBox *integrationType;
    QLabel *integrationTypeLabel;
    QCheckBox *limitMaxTimeStep;
    QLineEdit *maxStepLength;
    QLabel    *maxStepLengthLabel;
    QLineEdit *maxTimeStep;
    QLineEdit *relTol;
    QLabel    *relTolLabel;
    QLineEdit *absTol;
    QComboBox *absTolSizeType;
    QLabel    *absTolLabel;
    QLabel    *minPuncturesLabel;
    QSpinBox  *minPunctures;
    QLabel    *maxPuncturesLabel;
    QSpinBox  *maxPunctures;
    QLabel    *puncturePlaneLabel;

    QLabel       *puncturePlotTypeLabel;
    QWidget      *puncturePlotType;
    QButtonGroup *puncturePlotTypeButtonGroup;
    QLineEdit *puncturePeriodTolerance;
    QLabel    *puncturePeriodToleranceLabel;
    QLineEdit *maxSteps;
    QLabel    *maxStepsLabel;
    QCheckBox *limitMaxTime;
    QLineEdit *maxTime;

    QWidget      *puncturePlane;
    QButtonGroup *puncturePlaneButtonGroup;

    // QButtonGroup *coordinateButtonGroup;

    // Analysis
    QLabel       *analysisLabel;
    QWidget      *analysis;
    QButtonGroup *analysisButtonGroup;
    QLabel       *maximumToroidalWindingLabel;
    QSpinBox     *maximumToroidalWinding;
    QLabel       *overrideToroidalWindingLabel;
    QSpinBox     *overrideToroidalWinding;
    QLabel       *overridePoloidalWindingLabel;
    QSpinBox     *overridePoloidalWinding;
    QLabel       *windingPairConfidenceLabel;
    QLineEdit    *windingPairConfidence;
    QLabel       *rationalSurfaceFactorLabel;
    QLineEdit    *rationalSurfaceFactor;

    QCheckBox *showRationalSurfaces;
    QLabel    *rationalSurfaceMaxIterationsLabel;
    QSpinBox  *rationalSurfaceMaxIterations;

    QCheckBox *showOPoints;
    QLabel    *OPointMaxIterationsLabel;
    QSpinBox  *OPointMaxIterations;

    QCheckBox   *performOLineAnalysis;
    QLabel      *OLineToroidalWindingLabel;
    QSpinBox    *OLineToroidalWinding;
    QPushButton *OLineAxisFileDialogButton;
    QLineEdit   *OLineAxisFileName;

    QCheckBox *showChaotic;
    QCheckBox *showIslands;
    QCheckBox *show1DPlots;
    QCheckBox *summaryFlag;
    QCheckBox *verboseFlag;

    // Appearance
    QLabel    *dataValueLabel;
    QWidget   *dataValue;
    QComboBox *dataValueCombo;

    QLabel       *meshTypeLabel;
    QWidget      *meshType;
    QComboBox    *meshTypeCombo;
    QCheckBox    *showLines;
    QCheckBox    *showPoints;
    QLabel       *numberPlanesLabel;
    QSpinBox     *numberPlanes;
    QLabel       *singlePlaneLabel;
    QLineEdit    *singlePlane;
    QWidget      *overlaps;
    QLabel       *overlapsLabel;
    QButtonGroup *overlapsButtonGroup;

    // Advanced
    QLabel    *parallelAlgoLabel;
    QComboBox *parallelAlgo;
    QLabel    *maxSLCountLabel;
    QSpinBox  *maxSLCount;
    QLabel    *maxDomainCacheLabel;
    QSpinBox  *maxDomainCache;
    QLabel    *workGroupSizeLabel;
    QSpinBox  *workGroupSize;

    QButtonGroup *icButtonGroup;
    QCheckBox *pathlineOverrideStartingTimeFlag;
    QLineEdit *pathlineOverrideStartingTime;
    QLineEdit *pathlinePeriod;
    QButtonGroup *pathlineCMFEButtonGroup;

    QCheckBox *issueWarningForMaxSteps;
    QCheckBox *issueWarningForStepsize;
    QCheckBox *issueWarningForStiffness;
    QCheckBox *issueWarningForCriticalPoints;
    QLineEdit *criticalPointThreshold;
    QLabel    *criticalPointThresholdLabel;

    QGroupBox *sourceGroup;

    PoincareAttributes *atts;
};
#endif
