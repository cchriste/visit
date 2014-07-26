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

#ifndef QVISMOLECULEPLOTWINDOW_H
#define QVISMOLECULEPLOTWINDOW_H

#include <QvisPostableWindowObserver.h>
#include <AttributeSubject.h>

class MoleculeAttributes;
class QLabel;
class QCheckBox;
class QLineEdit;
class QSpinBox;
class QVBox;
class QButtonGroup;
class QvisColorTableButton;
class QvisOpacitySlider;
class QvisColorButton;
class QvisLineStyleWidget;
class QvisLineWidthWidget;
class QvisVariableButton;
class QComboBox;
class QGroupBox;

// ****************************************************************************
// Class: QvisMoleculePlotWindow
//
// Purpose: 
//   Defines QvisMoleculePlotWindow class.
//
// Notes:      This class was automatically generated!

// Programmer: Jeremy Meredith
// Creation:   March 23, 2006
//
// Modifications:
//   Cyrus Harrison, Fri Jul 18 14:38:14 PDT 2008
//   Initial Qt4 Port.
//
//   Allen Sanderson, Sun Mar  7 12:49:56 PST 2010
//   Change layout of window for 2.0 interface changes.
//
// ****************************************************************************

class QvisMoleculePlotWindow : public QvisPostableWindowObserver
{
    Q_OBJECT
  public:
    QvisMoleculePlotWindow(const int type,
                         MoleculeAttributes *subj,
                         const QString &caption = QString::null,
                         const QString &shortName = QString::null,
                         QvisNotepadArea *notepad = 0);
    virtual ~QvisMoleculePlotWindow();
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
    void drawAtomsAsChanged(int val);
    void scaleRadiusByChanged(int val);
    void drawBondsAsChanged(int val);
    void colorBondsChanged(int val);
    void bondSingleColorChanged(const QColor &color);
    void radiusVariableChanged(const QString &varName);
    void radiusScaleFactorProcessText();
    void radiusFixedProcessText();
    void atomSphereQualityChanged(int val);
    void bondCylinderQualityChanged(int val);
    void bondRadiusProcessText();
    void bondLineWidthChanged(int style);
    void bondLineStyleChanged(int style);
    void elementColorTableChanged(bool useDefault, const QString &ctName);
    void residueTypeColorTableChanged(bool useDefault, const QString &ctName);
    void residueSequenceColorTableChanged(bool useDefault, const QString &ctName);
    void continuousColorTableChanged(bool useDefault, const QString &ctName);
    void legendToggled(bool val);
    void minFlagChanged(bool val);
    void scalarMinProcessText();
    void maxFlagChanged(bool val);
    void scalarMaxProcessText();
  private:
    int                   plotType;
    QComboBox            *drawAtomsAs;
    QComboBox            *scaleRadiusBy;
    QComboBox            *drawBondsAs;
    QButtonGroup         *colorBondsGroup;
    QWidget              *colorBondsWidget;
    QvisColorButton      *bondSingleColor;
    QvisVariableButton   *radiusVariable;
    QLineEdit            *radiusScaleFactor;
    QLineEdit            *radiusFixed;
    QComboBox            *atomSphereQuality;
    QComboBox            *bondCylinderQuality;
    QLineEdit            *bondRadius;
    QvisLineWidthWidget  *bondLineWidth;
    QvisLineStyleWidget  *bondLineStyle;
    QvisColorTableButton *elementColorTable;
    QvisColorTableButton *residueTypeColorTable;
    QvisColorTableButton *residueSequenceColorTable;
    QvisColorTableButton *continuousColorTable;
    QCheckBox            *legendToggle;
    QCheckBox            *minFlag;
    QLineEdit            *scalarMin;
    QCheckBox            *maxFlag;
    QLineEdit            *scalarMax;
    QLabel               *drawAtomsAsLabel;
    QLabel               *scaleRadiusByLabel;
    QLabel               *drawBondsAsLabel;
    QLabel               *colorBondsLabel;
    QLabel               *radiusVariableLabel;
    QLabel               *radiusScaleFactorLabel;
    QLabel               *radiusFixedLabel;
    QLabel               *atomSphereQualityLabel;
    QLabel               *bondCylinderQualityLabel;
    QLabel               *bondRadiusLabel;
    QLabel               *bondLineWidthLabel;
    QLabel               *bondLineStyleLabel;
    QLabel               *elementColorTableLabel;
    QLabel               *residueTypeColorTableLabel;
    QLabel               *residueSequenceColorTableLabel;
    QLabel               *continuousColorTableLabel;
    QLabel               *minFlagLabel;
    QLabel               *maxFlagLabel;

    MoleculeAttributes   *atts;
};



#endif
