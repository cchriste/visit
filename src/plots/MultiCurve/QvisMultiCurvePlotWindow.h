/*****************************************************************************
*
* Copyright (c) 2000 - 2010, Lawrence Livermore National Security, LLC
* Produced at the Lawrence Livermore National Laboratory
* LLNL-CODE-400124
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

#ifndef QVISMULTICURVEPLOTWINDOW_H
#define QVISMULTICURVEPLOTWINDOW_H

#include <QvisPostableWindowObserver.h>
#include <AttributeSubject.h>

class MultiCurveAttributes;
class QButtonGroup;
class QCheckBox;
class QGroupBox;
class QLabel;
class QLineEdit;
class QSpinBox;
class QVBox;
class QvisColorButton;
class QvisColorManagerWidget;
class QvisColorTableButton;
class QvisLineStyleWidget;
class QvisLineWidthWidget;
class QvisOpacitySlider;
class QvisVariableButton;

// ****************************************************************************
// Class: QvisMultiCurvePlotWindow
//
// Purpose:
//    Defines QvisMultiCurvePlotWindow class.
//
// Notes:      Autogenerated by xml2window.
//
// Programmer: xml2window
// Creation:   omitted
//
// Modifications:
//   Eric Brugger, Wed Jan 21 08:20:18 PST 2009
//   I added yAxisTitleFormat, useYAxisRange, and yAxisRange.  I changed
//   markerVariable from a variable button to a text field.
//   
//   Eric Brugger, Wed Feb 18 07:56:29 PST 2009
//   I added displayIds and idVariable.
//
//   Eric Brugger, Fri Feb 20 16:13:15 PST 2009
//   I added displayLegend.
//
//   Eric Brugger, Thu Mar  5 17:46:26 PST 2009
//   I replaced useYAxisRange and yAxisRange with useYAxisTickSpacing and
//   yAxisTickSpacing.
//
// ****************************************************************************

class QvisMultiCurvePlotWindow : public QvisPostableWindowObserver
{
    Q_OBJECT
  public:
    QvisMultiCurvePlotWindow(const int type,
                         MultiCurveAttributes *subj,
                         const QString &caption = QString::null,
                         const QString &shortName = QString::null,
                         QvisNotepadArea *notepad = 0);
    virtual ~QvisMultiCurvePlotWindow();
    virtual void CreateWindowContents();
  public slots:
    virtual void apply();
    virtual void makeDefault();
    virtual void reset();
  protected:
    void UpdateWindow(bool doAll);
    void GetCurrentValues(int which_widget);
    void Apply(bool ignore = false);

    bool UpdateMultipleAreaColors();
  private slots:
    void colorModeChanged(int index);
    void singleColorChanged(const QColor &color);
    void singleColorOpacityChanged(int opacity);
    void multipleColorChanged(const QColor &color, int index);
    void opacityChanged(int opacity, int index);
    void lineStyleChanged(int style);
    void lineWidthChanged(int style);
    void yAxisTitleFormatProcessText();
    void useYAxisTickSpacingChanged(bool val);
    void yAxisTickSpacingProcessText();
    void displayMarkersChanged(bool val);
    void markerVariableProcessText();
    void displayIdsChanged(bool val);
    void idVariableProcessText();
    void displayLegendChanged(bool val);
  private:
    int                     plotType;
    QGroupBox              *curveColorGroup;
    QButtonGroup           *colorModeButtons;
    QvisColorButton        *singleColor;
    QvisOpacitySlider      *singleColorOpacity;
    QvisColorManagerWidget *multipleColors;
    QvisLineStyleWidget    *lineStyle;
    QvisLineWidthWidget    *lineWidth;
    QLineEdit              *yAxisTitleFormat;
    QCheckBox              *useYAxisTickSpacing;
    QLineEdit              *yAxisTickSpacing;
    QCheckBox              *displayMarkers;
    QLineEdit              *markerVariable;
    QCheckBox              *displayIds;
    QLineEdit              *idVariable;
    QCheckBox              *displayLegend;
    QLabel                 *lineStyleLabel;
    QLabel                 *lineWidthLabel;
    QLabel                 *yAxisTitleFormatLabel;
    QLabel                 *markerVariableLabel;
    QLabel                 *idVariableLabel;

    MultiCurveAttributes   *atts;
};



#endif
