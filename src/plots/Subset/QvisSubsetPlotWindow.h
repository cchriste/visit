/*****************************************************************************
*
* Copyright (c) 2000 - 2006, The Regents of the University of California
* Produced at the Lawrence Livermore National Laboratory
* All rights reserved.
*
* This file is part of VisIt. For details, see http://www.llnl.gov/visit/. The
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
*    documentation and/or materials provided with the distribution.
*  - Neither the name of the UC/LLNL nor  the names of its contributors may be
*    used to  endorse or  promote products derived from  this software without
*    specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT  HOLDERS AND CONTRIBUTORS "AS IS"
* AND ANY EXPRESS OR  IMPLIED WARRANTIES, INCLUDING,  BUT NOT  LIMITED TO, THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND  FITNESS FOR A PARTICULAR  PURPOSE
* ARE  DISCLAIMED.  IN  NO  EVENT  SHALL  THE  REGENTS  OF  THE  UNIVERSITY OF
* CALIFORNIA, THE U.S.  DEPARTMENT  OF  ENERGY OR CONTRIBUTORS BE  LIABLE  FOR
* ANY  DIRECT,  INDIRECT,  INCIDENTAL,  SPECIAL,  EXEMPLARY,  OR CONSEQUENTIAL
* DAMAGES (INCLUDING, BUT NOT  LIMITED TO, PROCUREMENT OF  SUBSTITUTE GOODS OR
* SERVICES; LOSS OF  USE, DATA, OR PROFITS; OR  BUSINESS INTERRUPTION) HOWEVER
* CAUSED  AND  ON  ANY  THEORY  OF  LIABILITY,  WHETHER  IN  CONTRACT,  STRICT
* LIABILITY, OR TORT  (INCLUDING NEGLIGENCE OR OTHERWISE)  ARISING IN ANY  WAY
* OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
* DAMAGE.
*
*****************************************************************************/

#ifndef QVIS_SUBSET_PLOT_WINDOW_H
#define QVIS_SUBSET_PLOT_WINDOW_H
#include <QvisPostableWindowObserver.h>

class SubsetAttributes;
class QButtonGroup;
class QCheckBox;
class QGroupBox;
class QLabel;
class QListBox;
class QSlider;
class QvisColorButton;
class QvisColorTableButton;
class QvisLineStyleWidget;
class QvisLineWidthWidget;
class QvisOpacitySlider;
class QvisPointControl;

// ****************************************************************************
// Class: QvisSubsetPlotWindow
//
// Purpose: 
//   This class is an observer window that watches subset plot attributes
//   and always represents their current state.
//
// Notes:  Copied from QvisMaterialPlotWindow with minor changes.
//
// Programmer: Kathleen Bonnell 
// Creation:   October 16, 2001 
//
// Modifications:
//   Brad Whitlock, Wed Jan 30 11:04:17 PDT 2002
//   Rewrote the way multiple colors are assigned.
//
//   Jeremy Meredith, Tue Mar 12 17:15:49 PST 2002
//   Added a toggle for wireframe mode.
//
//   Jeremy Meredith, Thu Mar 14 17:49:58 PST 2002
//   Added a toggle for internal surfaces.
//
//   Brad Whitlock, Mon Nov 25 17:04:48 PST 2002
//   I added a color table button.
//
//   Jeremy Meredith, Tue Dec 10 10:22:40 PST 2002
//   Added smoothing level.
//
//   Kathleen Bonnell, Fri Nov 12 11:35:11 PST 2004 
//   Added pointControl and its associated slots.  Also added GetCurrentValues. 
//
//   Brad Whitlock, Wed Jul 20 14:23:58 PST 2005
//   Added a new slot to handle a new signal from QvisPointControl.
//
// ****************************************************************************

class QvisSubsetPlotWindow : public QvisPostableWindowObserver
{
    Q_OBJECT
public:
    QvisSubsetPlotWindow(const int type, SubsetAttributes *subsetAtts_,
                           const char *caption = 0,
                           const char *shortName = 0,
                           QvisNotepadArea *notepad = 0);
    virtual ~QvisSubsetPlotWindow();
    virtual void CreateWindowContents();
public slots:
    virtual void apply();
    virtual void makeDefault();
    virtual void reset();
protected:
    void UpdateWindow(bool doAll);
    void Apply(bool ignore = false);
    void SetMultipleColorWidgets(int index);
    void UpdateMultipleArea();
    void UpdateItem(int i);
    void GetCurrentValues(int which_widget);
private slots:
    void lineStyleChanged(int newStyle);
    void lineWidthChanged(int newWidth);
    void legendToggled(bool val);
    void wireframeToggled(bool val);
    void drawInternalToggled(bool val);
    void singleColorChanged(const QColor &color);
    void singleColorOpacityChanged(int opacity);
    void multipleColorChanged(const QColor &color);
    void multipleColorOpacityChanged(int opacity);
    void colorModeChanged(int index);
    void colorTableClicked(bool useDefault, const QString &ctName);
    void subsetSelectionChanged();
    void overallOpacityChanged(int opacity);
    void smoothingLevelChanged(int index);

    void pointSizeChanged(double d);
    void pointSizePixelsChanged(int size);
    void pointTypeChanged(int index);
    void pointSizeVarToggled(bool on);
    void pointSizeVarChanged(const QString &);

private:
    int                     plotType;
    SubsetAttributes       *subsetAtts;
    QLabel                 *lineStyleLabel;
    QvisLineStyleWidget    *lineStyle;
    QLabel                 *lineWidthLabel;
    QvisLineWidthWidget    *lineWidth;
    QButtonGroup           *colorModeButtons;
    QCheckBox              *legendCheckBox;
    QGroupBox              *subsetColorGroup;
    QvisColorButton        *singleColor;
    QvisOpacitySlider      *singleColorOpacity;
    QLabel                 *multipleColorLabel;
    QListBox               *multipleColorList;
    QvisColorButton        *multipleColor;
    QvisOpacitySlider      *multipleColorOpacity;
    QvisColorTableButton   *colorTableButton;
    QvisOpacitySlider      *overallOpacity;
    QCheckBox              *wireframeCheckBox;
    QCheckBox              *drawInternalCheckBox;
    QButtonGroup           *smoothingLevelButtons;

    QvisPointControl       *pointControl;
};

#endif
