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

#ifndef QVIS_VECTOR_WINDOW_H
#define QVIS_VECTOR_WINDOW_H
#include <QvisPostableWindowObserver.h>

// Forward declarations
class QButtonGroup;
class QCheckBox;
class QComboBox;
class QGroupBox;
class QLabel;
class QLineEdit;
class QvisColorButton;
class QvisColorTableButton;
class QvisLineStyleWidget;
class QvisLineWidthWidget;
class QvisOpacitySlider;
class VectorAttributes;

// ****************************************************************************
// Class: QvisVectorPlotWindow
//
// Purpose:
//   This class is a postable window that watches vector plot attributes and
//   always represents their current state.
//
// Notes:      
//
// Programmer: Hank Childs & Brad Whitlock
// Creation:   Thu Mar 22 23:40:52 PST 2001
//
// Modifications:
//   Brad Whitlock, Sat Jun 16 18:21:34 PST 2001
//   I added color table stuff.
//
//   Brad Whitlock, Fri Aug 29 11:34:52 PDT 2003
//   I grouped related fields into group boxes.
//
//   Jeremy Meredith, Fri Nov 21 12:06:46 PST 2003
//   Added radio buttons for vector origin.
//
//   Eric Brugger, Tue Nov 23 10:13:05 PST 2004
//   Added scaleByMagnitude and autoScale.
//
//   Kathleen Bonnell, Wed Dec 22 16:42:35 PST 2004 
//   Added widgets for min/max and limitsSelection.
//
// ****************************************************************************

class QvisVectorPlotWindow : public QvisPostableWindowObserver
{
    Q_OBJECT
public:
    QvisVectorPlotWindow(const int type, VectorAttributes *_vecAtts,
                         const char *caption = 0,
                         const char *shortName = 0,
                         QvisNotepadArea *notepad = 0);
    virtual ~QvisVectorPlotWindow();
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
    void vectorColorChanged(const QColor &color);
    void processScaleText();
    void scaleByMagnitudeToggled();
    void autoScaleToggled();
    void processHeadSizeText();
    void reduceMethodChanged(int index);
    void processNVectorsText();
    void processStrideText();
    void legendToggled();
    void drawHeadToggled();
    void colorModeChanged(int);
    void colorTableClicked(bool useDefault, const QString &ctName);
    void originTypeChanged(int);

    void minToggled(bool on);
    void maxToggled(bool on);
    void processMaxLimitText();
    void processMinLimitText();
    void limitsSelectChanged(int);


private:
    int                  plotType;
    VectorAttributes     *vectorAtts;

    QvisLineStyleWidget  *lineStyle;
    QvisLineWidthWidget  *lineWidth;
    QGroupBox            *colorGroupBox;
    QvisColorButton      *vectorColor;
    QButtonGroup         *colorButtonGroup; 
    QvisColorTableButton *colorTableButton;
    QGroupBox            *scaleGroupBox;    
    QLineEdit            *scaleLineEdit;
    QCheckBox            *scaleByMagnitudeToggle;
    QCheckBox            *autoScaleToggle;
    QLineEdit            *headSizeLineEdit;
    QGroupBox            *reduceGroupBox;
    QButtonGroup         *reduceButtonGroup;
    QLineEdit            *nVectorsLineEdit;
    QLineEdit            *strideLineEdit;
    QCheckBox            *legendToggle;
    QCheckBox            *drawHeadToggle;
    QButtonGroup         *originButtonGroup;

    QGroupBox             *limitsGroupBox;
    QCheckBox             *minToggle;
    QCheckBox             *maxToggle;
    QComboBox             *limitsSelect;
    QLineEdit             *maxLineEdit;
    QLineEdit             *minLineEdit;
};

#endif
