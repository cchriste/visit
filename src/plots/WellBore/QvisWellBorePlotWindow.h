/*****************************************************************************
*
* Copyright (c) 2000 - 2008, Lawrence Livermore National Security, LLC
* Produced at the Lawrence Livermore National Laboratory
* LLNL-CODE-400142
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

#ifndef QVISWELLBOREPLOTWINDOW_H
#define QVISWELLBOREPLOTWINDOW_H

#include <QvisPostableWindowObserver.h>
#include <AttributeSubject.h>

class WellBoreAttributes;
class QButtonGroup;
class QCheckBox;
class QComboBox;
class QGroupBox;
class QLabel;
class QLineEdit;
class QListBox;
class QMultiLineEdit;
class QvisColorButton;
class QvisColorManagerWidget;
class QvisColorTableButton;
class QvisLineStyleWidget;
class QvisLineWidthWidget;
class QvisOpacitySlider;

// ****************************************************************************
// Class: QvisWellBorePlotWindow
//
// Purpose:
//    Defines QvisWellBorePlotWindow class.
//
// Notes:      Autogenerated by xml2window.
//
// Programmer: xml2window
// Creation:   omitted
//
// Modifications:
//     Eric Brugger, Mon Nov 10 13:16:16 PST 2008
//     Added the ability to display well bore names and stems.
//   
// ****************************************************************************

class QvisWellBorePlotWindow : public QvisPostableWindowObserver
{
    Q_OBJECT
  public:
    QvisWellBorePlotWindow(const int type,
                         WellBoreAttributes *subj,
                         const QString &caption = QString::null,
                         const QString &shortName = QString::null,
                         QvisNotepadArea *notepad = 0);
    virtual ~QvisWellBorePlotWindow();
    virtual void CreateWindowContents();
  public slots:
    virtual void apply();
    virtual void makeDefault();
    virtual void reset();
  protected:
    void UpdateWindow(bool doAll);
    void GetCurrentValues(int which_widget);
    void Apply(bool ignore = false);
    void BlockAllSignals(bool block);

    bool UpdateMultipleAreaColors();
    bool UpdateMultipleAreaNames();
    void UpdateWellName(int index);
    void UpdateWellDefinition(int index);

    void ReadWellDefinition();
  private slots:
    void readWellBoresButtonPressed();
    void writeWellBoresButtonPressed();
    void wellListSelectionChanged();
    void newWellButtonPressed();
    void deleteWellButtonPressed();
    void wellNameTextChanged(const QString &color);
    void wellDefinitionTextChanged();
    void colorModeChanged(int index);
    void singleColorChanged(const QColor &color);
    void singleColorOpacityChanged(int opacity);
    void multipleColorChanged(const QColor &color, int index);
    void opacityChanged(int opacity, int index);
    void colorTableClicked(bool useDefault, const QString &ctName);
    void drawWellsAsChanged(int val);
    void wellCylinderQualityChanged(int val);
    void wellRadiusProcessText();
    void wellLineWidthChanged(int style);
    void wellLineStyleChanged(int style);
    void wellAnnotationChanged(int val);
    void wellStemHeightProcessText();
    void wellNameScaleProcessText();
    void legendFlagChanged(bool val);
  private:
    int                     plotType;
    bool                    wellDefinitionChanged;
    int                     wellIndex;

    QPushButton            *readWellBoresButton;
    QPushButton            *writeWellBoresButton;
    QListBox               *wellListBox;
    QPushButton            *newWellButton;
    QPushButton            *deleteWellButton;
    QLabel                 *wellNameLabel;
    QLineEdit              *wellName;
    QLabel                 *wellDefinitionLabel;
    QMultiLineEdit         *wellDefinition;
    QGroupBox              *wellColorGroup;
    QButtonGroup           *colorModeButtons;
    QvisColorButton        *singleColor;
    QvisOpacitySlider      *singleColorOpacity;
    QvisColorManagerWidget *multipleColors;
    QvisColorTableButton   *colorTableButton;
    QLabel                 *drawWellsAsLabel;
    QComboBox              *drawWellsAs;
    QLabel                 *wellCylinderQualityLabel;
    QComboBox              *wellCylinderQuality;
    QLabel                 *wellRadiusLabel;
    QLineEdit              *wellRadius;
    QLabel                 *wellLineWidthLabel;
    QvisLineWidthWidget    *wellLineWidth;
    QLabel                 *wellLineStyleLabel;
    QvisLineStyleWidget    *wellLineStyle;
    QLabel                 *wellAnnotationLabel;
    QComboBox              *wellAnnotation;
    QLabel                 *wellStemHeightLabel;
    QLineEdit              *wellStemHeight;
    QLabel                 *wellNameScaleLabel;
    QLineEdit              *wellNameScale;
    QCheckBox              *legendFlag;

    WellBoreAttributes     *atts;

    // Routines for writing well bore files.
    bool                    GetPoint(int [3], const std::vector<int> &,
                                int &);
    bool                    WritePoint(FILE *, int [3], int [3]);

    // Routines for writing well bore files.
    enum TokenType
    {
        TOKEN_EOF,
        TOKEN_EQUAL,
        TOKEN_STRING,
        TOKEN_SLASH,
        TOKEN_INTEGER,
        TOKEN_WELLS,
        TOKEN_ENDWELLS,
        TOKEN_WELL,
        TOKEN_PERF,
        TOKEN_NAME,
        TOKEN_I,
        TOKEN_J,
        TOKEN_K,
        TOKEN_I_BOTTOM,
        TOKEN_J_BOTTOM,
        TOKEN_K_BOTTOM,
        TOKEN_I_TOP,
        TOKEN_J_TOP,
        TOKEN_K_TOP,
        TOKEN_ERROR
    };

    char                    buf[101];
    TokenType               GetToken(FILE *);
};



#endif
