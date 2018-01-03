/*****************************************************************************
*
* Copyright (c) 2000 - 2017, Lawrence Livermore National Security, LLC
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

#ifndef QVISEXPLODEWINDOW_H
#define QVISEXPLODEWINDOW_H

#include <QvisOperatorWindow.h>
#include <AttributeSubject.h>
#include <MapNode.h>

class ExplodeAttributes;
class QCheckBox;
class QLineEdit;
class QSpinBox;
class QButtonGroup;
class QvisVariableButton;
class QGridLayout;
class QTabWidget;
class QComboBox;
class QGroupBox;

// ****************************************************************************
// Class: QvisExplodeWindow
//
// Purpose:
//    Defines QvisExplodeWindow class.
//
// Notes:      Autogenerated by xml2window.
//
// Programmer: xml2window
// Creation:   omitted
//
// Modifications:
//   
//   Alister Maguire, Tue Dec 12 13:13:01 PST 2017
//   Added explosionTabsChangedIndex(), ClearLayout(), 
//   ExplosionType, and materialCheckBoxes. Also 
//   removed boundary names.
//
// ****************************************************************************

class QvisExplodeWindow : public QvisOperatorWindow
{
    Q_OBJECT
  public:
    QvisExplodeWindow(const int type,
                         ExplodeAttributes *subj,
                         const QString &caption = QString::null,
                         const QString &shortName = QString::null,
                         QvisNotepadArea *notepad = 0);
    virtual ~QvisExplodeWindow();
    virtual void CreateWindowContents();
  protected:
    void UpdateWindow(bool doAll);
    virtual void GetCurrentValues(int which_widget);
    void ClearLayout(QLayout *layout);
  private slots:
    void explosionTypeChanged(int val);
    void explosionPointProcessText();
    void planePointProcessText();
    void planeNormProcessText();
    void cylinderPoint1ProcessText();
    void cylinderPoint2ProcessText();
    void materialExplosionFactorProcessText();
    void cellExplosionFactorProcessText();
    void materialProcessText();
    void cylinderRadiusProcessText();
    void explosionTabsChangedIndex(int type);
    void explosionPatternChangedIndex(int type);
    void explodeMaterialCellsToggled(bool val);
    void explodeAllCellsToggled(bool val);
  private:
    enum ExplosionTypes
    {
        POINT,
        PLANE, 
        CYLINDER
    };
    enum ExplosionPatterns
    {
        IMPACT,
        SCATTER
    };

    std::string  boundaryString;

    QTabWidget   *explosionTabs;
    QLineEdit    *materialExplosionFactor;
    QLineEdit    *cellExplosionFactor;
    QLineEdit    *explosionPoint;
    QLineEdit    *planePoint;
    QLineEdit    *planeNorm;
    QLineEdit    *cylinderPoint1;
    QLineEdit    *cylinderPoint2;
    QLineEdit    *cylinderRadius;
    QWidget      *materialsWidget;
    QComboBox    *materialsCombo;
    QCheckBox    *explodeMaterialCells;
    QCheckBox    *explodeAllCells;
    QGroupBox    *materialGroup;
    QComboBox    *explosionPattern;

    ExplodeAttributes *atts;
};



#endif
