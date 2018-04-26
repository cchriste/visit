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

#ifndef QVISRELATVIEWWINDOW_H
#define QVISRELATVIEWWINDOW_H

#include <QvisOperatorWindow.h>
#include <AttributeSubject.h>
#include <map>

#define VARS 0
#define TUPS 1
#define SOME 2

class ModelFitAtts;
class QvisVariableButton;
class QGroupBox;
class QRadioButton;
 
class QTableWidget;

class tableEntry 
{
  protected:
    int row;
    int col;
    std::string contents;
  
  public:
    tableEntry(int inRow, int inCol, std::string inString)
    {
        setRow(inRow);
        setCol(inCol);
        setContents(inString);
    }
    int getRow()
    {
        return row;
    }
    int getCol()
    {
        return col;
    }
    std::string getContents()
    {
        return contents;
    }
    void setRow(int inRow)
    {
        row = inRow;
    }
    void setCol(int inCol)
    {
        col = inCol;
    }
    void setContents(std::string inContents)
    {
        contents = inContents;
    }
};

// ****************************************************************************
// Class: QvisModelFitWindow
//
// Purpose:
//    Defines QvisModelFitWindow class.
//
// Notes:      Autogenerated by xml2window.
//
// Programmer: xml2window
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

class QvisModelFitWindow : public QvisOperatorWindow
{
    Q_OBJECT
  public:
    QvisModelFitWindow(const int type,
                         ModelFitAtts *subj,
                         const QString &caption = QString::null,
                         const QString &shortName = QString::null,
                         QvisNotepadArea *notepad = 0);
    virtual ~QvisModelFitWindow();
    virtual void CreateWindowContents();
  protected:
    void UpdateWindow(bool doAll);
    virtual void GetCurrentValues(int which_widget);
  private slots:
    void addRelationship();
    //void addSpaces();
    void varAdded(const QString&);
    void addTuple();
    void deleteVariable();
    void deletePoint();
    void deleteRelationship();
    void fillTables();
    void selectTypeChanged(int buttonID);
    void inputSpaceChanged(int buttonID);
    void distanceTypeChanged(int buttonID);
    void prepareTable(int, int, int, int);
    void storeTableContents(int, int);
    void storeModelNames(int, int);

  private:
    QPushButton *addRelat, *addTup, *deleteMod,  *deleteVar, *deleteTup;
    QvisVariableButton *addVar;
    //QvisVariableButton *addModel;
    QRadioButton *vSpace,    *nSpace,   *lSpace,    *pSpace,    *varSpace, 
                 *normSpace, *logSpace, *probSpace, *euclidean, *manhattan, *chebyshev;
    QStringList vHeaderLbls;
    QStringList hHeaderLbls;
    QGroupBox *distanceBox;
    QGroupBox *scoreBox;

    bool inPrepareTable;
    int cur_column;
    int cur_row;
    int tuple_number;
    int num_relats;
    intVector modelNumbers;
    intVector selection_type;
    intVector input_space;
    intVector distance_type;
    intVector numVars;
    intVector numPoints;
    stringVector modelNames;
    std::multimap<int, std::string>guiVarNames;
    std::multimap<int, tableEntry>theTable;
    int score_norm_type;

    QTableWidget *ModelFit;
    QTableWidget *models;

    ModelFitAtts *atts;
};



#endif
