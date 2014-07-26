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

#ifndef QVIS_SIL_SET_SELECTOR_H
#define QVIS_SIL_SET_SELECTOR_H

#include <gui_exports.h>
#include <QWidget>
#include <SimpleObserver.h>
#include <GUIBase.h>
#include <vectortypes.h>

class QComboBox;
class QLabel;
class SILRestrictionAttributes;

// ****************************************************************************
// Class: QvisSILSetSelector
//
// Purpose: 
//   Defines QvisSILSetSelector class.
//
// Programmer: Kathleen Bonnell 
// Creation:   June 6, 2007 
//
// Modifications:
//   Kathleen Bonnell, Thu Jun 14 12:18:47 PDT 2007
//   Added userCategory, userSubset, so that options stored in sesisonfiles
//   can be restored.
//
//   Brad Whitlock, Fri Jul 18 08:35:26 PDT 2008
//   Qt 4.
//
// ****************************************************************************

class GUI_API QvisSILSetSelector : public QWidget, 
                                   public SimpleObserver, 
                                   public GUIBase
{
    Q_OBJECT
  public:
    QvisSILSetSelector(QWidget *parent,
            SILRestrictionAttributes *, intVector &);
    virtual ~QvisSILSetSelector();

    virtual void Update(Subject *);
    virtual void SubjectRemoved(Subject *);

    void SetCategoryName(const QString &name);
    QString GetCategoryName() const;
    void SetSubsetName(const QString &name);
    QString GetSubsetName() const;


  signals:
    void categoryChanged(const QString &);
    void subsetChanged(const QString &);

  private slots:
    void categoryNameChanged();
    void subsetNameChanged();

  private:
    void UpdateComboBoxes();
    void FillCategoryBox();
    void FillSubsetBox();

    QLabel    *categoryLabel;
    QComboBox *categoryName;
    QLabel    *subsetLabel;
    QComboBox *subsetName;

    SILRestrictionAttributes *silAtts;
    QString defaultItem;
    QString lastGoodCategory;
    QString lastGoodSubset;
    QString userCategory;
    QString userSubset;
    int silTopSet;
    int silNumSets;
    int silNumCollections;
    unsignedCharVector silUseSet;
    intVector allowedCategories;
};


#endif
