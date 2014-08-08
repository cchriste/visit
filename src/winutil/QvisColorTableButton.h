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

#ifndef QVIS_COLORTABLE_BUTTON_H
#define QVIS_COLORTABLE_BUTTON_H
#include <winutil_exports.h>
#include <QPushButton>

// Forward declarations.
class QAction;
class QActionGroup;
class QMenu;
class ColorTableAttributes;

// ****************************************************************************
// Class: QvisColorTableButton
//
// Purpose:
//   This is a type of push button that is aware of the different color tables
//   that can be used for plots.
//
// Notes:      
//
// Programmer: Brad Whitlock
// Creation:   Wed Jun 20 12:37:24 PDT 2001
//
// Modifications:
//   Brad Whitlock, Thu Feb 14 13:37:41 PST 2002
//   Added a counter.
//
//   Brad Whitlock, Tue Feb 20 11:47:37 PDT 2007
//   Changed API.
//
//   Brad Whitlock, Fri May  9 11:20:10 PDT 2008
//   Qt 4.
//
//   Brad Whitlock, Wed Apr 25 16:06:56 PDT 2012
//   Add color table icons.
//
//   Kathleen Biagas, Mon Aug  4 15:51:11 PDT 2014
//   Change colorTableNames to a QSringList, add mappedColorTableNames
//   to aid in grouping.  Add category argument to addColorTable.
//
// ****************************************************************************

class WINUTIL_API QvisColorTableButton : public QPushButton
{
    Q_OBJECT

    typedef std::vector<QvisColorTableButton *> ColorTableButtonVector;
public:
    QvisColorTableButton(QWidget *parent);
    virtual ~QvisColorTableButton();
    virtual QSize sizeHint() const;
    virtual QSizePolicy sizePolicy () const;

    void setColorTable(const QString &ctName);
    const QString &getColorTable() const;
    void useDefaultColorTable();

    // Methods to set the list of internal color tables.
    static void clearAllColorTables();
    static void addColorTable(const QString &ctName, const QString &ctCategory);
    static void updateColorTableButtons();
    static void setColorTableAttributes(ColorTableAttributes *cAtts);
signals:
    void selectedColorTable(bool useDefault, const QString &ctName);
private slots:
    void popupPressed();
    void colorTableSelected(QAction *);
private:
    static int  getColorTableIndex(const QString &ctName);
    static void regeneratePopupMenu();
    static QIcon getIcon(const QString &);
    static QIcon makeIcon(const QString &);

    QString                        colorTable;

    static int                     numInstances;
    static QMenu                  *colorTableMenu;
    static QActionGroup           *colorTableMenuActionGroup;
    static bool                    popupHasEntries;
    static ColorTableButtonVector  buttons;

    static QStringList             colorTableNames;
    static QMap<QString, QStringList> mappedColorTableNames;
    static ColorTableAttributes   *colorTableAtts;
};

#endif
