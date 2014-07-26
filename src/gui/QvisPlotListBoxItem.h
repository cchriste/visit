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

#ifndef QVIS_PLOT_LISTBOX_ITEM_H
#define QVIS_PLOT_LISTBOX_ITEM_H
#include <gui_exports.h>
#include <QListWidget>
#include <GUIBase.h>
#include <Plot.h>
#include <vector>

// ****************************************************************************
// Class: QvisPlotListBoxItem
//
// Purpose:
//   This is a list item that can be inserted into a QListWidget. The
//   PlotList that is displayed in the GUI is made of these.
//
// Notes:      
//
// Programmer: Brad Whitlock
// Creation:   Mon Sep 11 12:18:41 PDT 2000
//
// Modifications:
//   Brad Whitlock, Tue Apr 8 11:16:36 PDT 2003
//   I renamed the class to QvisPlotListBoxItem and made big changes to it
//   to allow expanded display of plots.
//
//   Brad Whitlock, Wed Jul 28 17:37:15 PST 2004
//   I added a method to return the prefix.
//
//   Brad Whitlock, Fri Feb 23 14:15:12 PST 2007
//   Added non const GetPlot method.
//
//   Brad Whitlock, Tue Jun 24 12:13:41 PDT 2008
//   Added GUIBase as a private base class so we can access the plugin managers.
//
//   Brad Whitlock, Fri Jul 23 14:55:01 PDT 2010
//   I added selectionName.
//
//   Brad Whitlock, Thu Mar 14 16:21:47 PDT 2013
//   Added applyOperators.
//
// ****************************************************************************
    
class GUI_API QvisPlotListBoxItem : public QListWidgetItem, private GUIBase
{
public:
    typedef enum {
        Expand,
        Subset,
        PlotAttributes,
        OperatorAttributes,
        Promote,
        Demote,
        Delete,
        ActiveOperator,
        Selection
    } ClickObjectType;

protected:
    class GUI_API ClickableRectangle
    {
    public:
        ClickableRectangle();
        ClickableRectangle(int, const QRect &, ClickObjectType);
        ClickableRectangle(const ClickableRectangle &obj);
        ~ClickableRectangle();
        void operator = (const ClickableRectangle &obj);

        QRect           r;
        ClickObjectType type;
        int             id;
    };

    typedef std::vector<ClickableRectangle> ClickableRectangleVector;
public:
    QvisPlotListBoxItem(const Plot &plot, const QString &prefix,  
                        const QString &selName);
   ~QvisPlotListBoxItem();

    int    height(const QListWidget *) const;
    int    width(const QListWidget *)  const;

    int    clicked(const QPoint &p, bool dc, int &id);
    bool   isExpanded() const;
    int    activeOperatorIndex() const;
    int    textX() const;

    void   setApplyOperators(bool);

    const Plot &GetPlot() const { return plot; };
    Plot &GetPlot() { return plot; };
    const QString &GetPrefix() const { return prefix; };
    const QString &GetSelectionName() const { return selectionName; }

    void paint(QPainter *);
private:
    QString GetDisplayString(const Plot &plot, const QString &prefix);
    void GetOperatorPixmap(int operatorType, QPixmap &pm);
    void GetPlotPixmap(int plotType, QPixmap &pm);
    void drawButtonSurface(QPainter *painter, const QRect &r) const;
    void drawUpDownButton(QPainter *painter, const QRect &r, bool up) const;
    void drawDeleteButton(QPainter *painter, const QRect &r) const;
    void AddClickableRectangle(int id, const QRect &r, ClickObjectType type);
    void setTextPen(QPainter *painter, bool highlightText) const;
private:
    Plot                     plot;
    QString                  prefix;
    QString                  selectionName;

    QPixmap                  plotIcon;
    QString                  plotName;
    QPixmap                 *operatorIcons;
    QString                 *operatorNames;
    int                      maxIconWidth, maxIconHeight;
    ClickableRectangleVector clickable;
    bool                     addClickableRects;
    bool                     applyOperators;

    static QPixmap          *subsetIcon; 
    static QPixmap          *selectionIcon;
};

#endif
