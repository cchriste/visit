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

#ifndef SILODIRVIEW_H
#define SILODIRVIEW_H

#include <qlistview.h>
class SiloDir;

// ****************************************************************************
//  Class:  SiloDirViewItem
//
//  Purpose:
//    An item for a SiloDirView.
//
//  Programmer:  Jeremy Meredith
//  Creation:    November 12, 2001
//
// ****************************************************************************
class SiloDirViewItem : public QListViewItem
{
  public:
    SiloDir *dir;
    QString  name;

    SiloDirViewItem(SiloDir *d, QListView *lv, const QString &n)
        : QListViewItem(lv,n), dir(d), name(n) { };
    SiloDirViewItem(SiloDir *d, QListViewItem *li, const QString &n) 
        : QListViewItem(li,n), dir(d), name(n) { };
};

// ****************************************************************************
//  Class:  SiloDirView
//
//  Purpose:
//    ListView-style container for viewing the TOC of a Silo directory.
//
//  Programmer:  Jeremy Meredith
//  Creation:    November 12, 2001
//
// ****************************************************************************
class SiloDirView : public QListView
{
    Q_OBJECT
  public:
    SiloDirView(QWidget *p, const QString &n);

    void Set(SiloDir *d);

    virtual QSize sizeHint() const;

  public slots:
    void ChangeDir(QListViewItem *i);

  protected:
    virtual void resizeEvent(QResizeEvent *re);

  private:
    int total_items;
    QPixmap *mesh_pixmap;
    QPixmap *mat_pixmap;
    QPixmap *var_pixmap;
    QPixmap *obj_pixmap;
    QPixmap *spec_pixmap;
    QPixmap *curve_pixmap;
    QPixmap *array_pixmap;
    QPixmap *silovar_pixmap;
};

#endif
