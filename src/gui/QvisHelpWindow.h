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

#ifndef QVIS_HELP_WINDOW_H
#define QVIS_HELP_WINDOW_H
#include <gui_exports.h>
#include <qmap.h>
#include <qpixmap.h>
#include <QvisDelayedWindow.h>

class QAction;
class QDomElement;
class QLineEdit;
class QListBox;
class QListView;
class QListViewItem;
class QPushButton;
class QSplitter;
class QTabWidget;
class QTextBrowser;
class QVBox;
class QvisHelpListViewItem;

// ****************************************************************************
// Class: QvisHelpWindow
//
// Purpose:
//   This class creates a help window.
//
// Notes:      
//
// Programmer: Brad Whitlock
// Creation:   Wed May 15 12:00:16 PDT 2002
//
// Modifications:
//   Brad Whitlock, Wed Jul 10 17:51:16 PST 2002
//   Finished it.
//
//   Brad Whitlock, Tue Sep 10 16:23:09 PST 2002
//   Added an internal convenience method.
//
//   Brad Whitlock, Thu Feb 17 12:14:33 PDT 2005
//   Added synchronizeContents.
//
// ****************************************************************************

class GUI_API QvisHelpWindow : public QvisDelayedWindow
{
    Q_OBJECT
public:
    QvisHelpWindow(const char *captionString);
    virtual ~QvisHelpWindow();

    virtual void CreateWindowContents();
    virtual void CreateNode(DataNode *);
    virtual void SetFromNode(DataNode *, const int *borders);
public slots:
    void displayCopyright();
    void displayReleaseNotes();
    void displayReleaseNotesIfAvailable();
    virtual void show();
private slots:
    void activeTabChanged(QWidget *);
    void activateContentsTab();
    void activateIndexTab();
    void activateBookmarkTab();
    void openHelp(QListViewItem *);
    void topicExpanded(QListViewItem *);
    void topicCollapsed(QListViewItem *);
    void displayNoHelp();
    void displayTitle(const QString &title);
    void displayHome();
    bool displayPage(const QString &page, bool reload = false);
    void increaseFontSize();
    void decreaseFontSize();
    void displayIndexTopic();
    void lookForIndexTopic(const QString &topic);
    void displayBookmarkTopic();
    void addBookmark();
    void removeBookmark();
private:
    typedef QMap<QString, QString> IndexMap;

    void LoadHelp(const QString &helpFile);
    void BuildIndex();
    void AddToIndex(const QString &topic, const QString &doc);
    void BuildContents(QListViewItem *parentItem,
                       const QDomElement &parentElement);
    void BuildBookmarks();
    QString TopicFromDoc(const QString &doc);
    bool TopicFromDocHelper(QString &str, const QString &doc,
                            QvisHelpListViewItem *item);
    QString CompleteFileName(const QString &page) const;
    void synchronizeContents(const QString &page);
    void displayReleaseNotesHelper(bool);

    QTabWidget   *helpTabs;
    QListView    *helpContents;
    QTextBrowser *helpBrowser;
    QSplitter    *splitter;
    QAction      *backAction;
    QAction      *forwardAction;

    QVBox        *helpIndexTab;
    QListBox     *helpIndex;
    QLineEdit    *helpIndexText;

    QVBox        *helpBookmarksTab;
    QPushButton  *addBookmarkButton;
    QPushButton  *removeBookmarkButton;
    QListBox     *helpBookMarks;

    QPixmap      closedBookIcon;
    QPixmap      openBookIcon;
    QPixmap      helpIcon;
    QString      helpFile;
    QString      helpPath;
    bool         firstShow;
    int          activeTab;
    IndexMap     index;
    IndexMap     bookmarks;
};

#endif
