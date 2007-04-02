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

#ifndef QVIS_NOTEPAD_AREA
#define QVIS_NOTEPAD_AREA
#include <gui_exports.h>
#include <qvbox.h>
#include <qsizepolicy.h>
#include <qmap.h>

// forward declarations
class QvisPostableWindow;
class QTabWidget;

// ****************************************************************************
// Class: QvisNotepadArea
//
// Purpose:
//   This class allows QvisPostableWindows to be posted into it.
//
// Notes:      
//
// Programmer: Brad Whitlock
// Creation:   Tue Jul 25 17:56:27 PST 2000
//
// Modifications:
//   Brad Whitlock, Thu Sep 6 15:15:59 PST 2001
//   Added a method to handle paletteChange events.
//
//   Brad Whitlock, Tue Sep 25 15:16:59 PST 2001
//   Changed an internal data structure and removed the method to handle
//   paletteChange events since I found a better way to do it.
//
// ****************************************************************************

class GUI_API QvisNotepadArea : public QVBox
{
    Q_OBJECT

    typedef struct
    {
        QWidget *parent;
        bool     parentIsScrollView;
        int      minWidth;
        int      minHeight;
    } PostedInfo;

    typedef QMap<QWidget*, PostedInfo> PostedInfoLookup;
public:
    QvisNotepadArea(QWidget *parent = 0, const char *name = 0);
    virtual ~QvisNotepadArea();
    void showPage(QvisPostableWindow *pw);
    void postWindow(QvisPostableWindow *pw);
private:
    int              numPosted;
    QTabWidget       *tabs;
    QWidget          *empty;
    PostedInfoLookup postedLookup;
};

#endif
