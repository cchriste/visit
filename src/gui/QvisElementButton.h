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

#ifndef QVIS_ELEMENT_BUTTON_H
#define QVIS_ELEMENT_BUTTON_H
#include <gui_exports.h>
#include <QColor>
#include <QPushButton>
#include <vector>

class QPainter;
class QMenu;
class QvisPeriodicTableWidget;
class QvisElementSelectionWidget;

// ****************************************************************************
// Class: QvisElementButton
//
// Purpose:
//   This class represents a color button widget that can be used to select
//   colors for materials, isocontours, and other items in the gui.
//
// Notes:      
//
// Programmer: Jeremy Meredith
// Creation:   August 29, 2006
//
// Notes: Taken largely from QvisColorButton
//
// Modifications:
//    Jeremy Meredith, Mon Feb 11 16:46:57 EST 2008
//    Changed to use the element selection widget instead of the
//    simple periodic table widget; the former was created to contain
//    both a periodic table widget and a "match any element" button
//    to allow wildcards.
//   
//    Jeremy Meredith, Tue Feb 12 14:01:52 EST 2008
//    Added support for hinting some elements to the user, e.g. to highlight
//    the elements that are actually in the database.
//
//    Brad Whitlock, Tue Jun  3 14:43:22 PDT 2008
//    Qt 4.
//
// ****************************************************************************

class GUI_API QvisElementButton : public QPushButton
{
    Q_OBJECT

    typedef std::vector<QvisElementButton *> ElementButtonVector;
public:
    QvisElementButton(QWidget *parent = 0, const void *userData = 0);
    virtual ~QvisElementButton();
    virtual QSize sizeHint() const;
    virtual QSizePolicy sizePolicy () const;

    int elementNumber() const;
    void setElementNumber(int);
    void setHintedElements(const std::vector<int>&);

signals:
    void selectedElement(int element);
protected:
private slots:
    void popupPressed();
    void elementSelected(int element);
private:
    int                                number;
    const void                        *userData;

    static QvisElementSelectionWidget *sharedpopup;
    static ElementButtonVector         buttons;
};

#endif
