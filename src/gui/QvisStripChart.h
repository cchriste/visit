
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
#ifndef QVIS_STRIPCHART_H
#define QVIS_STRIPCHART_H
#include <QVector>
#include <QWidget>

class QScrollArea;

// ****************************************************************************
// Class: VisItPointD
//
// Purpose:
//    Implements a simple point with a double as the base data type.
//    QT will have a QPointF in it's next release but for now I use this.
//
// Notes:
//
// Programmer: Shelly Prevost,
// Creation:   Wed Mar 21 16:35:30 PDT 2007.
//
// Modifications:
//
// ***************************************************************************
class VisItPointD
{
public:
           VisItPointD () { p_x = 0; p_y = 0; };
           VisItPointD ( double x, double y ) { p_x = x; p_y = y; };
    double x() { return p_x; };
    double y() { return p_y; };
    void   setX(double x) { p_x = x; };
    void   setY(double y) { p_y= y; };
private:
    double p_x;
    double p_y;  
};

typedef QVector<VisItPointD> Points;


// ****************************************************************************
// Class: VisItSimStripChart
//
// Purpose:
//    Implements a simple strip chart. It is intended 
//    for tracking rate of change of simulation variables.
//    The strip chart will auto range ( keep the full range
//    of data visible on the chart) and change the plot color
//    to red if the values are outside the set limits.
//
// Notes:      
//
// Programmer: Shelly Prevost, 
// Creation:   Friday Oct. 27, 2006
//
// Modifications:
//    Shelly Prevost, Wed Mar 21 16:35:30 PDT 2007.
//    Added Zoom In and out and focus ( center ) to controls
//
//    Shelly Prevost Fri Apr 13 14:03:03 PDT 2007
//    added Font variable to update font size as zoom changes. Also added variable
//    zoomOutLimit for zoom checks.
//
// ****************************************************************************

class VisItSimStripChart : public QWidget
{                                                          
    Q_OBJECT
    
public:
            VisItSimStripChart( QWidget *parent=0, int winX=4000, int winY=1000 );
            ~VisItSimStripChart();
    void    setOutOfBandLimits( double minY, double miaxY );
    void    getOutOfBandLimits( double &minY, double &maxY );
    void    enableOutOfBandLimits( bool enable);
    bool    getEnableOutOfBandLimits();
    bool    addDataPoint( double x, double y );
    void    setEnable( bool enable );
    bool    getEnable();
    void    getMinMaxData( double &minY, double &maxY);
    void    zoomOut();
    void    zoomIn();
    void    reset();
    void    focus(QScrollArea *sc);
    void    setFontSize();
    void    setEnableLogScale( bool enable );
    bool    getEnableLogScale();
    double  getCurrentData();
    int     getCurrentCycle();
       
protected:
    void    paintEvent( QPaintEvent * );
    void    mousePressEvent( QMouseEvent *);
    void    mouseReleaseEvent( QMouseEvent *);
    void    mouseMoveEvent( QMouseEvent *);
    void    paintGrid(QPainter *paint);

private:
    Points  points;           // point array
    int     timeShift;    // how far left to points
    bool    down;         // TRUE if mouse down
    float   delta;
    float   vdelta;
    float   middle;
    float   maxPoint;
    float   minPoint;
    double  minYLimit;
    double  maxYLimit;
    double  minData;
    double  maxData;
    double  currentData;
    int     currentCycle;
    int     winXSize;

    int     winYSize;
    bool    enabled;
    bool    enableLogScale;
    bool    outOfBandLimitsEnabled;
    float   zoom;
    bool    center;
    float   zoomOutLimit;
    QFont   *gridFont;
    int     pointSize;
    int     currentScaledY;
};
#endif /* QVIS_STRIPCHART_MGR */
