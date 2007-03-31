#ifndef QVISCURVEPLOTWINDOW_H
#define QVISCURVEPLOTWINDOW_H

#include <QvisPostableWindowObserver.h>
#include <AttributeSubject.h>

class CurveAttributes;
class QLabel;
class QNarrowLineEdit;
class QCheckBox;
class QButtonGroup;
class QvisColorButton;
class QvisLineStyleWidget;
class QvisLineWidthWidget;

// ****************************************************************************
// Class: QvisCurvePlotWindow
//
// Purpose: 
//   Defines QvisCurvePlotWindow class.
//
// Notes:      This class was automatically generated!

// Programmer: xml2window
// Creation:   Tue Jul 23 13:34:33 PST 2002
//
// Modifications:
//   Kathleen Bonnell, Tue Dec 23 13:27:22 PST 2003
//   Added showPoints, pointSize and pointSizeLabel widgets. 
//   Added showPointChanged and processPointSizeText slots.
//   
// ****************************************************************************

class QvisCurvePlotWindow : public QvisPostableWindowObserver
{
    Q_OBJECT
  public:
    QvisCurvePlotWindow(const int type,
                         CurveAttributes *subj,
                         const char *caption = 0,
                         const char *shortName = 0,
                         QvisNotepadArea *notepad = 0);
    virtual ~QvisCurvePlotWindow();
    virtual void CreateWindowContents();
  public slots:
    virtual void apply();
    virtual void makeDefault();
    virtual void reset();
  protected:
    void UpdateWindow(bool doAll);
    void GetCurrentValues(int which_widget);
    void Apply(bool ignore = false);
  private slots:
    void lineStyleChanged(int style);
    void lineWidthChanged(int style);
    void colorChanged(const QColor &color);
    void showLabelsChanged(bool val);
    void showPointsChanged(bool val);
    void processPointSizeText();
  private:
    int plotType;
    QvisLineStyleWidget *lineStyle;
    QvisLineWidthWidget *lineWidth;
    QvisColorButton     *color;
    QCheckBox           *showLabels;
    QCheckBox           *showPoints;
    QNarrowLineEdit     *pointSize;
    QLabel              *pointSizeLabel;

    CurveAttributes *atts;
};



#endif
