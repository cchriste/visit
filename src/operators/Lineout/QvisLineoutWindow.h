#ifndef QVISLINEOUTWINDOW_H
#define QVISLINEOUTWINDOW_H

#include <QvisOperatorWindow.h>
#include <AttributeSubject.h>

class LineoutAttributes;
class QCheckBox;
class QGroupBox;
class QLabel;
class QLineEdit;

// ****************************************************************************
// Class: QvisLineoutWindow
//
// Purpose: 
//   Defines QvisLineoutWindow class.
//
// Notes:      This class was automatically generated!

// Programmer: xml2window
// Creation:   Fri Nov 19 11:39:48 PDT 2004
//
// Modifications:
//   
// ****************************************************************************

class QvisLineoutWindow : public QvisOperatorWindow
{
    Q_OBJECT
  public:
    QvisLineoutWindow(const int type,
                         LineoutAttributes *subj,
                         const char *caption = 0,
                         const char *shortName = 0,
                         QvisNotepadArea *notepad = 0);
    virtual ~QvisLineoutWindow();
    virtual void CreateWindowContents();
  protected:
    void UpdateWindow(bool doAll);
    virtual void GetCurrentValues(int which_widget);
  private slots:
    void point1ProcessText();
    void point2ProcessText();
    void interactiveChanged(bool val);
    void ignoreGlobalChanged(bool val);
    void samplingOnChanged(bool val);
    void numberOfSamplePointsProcessText();
    void reflineLabelsChanged(bool val);
  private:
    QLineEdit *point1;
    QLineEdit *point2;
    QCheckBox *interactive;
    QGroupBox *ignoreGlobal;
    QCheckBox *samplingOn;
    QLineEdit *numberOfSamplePoints;
    QCheckBox *reflineLabels;
    QLabel *numberOfSamplePointsLabel;

    LineoutAttributes *atts;
};



#endif
