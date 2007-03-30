#ifndef QVIS_CLIP_WINDOW_H
#define QVIS_CLIP_WINDOW_H
#include <QvisOperatorWindow.h>

class QCheckBox;
class QComboBox;
class QButtonGroup;
class QLineEdit;
class QLabel;
class QTabWidget;
class QVBox;
class ClipAttributes;

// ****************************************************************************
// Class: QvisClipWindow
//
// Purpose:
//   This class is a postable window that watches clip operator
//   attributes and always represents their current state.
//
// Notes:      
//
// Programmer: Kathleen Bonnell 
// Creation:   May 7, 2001 
//
// Modifications:
//   Brad Whitlock, Fri Apr 12 13:04:35 PST 2002
//   Made it inherit from QvisOperatorWindow.
//
// ****************************************************************************

class QvisClipWindow : public QvisOperatorWindow
{
    Q_OBJECT
public:
    QvisClipWindow(const int type,
                    ClipAttributes *subj,
                    const char *caption = 0,
                    const char *shortName = 0,
                    QvisNotepadArea *notepad = 0);
    virtual ~QvisClipWindow();
protected:
    virtual void CreateWindowContents();
    void UpdateWindow(bool doAll);
    virtual void GetCurrentValues(int which_widget);
private slots:
    void processPlane1Origin();
    void processPlane2Origin();
    void processPlane3Origin();
    void processPlane1Normal();
    void processPlane2Normal();
    void processPlane3Normal();
    void processCenterText();
    void processRadiusText();
    void planeInverseToggled(bool);
    void sphereInverseToggled(bool);
    void plane1StatusClicked(int);
    void plane2StatusClicked(int);
    void plane3StatusClicked(int);
    void tabWidgetChanged(QWidget *);
private:
    QLineEdit    *plane1Origin;
    QLineEdit    *plane2Origin;
    QLineEdit    *plane3Origin;
    QLineEdit    *plane1Normal;
    QLineEdit    *plane2Normal;
    QLineEdit    *plane3Normal;
    QButtonGroup *plane1Status;
    QButtonGroup *plane2Status;
    QButtonGroup *plane3Status;
    QLineEdit    *centerLineEdit;
    QLineEdit    *radiusLineEdit;
    QCheckBox    *planeInverse;
    QCheckBox    *sphereInverse;
    QTabWidget   *tabWidget;
    QVBox        *planeBox;
    QVBox        *sphereBox;

    ClipAttributes *clipAtts;
};
#endif
