#ifndef QVIS_VIEW_WINDOW_H
#define QVIS_VIEW_WINDOW_H
#include <gui_exports.h>
#include <QvisPostableWindowSimpleObserver.h>

// Forward declarations.
class DataNode;
class QButtonGroup;
class QCheckBox;
class QComboBox;
class QGroupBox;
class QLabel;
class QLineEdit;
class QRadioButton;
class QSlider;
class QTabWidget;
class QVBox;
class ViewCurveAttributes;
class View2DAttributes;
class View3DAttributes;
class WindowInformation;
class QPushButton;

// ****************************************************************************
// Class: QvisViewWindow
//
// Purpose:
//   This class implements the VisIt view window.
//
// Notes:      
//
// Programmer: Brad Whitlock
// Creation:   Thu Jul 26 16:03:25 PST 2001
//
// Modifications:
//   Eric Brugger, Tue Aug 21 13:45:36 PDT 2001
//   I redesigned the window fairly extensively.
//
//   Brad Whitlock, Mon Aug 27 11:35:11 PDT 2001
//   I made the window postable and added some default view buttons for 3d.
//
//   Eric Brugger, Tue Aug  6 10:30:24 PDT 2002
//   I added a view command processor.
//
//   Brad Whitlock, Tue Sep 17 12:58:39 PDT 2002
//   I reorganized the window.
//
//   Eric Brugger, Mon Jan 13 14:50:38 PST 2003
//   Add the vp and wp commands to the cli.
//
//   Jeremy Meredith, Tue Feb  4 17:44:53 PST 2003
//   Added controls for the camera and view keyframes.
//
//   Eric Brugger, Fri Apr 18 11:45:10 PDT 2003
//   Removed auto center view.
//
//   Eric Brugger, Tue Jun 10 12:48:03 PDT 2003
//   I added image pan and image zoom fields to the 3d view.
//
//   Eric Brugger, Wed Aug 20 14:02:49 PDT 2003
//   I added support for curve views and split the view attributes into
//   2d and 3d parts.
//
//   Brad Whitlock, Thu Sep 11 09:33:08 PDT 2003
//   I added slots to reset or recenter the view.
//
//   Eric Brugger, Thu Oct 16 12:21:13 PDT 2003
//   I added full frame mode to the 2D view tab.
//
//   Eric Brugger, Tue Feb 10 10:29:21 PST 2004
//   I added center of rotation controls to the advanced tab.
//
//   Mark C. Miller, Thu Jul 21 12:52:42 PDT 2005
//   Added stuff for auto full frame mode
// ****************************************************************************

class GUI_API QvisViewWindow : public QvisPostableWindowSimpleObserver
{
    Q_OBJECT
public:
    QvisViewWindow(const char *caption = 0, const char *shortName = 0,
                   QvisNotepadArea *notepad = 0);
    virtual ~QvisViewWindow();
    virtual void CreateWindowContents();
    void SubjectRemoved(Subject *TheRemovedSubject);

    void ConnectCurveAttributes(ViewCurveAttributes *v);
    void Connect2DAttributes(View2DAttributes *v);
    void Connect3DAttributes(View3DAttributes *v);
    void ConnectWindowInformation(WindowInformation *);

    virtual void CreateNode(DataNode *parentNode);
    virtual void SetFromNode(DataNode *parentNode, const int *borders);
public slots:
    virtual void apply();
    virtual void show();
protected:
    void Apply(bool ignore = false);
    void GetCurrentValues(int which_widget);
    void GetCurrentValuesCurve(int which_widget);
    void GetCurrentValues2d(int which_widget);
    void GetCurrentValues3d(int which_widget);

    virtual void UpdateWindow(bool doAll);
    void UpdateCurve(bool doAll);
    void Update2D(bool doAll);
    void Update3D(bool doAll);
    void UpdateGlobal(bool doAll);
private slots:
    void processCommandText();

    void processViewportCurveText();
    void processDomainText();
    void processRangeText();

    void processViewportText();
    void processWindowText();
    void fullFrameActivationModeChanged(int);

    void processNormalText();
    void processFocusText();
    void processUpVectorText();
    void processViewAngleText();
    void processParallelScaleText();
    void processNearText();
    void processFarText();
    void processImagePanText();
    void processImageZoomText();
    void perspectiveToggled(bool val);
    void viewButtonClicked(int index);

    void processEyeAngleText();
    void eyeAngleSliderChanged(int val);
    void copyViewFromCameraChecked(bool);
    void makeViewKeyframe();
    void centerChecked(bool);
    void processCenterText();

    void lockedViewChecked(bool);
    void extentTypeChanged(int);
    void resetView();
    void recenterView();
    void undoView();
    void tabSelected(const QString &);
private:
    void ParseViewCommands(const char *str);
    void Pan(double panx, double pany);
    void RotateAxis(int axis, double angle);
    void Zoom(double zoom);
    void Viewport(const double *viewport);
    void Window(const double *window);
    void UpdateEyeAngleSliderFromAtts(void);

    ViewCurveAttributes *viewCurve;
    View2DAttributes    *view2d;
    View3DAttributes    *view3d;
    WindowInformation   *windowInfo;
    int                 activeTab;
    bool                activeTabSetBySlot;

    // Curve widgets
    QVBox       *pageCurve;
    QGroupBox   *viewCurveGroup;
    QLineEdit   *viewportCurveLineEdit;
    QLineEdit   *domainLineEdit;
    QLineEdit   *rangeLineEdit;

    // 2d widgets
    QVBox        *page2D;
    QGroupBox    *view2DGroup;
    QLineEdit    *viewportLineEdit;
    QLineEdit    *windowLineEdit;
    QLabel       *fullFrameLabel;
    QButtonGroup *fullFrameActivationMode;
    QRadioButton *fullFrameAuto;
    QRadioButton *fullFrameOn;
    QRadioButton *fullFrameOff;

    // 3d widgets
    QVBox       *page3D;
    QGroupBox   *view3DGroup;
    QLineEdit   *normalLineEdit;
    QLineEdit   *focusLineEdit;
    QLineEdit   *upvectorLineEdit;
    QLineEdit   *viewAngleLineEdit;
    QLineEdit   *parallelScaleLineEdit;
    QLineEdit   *nearLineEdit;
    QLineEdit   *farLineEdit;
    QLineEdit   *imagePanLineEdit;
    QLineEdit   *imageZoomLineEdit;
    QLineEdit   *eyeAngleLineEdit;
    QSlider     *eyeAngleSlider;
    QCheckBox   *perspectiveToggle;
    QComboBox   *alignComboBox;

    // Global and advanced option widgets
    QTabWidget  *tabs;
    QLineEdit   *commandLineEdit;
    QComboBox   *extentComboBox;
    QCheckBox   *lockedViewToggle;
    QCheckBox   *copyViewFromCameraToggle;
    QPushButton *makeViewKeyframeButton;
    QVBox       *pageAdvanced;
    QCheckBox   *centerToggle;
    QLineEdit   *centerLineEdit;
};

#endif
