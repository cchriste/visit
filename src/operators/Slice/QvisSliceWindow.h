#ifndef QVIS_SLICE_WINDOW_H
#define QVIS_SLICE_WINDOW_H
#include <QvisOperatorWindow.h>

#include <QNarrowLineEdit.h>

// Forward declarations.
class QCheckBox;
class QLineEdit;
class QLabel;
class QButtonGroup;
class SliceAttributes;
class QSlider;
class QGridLayout;
class QHBoxLayout;

// ****************************************************************************
// Class: QvisSliceWindow
//
// Purpose:
//   This class is a postable window that watches slice operator
//   attributes and always represents their current state.
//
// Notes:      
//
// Programmer: Brad Whitlock
// Creation:   Fri Aug 4 10:13:05 PDT 2000
//
// Modifications:
//   Jeremy Meredith, Wed Sep  5 15:50:26 PDT 2001
//   Renamed to match plugin name.
//
//   Brad Whitlock, Tue Oct 9 17:35:43 PST 2001
//   Added an "interactive" toggle.
//
//   Brad Whitlock, Tue Feb 26 14:18:00 PST 2002
//   Added a new slot function.
//
//   Brad Whitlock, Fri Apr 12 12:58:39 PDT 2002
//   Made it inherit from QvisOperatorWindow.
//
//   Brad Whitlock, Mon Nov 18 10:49:05 PDT 2002
//   Added upAxisLabel.
//
//   Jeremy Meredith, Mon May  5 14:47:53 PDT 2003
//   Reworked a lot of stuff for changes to how origin is handled.
//
// ****************************************************************************

class QvisSliceWindow : public QvisOperatorWindow
{
    Q_OBJECT
public:
    QvisSliceWindow(const int type,
                    SliceAttributes *subj,
                    const char *caption = 0,
                    const char *shortName = 0,
                    QvisNotepadArea *notepad = 0);
    virtual ~QvisSliceWindow();
protected:
    virtual void CreateWindowContents();
    void UpdateWindow(bool doAll);
    virtual void GetCurrentValues(int which_widget);
    void UpdateOriginArea();
private slots:
    void processNormalText();
    void processUpAxisText();
    void projectToggled(bool);
    void interactiveToggled(bool);
    void flipNormalToggled(bool);
    void normalTypeChanged(int index);
    void originTypeChanged(int index);

    void processOriginPointText();
    void processOriginInterceptText();
    void processOriginPercentText();
    void processOriginZoneText();
    void processOriginNodeText();

    void originPercentSliderPressed();
    void originPercentSliderReleased();
    void originPercentSliderChanged(int pct);
private:
    QCheckBox        *interactiveToggle;

    QButtonGroup     *normalTypeGroup;
    QLineEdit        *normalLineEdit;
    QCheckBox        *flipNormalToggle;

    QButtonGroup     *originTypeGroup;
    QHBoxLayout      *originPointLayout;
    QLineEdit        *originPointLineEdit;
    QHBoxLayout      *originInterceptLayout;
    QLineEdit        *originInterceptLineEdit;
    QHBoxLayout      *originPercentLayout;
    QNarrowLineEdit  *originPercentLineEdit;
    QHBoxLayout      *originZoneLayout;
    QLineEdit        *originZoneLineEdit;
    QHBoxLayout      *originNodeLayout;
    QLineEdit        *originNodeLineEdit;

    QSlider          *originPercentSlider;
    bool              sliderDragging;

    QLabel           *originPointLabel;
    QLabel           *originInterceptLabel;
    QLabel           *originPercentLabel;
    QLabel           *originZoneLabel;
    QLabel           *originNodeLabel;

    QCheckBox        *projectToggle;
    QLineEdit        *upAxisLineEdit;
    QLabel           *upAxisLabel;

    SliceAttributes  *sliceAtts;
};
#endif
