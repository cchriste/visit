#ifndef QVISINTERACTORWINDOW_H
#define QVISINTERACTORWINDOW_H

#include <AttributeSubject.h>
#include <QvisPostableWindowObserver.h>

class InteractorAttributes;
class QLabel;
class QCheckBox;
class QLineEdit;
class QSpinBox;
class QVBox;
class QButtonGroup;
class QvisColorTableButton;
class QvisOpacitySlider;
class QvisColorButton;
class QvisLineStyleWidget;
class QvisLineWidthWidget;

// ****************************************************************************
// Class: QvisInteractorWindow
//
// Purpose: 
//   Defines QvisInteractorWindow class.
//
// Notes:      This class was automatically generated!

// Programmer: xml2window
// Creation:   Mon Aug 16 15:29:28 PST 2004
//
// Modifications:
//   Eric Brugger, Thu Nov 11 11:47:28 PST 2004
//   I added the navigation mode toggle buttons.
//   
//   Eric Brugger, Thu Nov 18 13:01:35 PST 2004
//   I added the fill viewport on zoom toggle button.
//   
// ****************************************************************************

class QvisInteractorWindow : public QvisPostableWindowObserver
{
    Q_OBJECT
  public:
    QvisInteractorWindow(
                         InteractorAttributes *subj,
                         const char *caption = 0,
                         const char *shortName = 0,
                         QvisNotepadArea *notepad = 0);
    virtual ~QvisInteractorWindow();
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
    void showGuidelinesChanged(bool val);
    void clampSquareChanged(bool val);
    void fillViewportOnZoomChanged(bool val);
    void navigationModeChanged(int val);
  private:
    QCheckBox *showGuidelines;
    QCheckBox *clampSquare;
    QCheckBox *fillViewportOnZoom;
    QButtonGroup *navigationMode;

    InteractorAttributes *atts;
};



#endif
