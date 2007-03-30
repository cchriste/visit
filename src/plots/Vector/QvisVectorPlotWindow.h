#ifndef QVIS_VECTOR_WINDOW_H
#define QVIS_VECTOR_WINDOW_H
#include <QvisPostableWindowObserver.h>

// Forward declarations
class QButtonGroup;
class QCheckBox;
class QLabel;
class QLineEdit;
class QvisColorButton;
class QvisLineStyleWidget;
class QvisLineWidthWidget;
class VectorAttributes;
class QvisOpacitySlider;
class QvisColorTableButton;

// ****************************************************************************
// Class: QvisVectorPlotWindow
//
// Purpose:
//   This class is a postable window that watches vector plot attributes and
//   always represents their current state.
//
// Notes:      
//
// Programmer: Hank Childs & Brad Whitlock
// Creation:   Thu Mar 22 23:40:52 PST 2001
//
// Modifications:
//   Brad Whitlock, Sat Jun 16 18:21:34 PST 2001
//   I added color table stuff.
//   
// ****************************************************************************

class QvisVectorPlotWindow : public QvisPostableWindowObserver
{
    Q_OBJECT
public:
    QvisVectorPlotWindow(const int type, VectorAttributes *_vecAtts,
                         const char *caption = 0,
                         const char *shortName = 0,
                         QvisNotepadArea *notepad = 0);
    virtual ~QvisVectorPlotWindow();
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
    void lineStyleChanged(int newStyle);
    void lineWidthChanged(int newWidth);
    void vectorColorChanged(const QColor &color);
    void processScaleText();
    void processHeadSizeText();
    void reduceMethodChanged(int index);
    void processNVectorsText();
    void processStrideText();
    void legendToggled();
    void drawHeadToggled();
    void colorModeChanged(int);
    void colorTableClicked(bool useDefault, const QString &ctName);
private:
    int                  plotType;
    VectorAttributes     *vectorAtts;

    QvisLineStyleWidget  *lineStyle;
    QvisLineWidthWidget  *lineWidth;
    QvisColorButton      *vectorColor;
    QLineEdit            *scaleLineEdit;
    QLineEdit            *headSizeLineEdit;
    QButtonGroup         *reduceButtonGroup;
    QLineEdit            *nVectorsLineEdit;
    QLineEdit            *strideLineEdit;
    QCheckBox            *legendToggle;
    QCheckBox            *drawHeadToggle;
    QButtonGroup         *colorButtonGroup; 
    QvisColorTableButton *colorTableButton;
};

#endif
