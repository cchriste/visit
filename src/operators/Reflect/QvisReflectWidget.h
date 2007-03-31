#ifndef QVIS_REFLECT_WIDGET_H
#define QVIS_REFLECT_WIDGET_H
#include <qwidget.h>
#include <mini3D.h>

class QTimer;

// ****************************************************************************
// Class: QvisReflectWidget
//
// Purpose:
//   This is a 3D widget that allows us to set the octants for the reflect
//   operator.
//
// Notes:      
//
// Programmer: Brad Whitlock
// Creation:   Wed Mar 5 15:45:14 PST 2003
//
// Modifications:
//   Brad Whitlock, Mon Jun 23 16:49:01 PST 2003
//   I added a 2d interaction mode.
//
// ****************************************************************************

class QvisReflectWidget : public QWidget
{
    Q_OBJECT
public:
    QvisReflectWidget(QWidget *parent, const char *name=0);
    virtual ~QvisReflectWidget();
    virtual QSize sizeHint() const;
    virtual QSizePolicy sizePolicy() const;

    void setMode2D(bool val);
    bool getMode2D() const;
signals:
    void octantChanged(int);
    void valueChanged(bool *octants);
public slots:
    void setValues(bool *octants);
    void setOriginalOctant(int octant);
protected slots:
    void handleTimer();
protected:
    void drawOnOffActors(int n, float scale);
    void deleteBackingPixmap();
    void redrawScene(QPainter *painter);
    void redrawScene2D(QPainter *painter);
    void redrawScene3D(QPainter *painter);
    void setupAndDraw(QPainter *p);
    void setupCamera();
    void createSharedElements();
    void initializeAxes();
    void initializeAxes2D();
    void initializeArrow();
    void initializeSphere(m3d_complex_element &, int nx, int ny, float rad,
                          float r, float g, float b);
    void initializeCube(m3d_complex_element &, int nx, int ny,
                        float s, float r, float g, float b);
    void ScaleTranslateFromOriginToOctant(int octant, float s);

    virtual void mouseReleaseEvent(QMouseEvent *e);
    virtual void paintEvent(QPaintEvent *e);
    virtual void resizeEvent(QResizeEvent *e);

    bool          mode2D;

    QPixmap      *pixmap;
    m3d_renderer  renderer;
    bool          rendererCreated;

    int           originOctant;
    bool          octantOn[8];

    QTimer       *timer;
    int           activeCamera;
    float         cameraInterpolant;
    bool          switchingCameras;

    static bool                sharedElementsCreated;
    static m3d_complex_element axes;
    static m3d_complex_element axes2D;
    static m3d_complex_element onCube;
    static m3d_complex_element offCube;
    static m3d_complex_element onSphere;
    static m3d_complex_element offSphere;
    static m3d_complex_element arrow;
};


#endif
