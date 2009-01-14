#ifndef QVIS_CMAP2_WIDGET_H
#define QVIS_CMAP2_WIDGET_H

#include <QvisCMap2Display.h> // Must come before other Qt includes because of SLIVR's
                              // "emit" method conflicing with Qt's emit macro.
#include <QWidget>

class QLineEdit;
class QListWidget;
class QPushButton;
class QvisColorButton;
class QvisOpacitySlider;

// ****************************************************************************
// Class: QvisCMap2Widget
//
// Purpose:
//   Simple 2D transfer function editing widget based on SLIVR 2D widgets.
//
// Notes:
//
// Programmer: Brad Whitlock
// Creation:   Tue Mar 25 13:53:36 PDT 2008
//
// Modifications:
//   Josh Stratton, August 2008
//   Added additional SLIVR widget types.
//
//   Brad Whitlock, Tue Sep 30 09:48:58 PDT 2008
//   Qt 4.
//
//   Brad Whitlock, Mon Jan 12 14:35:56 PST 2009
//   I added selectWidget.
//
// ****************************************************************************

class QvisCMap2Widget : public QWidget
{
    Q_OBJECT
public:
    QvisCMap2Widget(QWidget *parent);
    virtual ~QvisCMap2Widget();

    void setHistogramTexture(const unsigned char *data, int size);

    WidgetID addTriangleWidget(const QString &wName,
                          float base,  // X-coordinate of bottom point
                          float top_x, // Offset in X from top left point
                          float top_y, // Y-coordinate of top left point
                          float w,     // Width of top of triangle.
                          float bottom // Distance along triangle median
                          );

    WidgetID addRectangleWidget(const QString &wName,
                           float left_x,
                           float left_y,
                           float w,
                           float h,
                           float offset
                           );

    WidgetID addEllipsoidWidget(const QString &wName,
                                float x,
                                float y,
                                float a, 
                                float b,
                                float rot
                                );

    WidgetID addParaboloidWidget(const QString &wName,
                                 float top_x, float top_y,
                                 float bottom_x, float bottom_y,
                                 float left_x, float left_y,
                                 float right_x, float right_y
                                 );

    int numWidgets() const;
    WidgetID getID(int index) const;
    void removeWidget(WidgetID id);
    void clear();

    void   setDefaultColor(const QColor &c);
    QColor getDefaultColor() const;

    void   setDefaultAlpha(float f);
    float  getDefaultAlpha() const;

    void   setColor(WidgetID id, const QColor &c);
    QColor getColor(WidgetID id) const;

    void   setAlpha(WidgetID id, float a);
    float  getAlpha(WidgetID id) const;

    void    setName(WidgetID id, const QString &name);
    QString getName(WidgetID id) const;

    QString getString(WidgetID id) const;

    WidgetID getSelectedWidget() const;
    void     selectWidget(WidgetID id);
signals:
    void widgetListChanged();
    void selectedWidget(WidgetID id);
    void widgetChanged(WidgetID id);
private slots:
    void updateList();
    void addNewTriangle();
    void addNewRectangle();
    void addNewEllipsoid();
    void addNewParaboloid();
    void deleteWidget();
    void showProperties(int);
    void selectLastItem();
    void updateWidget(WidgetID);
    void setWidgetColor(const QColor &);
    void setWidgetOpacity(int val);
    void setSizeLoc();
private:
    QString newName(const QString &base);

    int                nameIndex;
    QvisCMap2Display  *cmap2;
    QListWidget       *names;
    QPushButton       *deleteButton;
    QLineEdit         *wName;
    QLineEdit         *sizeLoc;
    QvisColorButton   *color;
    QvisOpacitySlider *opacity;
};

#endif
