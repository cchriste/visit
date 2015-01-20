#ifndef SIMPLE_VIS_APP_H
#define SIMPLE_VIS_APP_H
#include <QMainWindow>

class VisItViewer;
class CommandParser;

class QLabel;
class QListWidget;
class QButtonGroup;
class QSpinBox;
class QWidget;
class vtkQtRenderWindow;

// ****************************************************************************
// Class: SimpleVisApp
//
// Purpose:
//   Simple embedded viewer application that lets you choose from Pseudocolor
//   and Contour plots of variables from an open file.
//
// Notes:      
//
// Programmer: Brad Whitlock
// Creation:   Thu Aug 21 15:14:40 PDT 2008
//
// Modifications:
//   
// ****************************************************************************

class SimpleVisApp : public QMainWindow
{
    Q_OBJECT
public:
    SimpleVisApp(VisItViewer *, int *, char ***);
    virtual ~SimpleVisApp();

public slots:
    virtual void show();
    void execFile(int openClFile=0);

private slots:
    void selectFile();
    void openFile(const QString &);
    void changeVariable(const QString &);
    void changeVariableAndUpdate(const QString &);
    void changePlotType(int);
    void setNContours(int);
    void saveWindow();
    void openGUI();
    void quitApp();
private:
    static vtkQtRenderWindow *ReturnVisWin(void *);
    void resetWindow();

    QLabel            *scalarLabel;
    QListWidget       *variables;
    QWidget           *plotTypeWidget;
    QButtonGroup      *plotType;
    QWidget           *contourWidget;
    QSpinBox          *nContours;
    vtkQtRenderWindow *viswin;
    VisItViewer       *viewer;

    CommandParser     *cmd;

    QString            clfilename;
};

#endif
