#ifndef QVIS_GUI_APPLICATION_H
#define QVIS_GUI_APPLICATION_H
#include <gui_exports.h>
#include <vector>
#include <string>
#include <qobject.h>
#include <visit-config.h>
#include <ViewerProxy.h>
#include <ConfigManager.h>
#include <GUIBase.h>
#include <MessageAttributes.h>
#include <StatusSubject.h>
#include <QualifiedFilename.h>

// Long list of forward declarations.
class DataNode;
class ObserverToCallback;
class QApplication;
class QPrinter;
class QSocketNotifier;
class QvisAnimationWindow;
class QvisAnnotationWindow;
class QvisAppearanceWindow;
class QvisArbitrarySliceWindow;
class QvisColorTableWindow;
class QvisExpressionsWindow;
class QvisCommandLineWindow;
class QvisEngineWindow;
class QvisFileInformationWindow;
class QvisFileSelectionWindow;
class QvisGlobalLineoutWindow;
class QvisHelpWindow;
class QvisHostProfileWindow;
class QvisKeyframeWindow;
class QvisLightingWindow;
class QvisMainWindow;
class QvisMaterialWindow;
class QvisMessageWindow;
class QvisOnionPeelWindow;
class QvisOutputWindow;
class QvisPickWindow;
class QvisPluginWindow;
class QvisPreferencesWindow;
class QvisQueryWindow;
class QvisRenderingWindow;
class QvisSaveWindow;
class QvisSubsetWindow;
class QvisViewWindow;
class QvisWindowBase;
class SplashScreen;

// Create a type for a vector of postable windows.
typedef std::vector<QvisWindowBase *> WindowBaseVector;

// ****************************************************************************
// Class: QvisGUIApplication
//
// Purpose:
//    This is the main application object for the VisIt graphical user
//    interface. It contains all windows, proxies and important state
//    information.
//
// Notes:      The viewer proxy pointer in GUIBase is allocated in
//             this class's constructor.
//
// Programmer: Brad Whitlock
// Creation:   Thu Aug 31 14:12:07 PST 2000
//
// Modifications:
//    Brad Whitlock, Fri Feb 9 17:36:37 PST 2001
//    I added a save window.
//
//    Brad Whitlock, Fri Feb 16 14:47:49 PST 2001
//    I added a contour plot window.
//
//    Kathleen Bonnell, Tue Mar  6 16:03:13 PST 2001 
//    I added a surface plot window.
//
//    Eric Brugger, Wed Mar  7 16:20:47 PST 2001
//    I changed the types of the plot attribute windows to the more
//    generic QvisPostableWindowObserver type.
//
//    Brad Whitlock, Fri Mar 23 11:00:47 PDT 2001
//    Added the vector plot window.
//
//    Brad Whitlock, Fri Mar 23 16:16:26 PST 2001
//    Modified things so plot and operator windows are treated as plugins.
//
//    Brad Whitlock, Thu Mar 29 14:48:53 PST 2001
//    Added the localOnly member.
//
//    Brad Whitlock, Thu Apr 19 10:38:05 PDT 2001
//    Added slots to handle window iconification.
//
//    Brad Whitlock, Tue Apr 24 16:15:13 PST 2001
//    Added a flag called viewerIsAlive.
//
//    Brad Whitlock, Tue May 1 16:05:51 PST 2001
//    Added the engine window.
//
//    Brad Whitlock, Mon Jun 11 13:59:27 PST 2001
//    Added the color table window.
//
//    Brad Whitlock, Sun Jun 17 20:58:38 PST 2001
//    Added the annotation window.
//
//    Brad Whitlock, Thu May 24 14:54:44 PST 2001
//    Added the subset window.
//
//    Jeremy Meredith, Fri Jul 20 11:26:44 PDT 2001
//    Added CalculateShift and shiftX/Y.
//
//    Jeremy Meredith, Fri Jul 20 13:56:30 PDT 2001
//    Added CalculateScreen and screenX/Y/W/H.
//
//    Jeremy Meredith, Mon Jul 23 16:45:34 PDT 2001
//    Added smallScreen.
//
//    Brad Whitlock, Thu Jul 26 15:42:17 PST 2001
//    Added the view window.
//
//    Sean Ahern, Fri Aug 31 23:59:54 PDT 2001
//    Added the splash screen.
//
//    Jeremy Meredith, Wed Sep  5 13:59:37 PDT 2001
//    Added plugin manager window.
//
//    Brad Whitlock, Tue Sep 4 22:48:48 PST 2001
//    Added infrastructure for gui color/style support.
//
//    Jeremy Meredith, Fri Sep 14 13:54:13 PDT 2001
//    Added preshiftX/Y.
//
//    Jeremy Meredith, Tue Sep 25 15:23:11 PDT 2001
//    Removed functions to calculate window, screen metrics.  These have
//    been encapsulated into WindowMetrics.  All platform specific code
//    has been removed from this class and its source file.
//
//    Jeremy Meredith, Fri Sep 28 13:57:44 PDT 2001
//    Added LoadPlugins method.
//
//    Brad Whitlock, Wed Oct 17 09:19:20 PDT 2001
//    Added the lighting window.
//
//    Eric Brugger, Mon Nov 19 14:58:56 PST 2001
//    Added the animation window.
//
//    Kathleen Bonnell, Wed Dec 12 12:04:58 PST 2001 
//    Added the pick window.
//
//    Jeremy Meredith, Fri Feb  1 16:15:28 PST 2002
//    Added members for QT's commandline arguments.
//
//    Brad Whitlock, Tue Feb 19 10:33:47 PDT 2002
//    Added copyright window and added mechanism for reading multiple
//    config files.
//
//    Brad Whitlock, Wed Feb 20 13:47:42 PST 2002
//    Added printer support.
//
//    Sean Ahern, Tue Apr 16 11:50:18 PDT 2002
//    Renamed MoveMainWindow to MoveAndResizeMainWindow to better suit its
//    purpose.
//
//    Jeremy Meredith, Wed May  8 15:23:03 PDT 2002
//    Added keyframe window.
//
//    Brad Whitlock, Thu May 9 17:02:26 PST 2002
//    Got rid of fileServer since it is now a static member of the base class.
//
//    Brad Whitlock, Thu Jul 11 16:49:40 PST 2002
//    Added the help window.
//
//    Brad Whitlock, Tue Aug 20 13:53:31 PST 2002
//    Added the file information window.
//
//    Brad Whitlock, Wed Sep 18 17:29:21 PST 2002
//    Added the rendering window.
//
//    Brad Whitlock, Fri Sep 6 16:04:40 PST 2002
//    Added the query window.
//
//    Jeremy Meredith, Thu Oct 24 16:15:11 PDT 2002
//    Added material options window.
//
//    Brad Whitlock, Thu Dec 26 16:47:18 PST 2002
//    I changed the interface for the StartMDServer method.
//
//    Kathleen Bonnell, Wed Feb 19 13:13:24 PST 2003
//    Added Global Lineout Window. 
//
//    Eric Brugger, Thu Mar 13 11:34:20 PST 2003
//    Added preferences window.
//
//    Brad Whitlock, Mon May 5 14:02:12 PST 2003
//    I changed the interface to StartMDServer.
//
//    Brad Whitlock, Mon Jun 23 11:42:54 PDT 2003
//    I added RefreshFileList.
//
//    Brad Whitlock, Mon Jul 14 11:50:34 PDT 2003
//    I added RestoreSession and SaveSession.
//
// ****************************************************************************

class GUI_API QvisGUIApplication : public QObject, public ConfigManager, public GUIBase
{
    Q_OBJECT
public:
    QvisGUIApplication(int &argc, char **argv);
    ~QvisGUIApplication();
    int Exec();

private:
    void AddViewerArguments(int argc, char **argv);
    void AddViewerSpaceArguments();
    void SplashScreenProgress(const char *, int);
    void CalculateViewerArea(int orientation, int &x, int &y,
                             int &width, int &height);
    void CreateMainWindow();
    bool CreateWindows(int startPercent, int endPercent);
    void CreatePluginWindows();
    void LaunchViewer();
    void InitializeFileServer(DataNode *);
    void LoadFile();
    void MoveAndResizeMainWindow(int orientation);
    void ProcessArguments(int &argc, char **argv);
    virtual DataNode *ReadConfigFile(const char *filename);
    void ProcessConfigSettings(DataNode *settings, bool systemConfig);
    void ProcessWindowConfigSettings(DataNode *settings);
    void SetOrientation(int orientation);
    virtual void WriteConfigFile(const char *filename);

    void ReadPluginWindowConfigs(DataNode *parentNode, const char *configVersion);
    void WritePluginWindowConfigs(DataNode *parentNode);
    void Synchronize(int tag);
    void HandleSynchronize(int val);

    // Internal callbacks
    static void StartMDServer(const std::string &hostName,
                              const stringVector &args, void *data);
    static void UpdatePrinterAttributes(Subject *subj, void *data);
    static void SyncCallback(Subject *s, void *data);
private slots:
    void HeavyInitialization();
    void ReadFromViewer(int sock_id);
    void SaveSettings();
    void ActivatePlotWindow(int index);
    void ActivateOperatorWindow(int index);
    void IconifyWindows();
    void DeIconifyWindows();
    void AboutVisIt();
    void CustomizeAppearance(bool notify);
    void LoadPlugins();
    void FinalInitialization();

    void SaveWindow();
    void SetPrinterOptions();
    void PrintWindow();
    void RefreshFileList();
    void RestoreSession();
    void SaveSession();
private:
    int                          completeInit;
    bool                         viewerIsAlive;

    MessageAttributes            message;

    // The Application
    QApplication                 *mainApp;

    // A socket notifier to tell us when the viewer proxy has input.
    QSocketNotifier              *fromViewer;
    bool                          allowSocketRead;

    QPrinter                     *printer;
    ObserverToCallback           *printerObserver;

    SplashScreen                 *splash;
    bool                          showSplash;

    // Commandline arguments for QT
    int    qt_argc;
    char **qt_argv;

    // Windows.
    QvisMainWindow               *mainWin;
    QvisColorTableWindow         *colorTableWin;
    QvisExpressionsWindow        *exprWin;
    QvisCommandLineWindow        *cliWin;
    QvisFileSelectionWindow      *fileWin;
    QvisFileInformationWindow    *fileInformationWin;
    QvisMessageWindow            *messageWin;
    QvisOutputWindow             *outputWin;
    QvisHostProfileWindow        *hostWin;
    QvisSaveWindow               *saveWin;
    QvisEngineWindow             *engineWin;
    QvisAnnotationWindow         *annotationWin;
    QvisSubsetWindow             *subsetWin;
    QvisViewWindow               *viewWin;
    QvisPluginWindow             *pluginWin;
    QvisAppearanceWindow         *appearanceWin;
    QvisLightingWindow           *lightingWin;
    QvisMaterialWindow           *materialWin;
    QvisPickWindow               *pickWin;
    QvisAnimationWindow          *animationWin;
    QvisKeyframeWindow           *keyframeWin;
    QvisHelpWindow               *helpWin;
    QvisQueryWindow              *queryWin;
    QvisPreferencesWindow        *preferencesWin;
    QvisRenderingWindow          *renderingWin;
    QvisGlobalLineoutWindow      *globalLineoutWin;

    // Contains pointers to all of the plot windows.
    WindowBaseVector             plotWindows;

    // Contains pointers to all of the operator windows.
    WindowBaseVector             operatorWindows;

    // Contains pointers to all of the other windows.
    WindowBaseVector             otherWindows;

    // Options that can be set via the command line
    bool                         localOnly;
    bool                         readConfig;
    // These hold the border sizes of the main window.
    int                          borders[4];
    // Other sizing parameters for windows
    int                          shiftX;
    int                          shiftY;
    int                          preshiftX;
    int                          preshiftY;
    int                          screenX;
    int                          screenY;
    int                          screenW;
    int                          screenH;
    // Some appearance attributes
    QString                      foregroundColor;
    QString                      backgroundColor;
    QString                      applicationStyle;



    // File to load on startup.
    QualifiedFilename            loadFile;

    int                          sessionCount;

    // Synchronization attributes
    ObserverToCallback           *syncObserver;
    int                          initStage;
    int                          heavyInitStage;
    int                          windowInitStage;
    int                          windowTimeId;

    // Config file storage
    DataNode                     *systemSettings;
    DataNode                     *localSettings;
};

#endif
