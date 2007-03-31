// ************************************************************************* //
//                            ViewerWindowManager.h                          //
// ************************************************************************* //

#ifndef VIEWER_WINDOW_MANAGER_H
#define VIEWER_WINDOW_MANAGER_H
#include <viewer_exports.h>
#include <qobject.h>
#include <VisWindowTypes.h>
#include <avtTypes.h>
#include <avtImage.h>
#include <avtDataset.h>
#include <avtView2D.h>
#include <avtView3D.h>
#include <vectortypes.h>

#define VIEWER_WINDOW_MANAGER_VSTACK 10

class AnimationAttributes;
class AnnotationAttributes;
class GlobalAttributes;
class KeyframeAttributes;
class LightList;
class Line;
class PrinterAttributes;
class RenderingAttributes;
class SaveWindowAttributes;
class VisWindow;
class ViewerWindow;
class ViewerAnimation;
class ViewAttributes;
class ViewerWindowManagerAttributes;
class WindowInformation;
class avtFileWriter;
class avtToolInterface;
class QTimer;

typedef struct {
   int       x;
   int       y;
   int       width;
   int       height;
} WindowLimits;

// ****************************************************************************
//  Class: ViewerWindowManager
//
//  Purpose:
//    ViewerWindowManager is a class that manages viewer windows.  It has
//    methods for adding and deleting windows, specifying window layouts,
//    and setting the active window.  It is modeled after the singleton
//    class which means that applications can't directly instantiate the
//    class but can only get a reference to the one instance through the
//    Instance method.
//
//  Note:  The class currently supports a maximum of 16 windows.
//
//  Programmer: Eric Brugger
//  Creation:   September 7, 2000
//
//  Modifications:
//    Brad Whitlock, Tue Nov 7 10:08:37 PDT 2000
//    Changed to use ViewerWindow.
//
//    Hank Childs, Fri Nov 10 11:20:19 PST 2000
//    Added TogglePerspective.
//
//    Brad Whitlock, Mon Nov 27 14:28:41 PST 2000
//    I made it inherit from QObject so I could add some signals that are
//    emitted when creating and deleting windows.
//
//    Eric Brugger, Fri Feb 23 12:33:00 PST 2001
//    I added RecenterView and ToggleAutoCenterMode.
//
//    Brad Whitlock, Thu Apr 19 11:10:12 PDT 2001
//    I added methods to handle window iconification.
//
//    Brad Whitlock, Thu Jun 14 16:22:21 PST 2001
//    Added a method to update the color table.
//
//    Jeremy Meredith, Fri Jul 20 11:23:42 PDT 2001
//    Added SetShift and shiftX/Y for correcting some window managers.
//
//    Brad Whitlock, Thu Jul 19 16:13:24 PST 2001
//    Added methods to set the view.
//
//    Eric Brugger, Tue Aug 21 10:16:06 PDT 2001
//    Added the ViewCallback method.  Removed an argument from UpdateVeiwAtts.
//
//    Brad Whitlock, Wed Aug 22 11:25:07 PDT 2001
//    Added a second DeleteWindow method.
//
//    Brad Whitlock, Mon, Aug 26 23:44:23 PST 2001
//    I added a new slot function and converted some functions to slots. I
//    also moved the annotation attributes into this class.
//
//    Jeremy Meredith, Fri Sep 14 13:49:38 PDT 2001
//    Added SetPreshift and preshiftX/Y for further correction of some
//    window managers.
//
//    Kathleen Bonnell, Wed Sep 26 11:49:18 PDT 2001 
//    Added SetPlotColors to send foreground/background color information
//    to the plots.
//
//    Brad Whitlock, Mon Oct 1 16:11:38 PST 2001
//    Added methods to handle tools.
//
//    Brad Whitlock, Fri Sep 14 14:19:57 PST 2001
//    Added methods to deal with lighting.
//
//    Eric Brugger, Fri Oct 26 14:31:59 PDT 2001
//    I added a timer for playing animations along with the method
//    UpdateAnimationTimer to turn it on and off appropriately.  I also
//    added an animation argument to UpdateAnimationState so that it would
//    only update the state if the animation was the active animation.
//
//    Eric Brugger, Fri Nov  2 12:37:44 PST 2001
//    I added a data member to keep track of whether the window has been
//    referenced before.
//
//    Brad Whitlock, Fri Nov 2 10:28:46 PDT 2001
//    Added a method to get the viewer window manager attributes.
//
//    Eric Brugger, Wed Nov 21 12:11:35 PST 2001
//    I added methods to deal with animation attributes.
//
//    Kathleen Bonnell, Tue Nov 27 16:03:00 PST 2001 
//    Added pick attributes and related methods.
//    Added method 'ClearPickPoints'. 
//
//    Brad Whitlock, Mon Jan 28 16:27:14 PST 2002
//    Added support for a stack of views.
//
//    Brad Whitlock, Fri Feb 8 10:29:32 PDT 2002
//    Added a method to query all plots in all windows to see if a file
//    is being used for the plot.
//
//    Brad Whitlock, Mon Feb 11 17:26:50 PST 2002
//    Added a method to set the active plots for a window.
//
//    Brad Whitlock, Wed Feb 20 14:34:42 PST 2002
//    I added printing support.
//
//    Sean Ahern, Tue Apr 16 12:36:24 PDT 2002
//    Added the ability to show/hide all windows and have deferred creation.
//
//    Kathleen Bonnell, Thu May  9 09:57:56 PDT 2002  
//    Added methods in support of Lineout: ClearRefLines, AddToLineoutWindow,
//    ResetLineoutDesignation.
//
//    Hank Childs, Thu May 23 17:33:41 PDT 2002
//    Added a more generic file writer that can take care of files and images.
//
//    Hank Childs, Wed May 29 08:42:29 PDT 2002
//    Added ToggleSpinMode.
//
//    Kathleen Bonnell, Wed May 29 11:24:14 PDT 2002  
//    Changed arguments to method AddtoLineoutWindow. 
//
//    Kathleen Bonnell, Mon Jun 10 17:22:57 PDT 2002 
//    Replace 'AddToLineoutWindow' with 'GetLineoutWindow'. 
//
//    Brad Whitlock, Thu Jun 27 16:37:21 PST 2002
//    Added methods to copy attributes from one window to another window.
//
//    Hank Childs, Mon Jul 15 09:58:27 PDT 2002
//    Added SetViewExtentsType.
//
//    Brad Whitlock, Tue Jul 23 16:15:07 PST 2002
//    I added a method to disable all tools.
//
//    Brad Whitlock, Thu Sep 19 13:27:06 PST 2002
//    I added a new state object.
//
//    Brad Whitlock, Fri Sep 6 15:17:26 PST 2002
//    I added methods to do pick and lineout.
//
//    Brad Whitlock, Tue Oct 15 16:36:41 PST 2002
//    I added methods to clone a window and copy plots to a window.
//
//    Brad Whitlock, Wed Oct 30 14:53:55 PST 2002
//    I added methods to tell the client and the viewer windows when an
//    engine is busy executing.
//
//    Brad Whitlock, Mon Nov 11 11:55:57 PDT 2002
//    I added methods to toggle locktime and locktools.
//
//    Kathleen Bonnell, Fri Nov 15 09:07:36 PST 2002 
//    Removed PickAttributes and related Set/Get methods. Now handled
//    by ViewerQueryManager. 
//
//    Brad Whitlock, Tue Nov 19 14:44:36 PST 2002
//    I added more default arguments to UpdateViewAtts.
//
//    Eric Brugger, Mon Nov 25 09:58:47 PST 2002
//    I added GetKeyframeClientAtts and SetKeyframeAttsFromClient.
//
//    Kathleen Bonnell, Mon Dec 23 13:14:18 PST 2002 
//    Remove Lineout.
//
//    Eric Brugger, Mon Jan  6 12:27:06 PST 2003
//    I added ClearViewKeyframes, DeleteViewKeyframe, SetViewKeyframe and
//    ToggleCameraViewMode.
//
//    Eric Brugger, Wed Jan 29 11:08:28 PST 2003
//    I added MoveViewKeyframe.
//
//    Eric Brugger, Fri Jan 31 13:43:37 PST 2003 
//    I removed the state argument from UpdateAnimationState.
//
//    Brad Whitlock, Fri Jan 31 12:39:24 PDT 2003
//    I added a method to get the number of windows and another method to
//    update actions.
//
//    Brad Whitlock, Thu Feb 27 12:01:39 PDT 2003
//    I removed an overloaded version of UpdateAnimationState.
//
//    Brad Whitlock, Tue Apr 1 10:40:36 PDT 2003
//    I added a method to hide toolbars for all windows.
//
//    Eric Brugger, Fri Apr 11 14:01:53 PDT 2003
//    I added the method SimpleAddWindow.
//
//    Eric Brugger, Fri Apr 18 12:36:15 PDT 2003
//    I replaced ToggleAutoCenterMode with ToggleMaintainViewMode.
//
//    Kathleen Bonnell, Wed May 14 17:07:07 PDT 2003  
//    Removed 'Pick' method. 
//
//    Kathleen Bonnell, Thu May 15 10:00:02 PDT 2003 
//    Added ToggleFullFrameMode. 
//
// ****************************************************************************

class VIEWER_API ViewerWindowManager : public QObject
{
    Q_OBJECT
  public:
    virtual ~ViewerWindowManager();

    static ViewerWindowManager *Instance();

    void SetBorders(const char *windowBorders);
    void SetShift(const char *windowShift);
    void SetPreshift(const char *windowPreshift);
    void SetGeometry(const char *windowGeometry);

    void AddInitialWindows();
    void AddWindow(bool copyAtts = false);
    void CloneWindow();
    void ClearAllWindows();
    void ClearWindow(int windowIndex = -1);
    void CopyAnnotationsToWindow(int from, int to);
    void CopyLightingToWindow(int from, int to);
    void CopyViewToWindow(int from, int to);
    void CopyAnimationToWindow(int from, int to);
    void DeleteWindow();
    void DeleteWindow(ViewerWindow *win);
    void DisableRedraw(int windowIndex = -1);
    bool FileInUse(const char *host, const char *dbName) const;
    void IconifyAllWindows();
    void DeIconifyAllWindows();
    void ShowAllWindows();
    void HideAllWindows();
    bool GetWindowsHidden() const { return windowsHidden; };
    void PrintWindow(int windowIndex = -1);
    void RecenterView(int windowIndex = -1);
    void RedrawWindow(int windowIndex = -1);
    void ResetView(int windowIndex = -1);
    void SaveWindow(int windowIndex = -1);
    void SetActiveWindow(const int windowId);
    void SetAnimationAttsFromClient();
    void SetAnnotationAttsFromClient();
    void SetAnnotationAttsFromDefault();
    void SetKeyframeAttsFromClient();
    void SetFrameIndex(int frame, int windowIndex = -1);
    void NextFrame(int windowIndex = -1);
    void PrevFrame(int windowIndex = -1);
    void Stop(int windowIndex = -1);
    void Play(int windowIndex = -1);
    void ReversePlay(int windowIndex = -1);
    void SetInteractionMode(INTERACTION_MODE m, int windowIndex = -1);
    void SetLightListFromClient();
    void SetLightListFromDefault();
    void SetView2DFromClient();
    void SetView3DFromClient();
    void ClearViewKeyframes();
    void DeleteViewKeyframe(const int frame);
    void MoveViewKeyframe(int oldFrame, int newFrame);
    void SetViewKeyframe();
    void SetWindowLayout(const int windowLayout);

    void SetRenderingAttributes(int windowIndex = -1);
    void ToggleBoundingBoxMode(int windowIndex = -1);
    void TogglePerspective(int windowIndex = -1);
    void ToggleMaintainViewMode(int windowIndex = -1);
    void ToggleLockViewMode(int windowIndex = -1);
    void ToggleSpinMode(int windowIndex = -1);
    void ToggleCameraViewMode(int windowIndex = -1);
    void ToggleLockTime(int windowIndex = -1);
    void ToggleLockTools(int windowIndex = -1);
    void ToggleFullFrameMode(int windowIndex = -1);
    void SetViewExtentsType(avtExtentType, int windowIndex = -1);
    void UpdateColorTable(const char *ctName);
    void ResetLineoutDesignation(int windowIndex = -1); 
    int             *GetWindowIndices(int *nwin) const;
    int              GetNumWindows() const;
    ViewerWindow    *GetActiveWindow() const;
    ViewerWindow    *GetLineoutWindow() ;
    ViewerAnimation *GetActiveAnimation() const;
    int              GetWindowLayout() const { return layout; };
    void UpdateActions();
    void HideToolbarsForAllWindows();
    void ShowToolbarsForAllWindows();

    void UpdateAnimationTimer();
    void StopTimer();

    void UpdateAnimationState(const ViewerAnimation *, const int mode) const;
    void UpdateGlobalAtts() const;
    void UpdateViewAtts(int windowIndex = -1, bool update2d = true,
                        bool update3d = true);
    void UpdateAnimationAtts();
    void UpdateAnnotationAtts();
    void UpdateLightListAtts();
    void UpdateWindowAtts();
    void UpdateWindowInformation(int windowIndex = -1,
                                 bool reportTimes = false);
    void UpdateRenderingAtts(int windowIndex = -1);
    void UpdateAllAtts();

    void BeginEngineExecute();
    void EndEngineExecute();

    static GlobalAttributes              *GetClientAtts();
    static SaveWindowAttributes          *GetSaveWindowClientAtts();
    static ViewAttributes                *GetView2DClientAtts();
    static ViewAttributes                *GetView3DClientAtts();
    static AnimationAttributes           *GetAnimationClientAtts();
    static AnnotationAttributes          *GetAnnotationClientAtts();
    static AnnotationAttributes          *GetAnnotationDefaultAtts();
    static PrinterAttributes             *GetPrinterClientAtts();
    static ViewerWindowManagerAttributes *GetWindowAtts();
    static WindowInformation             *GetWindowInformation();
    static void                          SetClientAnnotationAttsFromDefault();
    static void                          SetDefaultAnnotationAttsFromClient();
    static KeyframeAttributes            *GetKeyframeClientAtts();
    static LightList                     *GetLightListClientAtts();
    static LightList                     *GetLightListDefaultAtts();
    static void                          SetClientLightListFromDefault();
    static void                          SetDefaultLightListFromClient();
    static void                          ViewCallback(VisWindow *);
    static RenderingAttributes           *GetRenderingAttributes();
    static void                          RenderInformationCallback(void *data);
  signals:
    void createWindow(ViewerWindow *);
    void deleteWindow(ViewerWindow *);

  public slots:
    void InvertBackgroundColor(int windowIndex = -1);
    void SetToolEnabled(int toolId, bool enabled, int windowIndex = -1);
    void DisableAllTools(int windowIndex = -1);
    void UndoView(int windowIndex = -1);

  protected:
    ViewerWindowManager();
    avtImage_p CreateSingleImage(int windowIndex, int width, int height,
                                 bool screenCapture);
    avtImage_p CreateTiledImage(int width, int height);
    avtDataset_p GetDataset(int windowIndex);
    static void ToolCallback(const avtToolInterface &);

  private:
    void InitWindowLimits();
    int  SimpleAddWindow();
    void CreateVisWindow(const int windowIndex, const int width,
                         const int height,
                         const int x, const int y);
    void SetWindowAttributes(int windowIndex, bool copyAtts);

    void SetPlotColors(const AnnotationAttributes *);

    static ViewerWindowManager           *instance;
    static ViewerWindowManagerAttributes *windowAtts;
    static WindowInformation             *windowInfo;
    static const int                     maxWindows;
    static const int                     maxLayouts;
    static const int                     validLayouts[];
    static GlobalAttributes              *clientAtts;
    static SaveWindowAttributes          *saveWindowClientAtts;
    static ViewAttributes                *view2DClientAtts;
    static ViewAttributes                *view3DClientAtts;
    static AnimationAttributes           *animationClientAtts;
    static AnnotationAttributes          *annotationClientAtts;
    static AnnotationAttributes          *annotationDefaultAtts;
    static KeyframeAttributes            *keyframeClientAtts;
    static LightList                     *lightListClientAtts;
    static LightList                     *lightListDefaultAtts;
    static PrinterAttributes             *printerAtts;
    static RenderingAttributes           *renderAtts;

    int               layout;
    int               layoutIndex;

    int               activeWindow;

    int               screenWidth;
    int               screenHeight;
    int               screenX;
    int               screenY;
    int               borderTop;
    int               borderBottom;
    int               borderLeft;
    int               borderRight;
    int               shiftX;
    int               shiftY;
    int               preshiftX;
    int               preshiftY;

    int               filesWritten;
    avtFileWriter    *fileWriter;

    int               nWindows;
    ViewerWindow      **windows;
    int               *x_locations;
    int               *y_locations;
    bool              *referenced;

    bool              windowsHidden;
    bool              windowsIconified;

    int               lineoutWindow;

    WindowLimits      **windowLimits;

    QTimer            *timer;
    int               animationTimeout;
    int               lastAnimation;

    bool              viewStacking;
    avtView2D         view2DStack[VIEWER_WINDOW_MANAGER_VSTACK];
    int               view2DStackTop;
    avtView3D         view3DStack[VIEWER_WINDOW_MANAGER_VSTACK];
    int               view3DStackTop;


private slots:
    void AnimationCallback();

};

#endif
