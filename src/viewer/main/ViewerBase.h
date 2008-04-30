#ifndef VIEWER_BASE_H
#define VIEWER_BASE_H
#include <qobject.h>
#include <viewer_exports.h>

class ViewerMethods;
class ViewerState;

// ****************************************************************************
// Class: ViewerBase
//
// Purpose:
//   Base class for viewer classes
//
// Notes:      
//
// Programmer: Brad Whitlock
// Creation:   Mon Feb 12 09:25:39 PDT 2007
//
// Modifications:   
//    Cyrus Harrison, Thu Feb 21 15:04:14 PST 2008
//    Added control to suppress messages. 
//
//    Brad Whitlock, Tue Apr 29 11:11:04 PDT 2008
//    Converted most messaging to QString.
//
// ****************************************************************************

class VIEWER_API ViewerBase : public QObject
{
public:
    ViewerBase(QObject *parent = 0, const char *name = 0);
    virtual ~ViewerBase();

    // Methods to get pointers to the client state and the methods.
    static ViewerState   *GetViewerState();
    static ViewerMethods *GetViewerMethods();

    // Methods to send messages
    static void Error(const QString &message, bool = true);
    static void Warning(const QString &message, bool = true);
    static void Message(const QString &message, bool = true);
    static void ErrorClear();

    // Methods to send status
    static void Status(const QString &message);
    static void Status(const QString &message, int milliseconds);
    static void Status(const char *sender, const QString &message);
    static void Status(const char *sender, int percent, int curStage,
                const char *curStageName, int maxStage);
    static void ClearStatus(const char *sender = 0);
    
    // Message Suppression Control
    static void EnableMessageSuppression()  { suppressMessages = true;}
    static void DisableMessageSuppression() { suppressMessages = false;}
    static bool SuppressMessages()          { return suppressMessages;}
    
private:
    static ViewerState   *base_viewerState;
    static ViewerMethods *base_viewerMethods;
    static bool           suppressMessages;
};

#endif
