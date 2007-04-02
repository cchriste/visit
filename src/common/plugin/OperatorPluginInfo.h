// ************************************************************************* //
//                            OperatorPluginInfo.h                           //
// ************************************************************************* //

#ifndef OPERATOR_PLUGIN_INFO_H
#define OPERATOR_PLUGIN_INFO_H
#include <plugin_exports.h>
#include <stdio.h>

// Forward declarations.
class AttributeSubject;
class QvisNotepadArea;
class QvisPostableWindowObserver;
class QvisWizard;
class QWidget;
class avtPluginFilter;
class ViewerPlot;

// ****************************************************************************
//  Class: *OperatorPluginInfo
//
//  Purpose:
//    Four classes that provide all the information about an operator plugin.
//    The information is broken up into four classes relevant to the gui, the
//    viewer, the engine, or all three.
//
//  Programmer: Jeremy Meredith
//  Creation:   July 25, 2001
//
//  Modifications:
//    Jeremy Meredith, Fri Sep 28 13:31:21 PDT 2001
//    Added Common info.  Made the specific ones inherit from the generic ones.
//
//    Brad Whitlock, Thu Jul 26 12:52:12 PDT 2001
//    Added scripting plugins.
//
//    Brad Whitlock, Mon May 6 14:18:15 PST 2002
//    I added the SetDefaults method for scripting plugins.
//
//    Brad Whitlock, Wed Mar 12 14:14:57 PST 2003
//    I added the XPMIconData method for GUI and Viewer plugins.
//
//    Brad Whitlock, Thu Apr 10 13:14:52 PST 2003
//    I added Removeable, Moveable, and AllowsSubsequentOperators.
//
//    Kathleen Bonnell, Thu Sep 11 10:18:54 PDT 2003
//    I added a bool argument to InitializeOperatorAtts. 
//
//    Jeremy Meredith, Wed Nov  5 10:28:29 PST 2003
//    Added ability to disable plugins by default.
//
//    Brad Whitlock, Tue Dec 14 09:44:04 PDT 2004
//    Added support for operator wizards.
//
//    Hank Childs, Tue Mar 22 16:06:15 PST 2005
//    Made destructor virtual.
//
//    Brad Whitlock, Thu Jan 5 14:40:42 PST 2006
//    I removed explicit support for logging from the scripting plugins.
//
// ****************************************************************************

class PLUGIN_API GeneralOperatorPluginInfo
{
  public:
    virtual ~GeneralOperatorPluginInfo() {;};
    virtual char *GetName() const = 0;
    virtual char *GetVersion() const = 0;
    virtual char *GetID() const = 0;
    virtual bool  EnabledByDefault() const { return true; }
};

class PLUGIN_API CommonOperatorPluginInfo : public virtual GeneralOperatorPluginInfo
{
  public:
    virtual AttributeSubject *AllocAttributes() = 0;
    virtual void CopyAttributes(AttributeSubject *to,
                                AttributeSubject *from) = 0;
};

class PLUGIN_API GUIOperatorPluginInfo : public virtual CommonOperatorPluginInfo
{
  public:
    virtual const char *GetMenuName() const = 0;
    virtual QvisPostableWindowObserver *CreatePluginWindow(int type,
        AttributeSubject *attr, QvisNotepadArea *notepad) = 0;
    virtual QvisWizard *CreatePluginWizard(AttributeSubject *attr,
        QWidget *parent, const char *name =0)
    {
        return 0;
    }
    virtual const char **XPMIconData() const { return 0; }
};

class PLUGIN_API ViewerOperatorPluginInfo : public virtual CommonOperatorPluginInfo
{
  public:
    virtual AttributeSubject *GetClientAtts() = 0;
    virtual AttributeSubject *GetDefaultAtts() = 0;
    virtual void SetClientAtts(AttributeSubject *atts) = 0;
    virtual void GetClientAtts(AttributeSubject *atts) = 0;

    virtual void InitializeOperatorAtts(AttributeSubject *atts,
                                        const ViewerPlot *plot,
                                        const bool fromDefault) = 0;
    virtual const char **XPMIconData() const { return 0; }
    virtual bool Removeable() const { return true; }
    virtual bool Moveable() const { return true; }
    virtual bool AllowsSubsequentOperators() const { return true; }
};

class PLUGIN_API EngineOperatorPluginInfo : public virtual CommonOperatorPluginInfo
{
  public:
    virtual avtPluginFilter *AllocAvtPluginFilter() = 0;
};

class PLUGIN_API ScriptingOperatorPluginInfo : public virtual CommonOperatorPluginInfo
{
  public:
    virtual void InitializePlugin(AttributeSubject *subj, void *data) = 0;
    virtual void *GetMethodTable(int *nMethods) = 0;
    virtual char *GetLogString() = 0;
    virtual bool TypesMatch(void *pyobject) { return false; }
    virtual void SetDefaults(const AttributeSubject *) = 0;
};

#endif
