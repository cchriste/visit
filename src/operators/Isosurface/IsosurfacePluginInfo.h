// ************************************************************************* //
//  File: IsosurfacePluginInfo.h
// ************************************************************************* //

#ifndef ISOSURFACE_PLUGIN_INFO_H
#define ISOSURFACE_PLUGIN_INFO_H
#include <OperatorPluginInfo.h>
#include <operator_plugin_exports.h>

class IsosurfaceAttributes;

// ****************************************************************************
//  Class: IsosurfacePluginInfo
//
//  Purpose:
//    Five classes that provide all the information about an Isosurface operator
//
//  Programmer: whitlocb -- generated by xml2info
//  Creation:   Fri Jan 6 18:13:07 PST 2006
//
//  Modifications:
//
// ****************************************************************************

class IsosurfaceGeneralPluginInfo : public virtual GeneralOperatorPluginInfo
{
  public:
    virtual char *GetName() const;
    virtual char *GetVersion() const;
    virtual char *GetID() const;
    virtual bool  EnabledByDefault() const;
};

class IsosurfaceCommonPluginInfo : public virtual CommonOperatorPluginInfo, public virtual IsosurfaceGeneralPluginInfo
{
  public:
    virtual AttributeSubject *AllocAttributes();
    virtual void CopyAttributes(AttributeSubject *to, AttributeSubject *from);
};

class IsosurfaceGUIPluginInfo : public virtual GUIOperatorPluginInfo, public virtual IsosurfaceCommonPluginInfo
{
  public:
    virtual const char *GetMenuName() const;
    virtual QvisPostableWindowObserver *CreatePluginWindow(int type,
        AttributeSubject *attr, QvisNotepadArea *notepad);
    virtual const char **XPMIconData() const;
};

class IsosurfaceViewerPluginInfo : public virtual ViewerOperatorPluginInfo, public virtual IsosurfaceCommonPluginInfo
{
  public:
    virtual AttributeSubject *GetClientAtts();
    virtual AttributeSubject *GetDefaultAtts();
    virtual void SetClientAtts(AttributeSubject *atts);
    virtual void GetClientAtts(AttributeSubject *atts);

    virtual void InitializeOperatorAtts(AttributeSubject *atts,
                                        const ViewerPlot *plot,
                                        const bool fromDefault);
    virtual const char **XPMIconData() const;

    static void InitializeGlobalObjects();
  private:
    static IsosurfaceAttributes *defaultAtts;
    static IsosurfaceAttributes *clientAtts;
};

class IsosurfaceEnginePluginInfo : public virtual EngineOperatorPluginInfo, public virtual IsosurfaceCommonPluginInfo
{
  public:
    virtual avtPluginFilter *AllocAvtPluginFilter();
};

class IsosurfaceScriptingPluginInfo : public virtual ScriptingOperatorPluginInfo, public virtual IsosurfaceCommonPluginInfo
{
  public:
    virtual void InitializePlugin(AttributeSubject *subj, void *data);
    virtual void *GetMethodTable(int *nMethods);
    virtual bool TypesMatch(void *pyobject);
    virtual char *GetLogString();
    virtual void SetDefaults(const AttributeSubject *atts);
};

#endif
