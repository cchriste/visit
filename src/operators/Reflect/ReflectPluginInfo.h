// ************************************************************************* //
//  File: ReflectPluginInfo.h
// ************************************************************************* //

#ifndef REFLECT_PLUGIN_INFO_H
#define REFLECT_PLUGIN_INFO_H
#include <OperatorPluginInfo.h>
#include <operator_plugin_exports.h>

class ReflectAttributes;

// ****************************************************************************
//  Class: ReflectPluginInfo
//
//  Purpose:
//    Five classes that provide all the information about an Reflect operator
//
//  Programmer: whitlocb -- generated by xml2info
//  Creation:   Fri Jan 6 18:13:35 PST 2006
//
//  Modifications:
//
// ****************************************************************************

class ReflectGeneralPluginInfo : public virtual GeneralOperatorPluginInfo
{
  public:
    virtual char *GetName() const;
    virtual char *GetVersion() const;
    virtual char *GetID() const;
    virtual bool  EnabledByDefault() const;
};

class ReflectCommonPluginInfo : public virtual CommonOperatorPluginInfo, public virtual ReflectGeneralPluginInfo
{
  public:
    virtual AttributeSubject *AllocAttributes();
    virtual void CopyAttributes(AttributeSubject *to, AttributeSubject *from);
};

class ReflectGUIPluginInfo : public virtual GUIOperatorPluginInfo, public virtual ReflectCommonPluginInfo
{
  public:
    virtual const char *GetMenuName() const;
    virtual QvisPostableWindowObserver *CreatePluginWindow(int type,
        AttributeSubject *attr, QvisNotepadArea *notepad);
    virtual const char **XPMIconData() const;
};

class ReflectViewerPluginInfo : public virtual ViewerOperatorPluginInfo, public virtual ReflectCommonPluginInfo
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
    static ReflectAttributes *defaultAtts;
    static ReflectAttributes *clientAtts;
};

class ReflectEnginePluginInfo : public virtual EngineOperatorPluginInfo, public virtual ReflectCommonPluginInfo
{
  public:
    virtual avtPluginFilter *AllocAvtPluginFilter();
};

class ReflectScriptingPluginInfo : public virtual ScriptingOperatorPluginInfo, public virtual ReflectCommonPluginInfo
{
  public:
    virtual void InitializePlugin(AttributeSubject *subj, void *data);
    virtual void *GetMethodTable(int *nMethods);
    virtual bool TypesMatch(void *pyobject);
    virtual char *GetLogString();
    virtual void SetDefaults(const AttributeSubject *atts);
};

#endif
