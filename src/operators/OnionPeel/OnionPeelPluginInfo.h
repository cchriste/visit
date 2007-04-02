// ************************************************************************* //
//  File: OnionPeelPluginInfo.h
// ************************************************************************* //

#ifndef ONIONPEEL_PLUGIN_INFO_H
#define ONIONPEEL_PLUGIN_INFO_H
#include <OperatorPluginInfo.h>
#include <operator_plugin_exports.h>

class OnionPeelAttributes;

// ****************************************************************************
//  Class: OnionPeelPluginInfo
//
//  Purpose:
//    Five classes that provide all the information about an OnionPeel operator
//
//  Programmer: whitlocb -- generated by xml2info
//  Creation:   Fri Jan 6 18:13:24 PST 2006
//
//  Modifications:
//
// ****************************************************************************

class OnionPeelGeneralPluginInfo : public virtual GeneralOperatorPluginInfo
{
  public:
    virtual char *GetName() const;
    virtual char *GetVersion() const;
    virtual char *GetID() const;
    virtual bool  EnabledByDefault() const;
};

class OnionPeelCommonPluginInfo : public virtual CommonOperatorPluginInfo, public virtual OnionPeelGeneralPluginInfo
{
  public:
    virtual AttributeSubject *AllocAttributes();
    virtual void CopyAttributes(AttributeSubject *to, AttributeSubject *from);
};

class OnionPeelGUIPluginInfo : public virtual GUIOperatorPluginInfo, public virtual OnionPeelCommonPluginInfo
{
  public:
    virtual const char *GetMenuName() const;
    virtual QvisPostableWindowObserver *CreatePluginWindow(int type,
        AttributeSubject *attr, QvisNotepadArea *notepad);
    virtual const char **XPMIconData() const;
};

class OnionPeelViewerPluginInfo : public virtual ViewerOperatorPluginInfo, public virtual OnionPeelCommonPluginInfo
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
    static OnionPeelAttributes *defaultAtts;
    static OnionPeelAttributes *clientAtts;
};

class OnionPeelEnginePluginInfo : public virtual EngineOperatorPluginInfo, public virtual OnionPeelCommonPluginInfo
{
  public:
    virtual avtPluginFilter *AllocAvtPluginFilter();
};

class OnionPeelScriptingPluginInfo : public virtual ScriptingOperatorPluginInfo, public virtual OnionPeelCommonPluginInfo
{
  public:
    virtual void InitializePlugin(AttributeSubject *subj, void *data);
    virtual void *GetMethodTable(int *nMethods);
    virtual bool TypesMatch(void *pyobject);
    virtual char *GetLogString();
    virtual void SetDefaults(const AttributeSubject *atts);
};

#endif
