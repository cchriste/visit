// ************************************************************************* //
//  File: DisplacePluginInfo.h
// ************************************************************************* //

#ifndef DISPLACE_PLUGIN_INFO_H
#define DISPLACE_PLUGIN_INFO_H
#include <OperatorPluginInfo.h>
#include <operator_plugin_exports.h>

class DisplaceAttributes;

// ****************************************************************************
//  Class: DisplacePluginInfo
//
//  Purpose:
//    Five classes that provide all the information about an Displace operator
//
//  Programmer: whitlocb -- generated by xml2info
//  Creation:   Fri Jan 6 18:12:48 PST 2006
//
//  Modifications:
//
// ****************************************************************************

class DisplaceGeneralPluginInfo : public virtual GeneralOperatorPluginInfo
{
  public:
    virtual char *GetName() const;
    virtual char *GetVersion() const;
    virtual char *GetID() const;
    virtual bool  EnabledByDefault() const;
};

class DisplaceCommonPluginInfo : public virtual CommonOperatorPluginInfo, public virtual DisplaceGeneralPluginInfo
{
  public:
    virtual AttributeSubject *AllocAttributes();
    virtual void CopyAttributes(AttributeSubject *to, AttributeSubject *from);
};

class DisplaceGUIPluginInfo : public virtual GUIOperatorPluginInfo, public virtual DisplaceCommonPluginInfo
{
  public:
    virtual const char *GetMenuName() const;
    virtual QvisPostableWindowObserver *CreatePluginWindow(int type,
        AttributeSubject *attr, QvisNotepadArea *notepad);
    virtual const char **XPMIconData() const;
};

class DisplaceViewerPluginInfo : public virtual ViewerOperatorPluginInfo, public virtual DisplaceCommonPluginInfo
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
    static DisplaceAttributes *defaultAtts;
    static DisplaceAttributes *clientAtts;
};

class DisplaceEnginePluginInfo : public virtual EngineOperatorPluginInfo, public virtual DisplaceCommonPluginInfo
{
  public:
    virtual avtPluginFilter *AllocAvtPluginFilter();
};

class DisplaceScriptingPluginInfo : public virtual ScriptingOperatorPluginInfo, public virtual DisplaceCommonPluginInfo
{
  public:
    virtual void InitializePlugin(AttributeSubject *subj, void *data);
    virtual void *GetMethodTable(int *nMethods);
    virtual bool TypesMatch(void *pyobject);
    virtual char *GetLogString();
    virtual void SetDefaults(const AttributeSubject *atts);
};

#endif
