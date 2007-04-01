// ************************************************************************* //
//  File: InverseGhostZonePluginInfo.h
// ************************************************************************* //

#ifndef INVERSEGHOSTZONE_PLUGIN_INFO_H
#define INVERSEGHOSTZONE_PLUGIN_INFO_H
#include <OperatorPluginInfo.h>
#include <operator_plugin_exports.h>

class InverseGhostZoneAttributes;

// ****************************************************************************
//  Class: InverseGhostZonePluginInfo
//
//  Purpose:
//    Five classes that provide all the information about an InverseGhostZone operator
//
//  Programmer: childs -- generated by xml2info
//  Creation:   Thu Jan 8 09:57:12 PDT 2004
//
//  Modifications:
//
// ****************************************************************************

class InverseGhostZoneGeneralPluginInfo : public virtual GeneralOperatorPluginInfo
{
  public:
    virtual char *GetName() const;
    virtual char *GetVersion() const;
    virtual char *GetID() const;
    virtual bool  EnabledByDefault() const;
};

class InverseGhostZoneCommonPluginInfo : public virtual CommonOperatorPluginInfo, public virtual InverseGhostZoneGeneralPluginInfo
{
  public:
    virtual AttributeSubject *AllocAttributes();
    virtual void CopyAttributes(AttributeSubject *to, AttributeSubject *from);
};

class InverseGhostZoneGUIPluginInfo : public virtual GUIOperatorPluginInfo, public virtual InverseGhostZoneCommonPluginInfo
{
  public:
    virtual const char *GetMenuName() const;
    virtual QvisPostableWindowObserver *CreatePluginWindow(int type,
        AttributeSubject *attr, QvisNotepadArea *notepad);
    virtual const char **XPMIconData() const;
};

class InverseGhostZoneViewerPluginInfo : public virtual ViewerOperatorPluginInfo, public virtual InverseGhostZoneCommonPluginInfo
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
    static InverseGhostZoneAttributes *defaultAtts;
    static InverseGhostZoneAttributes *clientAtts;
};

class InverseGhostZoneEnginePluginInfo : public virtual EngineOperatorPluginInfo, public virtual InverseGhostZoneCommonPluginInfo
{
  public:
    virtual avtPluginFilter *AllocAvtPluginFilter();
};

class InverseGhostZoneScriptingPluginInfo : public virtual ScriptingOperatorPluginInfo, public virtual InverseGhostZoneCommonPluginInfo
{
  public:
    virtual void InitializePlugin(AttributeSubject *subj, FILE *log);
    virtual void *GetMethodTable(int *nMethods);
    virtual bool TypesMatch(void *pyobject);
    virtual void SetLogging(bool val);
    virtual void SetDefaults(const AttributeSubject *atts);
};

#endif
