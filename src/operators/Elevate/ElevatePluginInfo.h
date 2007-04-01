// ************************************************************************* //
//  File: ElevatePluginInfo.h
// ************************************************************************* //

#ifndef ELEVATE_PLUGIN_INFO_H
#define ELEVATE_PLUGIN_INFO_H
#include <OperatorPluginInfo.h>
#include <operator_plugin_exports.h>

class ElevateAttributes;

// ****************************************************************************
//  Class: ElevatePluginInfo
//
//  Purpose:
//    Five classes that provide all the information about an Elevate operator
//
//  Programmer: childs -- generated by xml2info
//  Creation:   Tue Feb 1 11:37:30 PDT 2005
//
//  Modifications:
//
// ****************************************************************************

class ElevateGeneralPluginInfo : public virtual GeneralOperatorPluginInfo
{
  public:
    virtual char *GetName() const;
    virtual char *GetVersion() const;
    virtual char *GetID() const;
    virtual bool  EnabledByDefault() const;
};

class ElevateCommonPluginInfo : public virtual CommonOperatorPluginInfo, public virtual ElevateGeneralPluginInfo
{
  public:
    virtual AttributeSubject *AllocAttributes();
    virtual void CopyAttributes(AttributeSubject *to, AttributeSubject *from);
};

class ElevateGUIPluginInfo : public virtual GUIOperatorPluginInfo, public virtual ElevateCommonPluginInfo
{
  public:
    virtual const char *GetMenuName() const;
    virtual QvisPostableWindowObserver *CreatePluginWindow(int type,
        AttributeSubject *attr, QvisNotepadArea *notepad);
    virtual const char **XPMIconData() const;
};

class ElevateViewerPluginInfo : public virtual ViewerOperatorPluginInfo, public virtual ElevateCommonPluginInfo
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
    static ElevateAttributes *defaultAtts;
    static ElevateAttributes *clientAtts;
};

class ElevateEnginePluginInfo : public virtual EngineOperatorPluginInfo, public virtual ElevateCommonPluginInfo
{
  public:
    virtual avtPluginFilter *AllocAvtPluginFilter();
};

class ElevateScriptingPluginInfo : public virtual ScriptingOperatorPluginInfo, public virtual ElevateCommonPluginInfo
{
  public:
    virtual void InitializePlugin(AttributeSubject *subj, FILE *log);
    virtual void *GetMethodTable(int *nMethods);
    virtual bool TypesMatch(void *pyobject);
    virtual void SetLogging(bool val);
    virtual void SetDefaults(const AttributeSubject *atts);
};

#endif
