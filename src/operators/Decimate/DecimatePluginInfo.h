// ************************************************************************* //
//  File: DecimatePluginInfo.h
// ************************************************************************* //

#ifndef DECIMATE_PLUGIN_INFO_H
#define DECIMATE_PLUGIN_INFO_H
#include <OperatorPluginInfo.h>
#include <operator_plugin_exports.h>

class DecimateAttributes;

// ****************************************************************************
//  Class: DecimatePluginInfo
//
//  Purpose:
//    Five classes that provide all the information about an Decimate operator
//
//  Programmer: whitlocb -- generated by xml2info
//  Creation:   Fri Jan 6 18:12:40 PST 2006
//
//  Modifications:
//
// ****************************************************************************

class DecimateGeneralPluginInfo : public virtual GeneralOperatorPluginInfo
{
  public:
    virtual char *GetName() const;
    virtual char *GetVersion() const;
    virtual char *GetID() const;
    virtual bool  EnabledByDefault() const;
};

class DecimateCommonPluginInfo : public virtual CommonOperatorPluginInfo, public virtual DecimateGeneralPluginInfo
{
  public:
    virtual AttributeSubject *AllocAttributes();
    virtual void CopyAttributes(AttributeSubject *to, AttributeSubject *from);
};

class DecimateGUIPluginInfo : public virtual GUIOperatorPluginInfo, public virtual DecimateCommonPluginInfo
{
  public:
    virtual const char *GetMenuName() const;
    virtual QvisPostableWindowObserver *CreatePluginWindow(int type,
        AttributeSubject *attr, QvisNotepadArea *notepad);
};

class DecimateViewerPluginInfo : public virtual ViewerOperatorPluginInfo, public virtual DecimateCommonPluginInfo
{
  public:
    virtual AttributeSubject *GetClientAtts();
    virtual AttributeSubject *GetDefaultAtts();
    virtual void SetClientAtts(AttributeSubject *atts);
    virtual void GetClientAtts(AttributeSubject *atts);

    virtual void InitializeOperatorAtts(AttributeSubject *atts,
                                        const ViewerPlot *plot,
                                        const bool fromDefault);

    static void InitializeGlobalObjects();
  private:
    static DecimateAttributes *defaultAtts;
    static DecimateAttributes *clientAtts;
};

class DecimateEnginePluginInfo : public virtual EngineOperatorPluginInfo, public virtual DecimateCommonPluginInfo
{
  public:
    virtual avtPluginFilter *AllocAvtPluginFilter();
};

class DecimateScriptingPluginInfo : public virtual ScriptingOperatorPluginInfo, public virtual DecimateCommonPluginInfo
{
  public:
    virtual void InitializePlugin(AttributeSubject *subj, void *data);
    virtual void *GetMethodTable(int *nMethods);
    virtual bool TypesMatch(void *pyobject);
    virtual char *GetLogString();
    virtual void SetDefaults(const AttributeSubject *atts);
};

#endif
