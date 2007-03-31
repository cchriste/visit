// ************************************************************************* //
//  File: BoxPluginInfo.h
// ************************************************************************* //

#ifndef BOX_PLUGIN_INFO_H
#define BOX_PLUGIN_INFO_H
#include <OperatorPluginInfo.h>
#include <operator_plugin_exports.h>

class BoxAttributes;

// ****************************************************************************
//  Class: BoxPluginInfo
//
//  Purpose:
//    Five classes that provide all the information about an Box operator
//
//  Programmer: kbonnell -- generated by xml2info
//  Creation:   Tue Sep 9 15:59:05 PST 2003
//
//  Modifications:
//
// ****************************************************************************

class BoxGeneralPluginInfo : public virtual GeneralOperatorPluginInfo
{
  public:
    virtual char *GetName() const;
    virtual char *GetVersion() const;
    virtual char *GetID() const;
};

class BoxCommonPluginInfo : public virtual CommonOperatorPluginInfo, public virtual BoxGeneralPluginInfo
{
  public:
    virtual AttributeSubject *AllocAttributes();
    virtual void CopyAttributes(AttributeSubject *to, AttributeSubject *from);
};

class BoxGUIPluginInfo : public virtual GUIOperatorPluginInfo, public virtual BoxCommonPluginInfo
{
  public:
    virtual const char *GetMenuName() const;
    virtual QvisPostableWindowObserver *CreatePluginWindow(int type,
        AttributeSubject *attr, QvisNotepadArea *notepad);
    virtual const char **XPMIconData() const;
};

class BoxViewerPluginInfo : public virtual ViewerOperatorPluginInfo, public virtual BoxCommonPluginInfo
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
    static BoxAttributes *defaultAtts;
    static BoxAttributes *clientAtts;
};

class BoxEnginePluginInfo : public virtual EngineOperatorPluginInfo, public virtual BoxCommonPluginInfo
{
  public:
    virtual avtPluginFilter *AllocAvtPluginFilter();
};

class BoxScriptingPluginInfo : public virtual ScriptingOperatorPluginInfo, public virtual BoxCommonPluginInfo
{
  public:
    virtual void InitializePlugin(AttributeSubject *subj, FILE *log);
    virtual void *GetMethodTable(int *nMethods);
    virtual bool TypesMatch(void *pyobject);
    virtual void SetLogging(bool val);
    virtual void SetDefaults(const AttributeSubject *atts);
};

#endif
