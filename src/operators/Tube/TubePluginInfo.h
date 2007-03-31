// ************************************************************************* //
//  File: TubePluginInfo.h
// ************************************************************************* //

#ifndef TUBE_PLUGIN_INFO_H
#define TUBE_PLUGIN_INFO_H
#include <OperatorPluginInfo.h>
#include <operator_plugin_exports.h>

class TubeAttributes;

// ****************************************************************************
//  Class: TubePluginInfo
//
//  Purpose:
//    Five classes that provide all the information about an Tube operator
//
//  Programmer: kbonnell -- generated by xml2info
//  Creation:   Tue Sep 9 16:06:21 PST 2003
//
//  Modifications:
//
// ****************************************************************************

class TubeGeneralPluginInfo : public virtual GeneralOperatorPluginInfo
{
  public:
    virtual char *GetName() const;
    virtual char *GetVersion() const;
    virtual char *GetID() const;
};

class TubeCommonPluginInfo : public virtual CommonOperatorPluginInfo, public virtual TubeGeneralPluginInfo
{
  public:
    virtual AttributeSubject *AllocAttributes();
    virtual void CopyAttributes(AttributeSubject *to, AttributeSubject *from);
};

class TubeGUIPluginInfo : public virtual GUIOperatorPluginInfo, public virtual TubeCommonPluginInfo
{
  public:
    virtual const char *GetMenuName() const;
    virtual QvisPostableWindowObserver *CreatePluginWindow(int type,
        AttributeSubject *attr, QvisNotepadArea *notepad);
    virtual const char **XPMIconData() const;
};

class TubeViewerPluginInfo : public virtual ViewerOperatorPluginInfo, public virtual TubeCommonPluginInfo
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
    static TubeAttributes *defaultAtts;
    static TubeAttributes *clientAtts;
};

class TubeEnginePluginInfo : public virtual EngineOperatorPluginInfo, public virtual TubeCommonPluginInfo
{
  public:
    virtual avtPluginFilter *AllocAvtPluginFilter();
};

class TubeScriptingPluginInfo : public virtual ScriptingOperatorPluginInfo, public virtual TubeCommonPluginInfo
{
  public:
    virtual void InitializePlugin(AttributeSubject *subj, FILE *log);
    virtual void *GetMethodTable(int *nMethods);
    virtual bool TypesMatch(void *pyobject);
    virtual void SetLogging(bool val);
    virtual void SetDefaults(const AttributeSubject *atts);
};

#endif
