// ************************************************************************* //
//  File: RevolvePluginInfo.h
// ************************************************************************* //

#ifndef REVOLVE_PLUGIN_INFO_H
#define REVOLVE_PLUGIN_INFO_H
#include <OperatorPluginInfo.h>
#include <operator_plugin_exports.h>

class RevolveAttributes;

// ****************************************************************************
//  Class: RevolvePluginInfo
//
//  Purpose:
//    Five classes that provide all the information about an Revolve operator
//
//  Programmer: childs -- generated by xml2info
//  Creation:   Wed Dec 11 14:17:27 PST 2002
//
//  Modifications:
//
// ****************************************************************************

class RevolveGeneralPluginInfo : public virtual GeneralOperatorPluginInfo
{
  public:
    virtual char *GetName() const;
    virtual char *GetVersion() const;
    virtual char *GetID() const;
};

class RevolveCommonPluginInfo : public virtual CommonOperatorPluginInfo, public virtual RevolveGeneralPluginInfo
{
  public:
    virtual AttributeSubject *AllocAttributes();
    virtual void CopyAttributes(AttributeSubject *to, AttributeSubject *from);
};

class RevolveGUIPluginInfo : public virtual GUIOperatorPluginInfo, public virtual RevolveCommonPluginInfo
{
  public:
    virtual const char *GetMenuName() const;
    virtual QvisPostableWindowObserver *CreatePluginWindow(int type,
        AttributeSubject *attr, QvisNotepadArea *notepad);
};

class RevolveViewerPluginInfo : public virtual ViewerOperatorPluginInfo, public virtual RevolveCommonPluginInfo
{
  public:
    virtual AttributeSubject *GetClientAtts();
    virtual AttributeSubject *GetDefaultAtts();
    virtual void SetClientAtts(AttributeSubject *atts);
    virtual void GetClientAtts(AttributeSubject *atts);

    virtual void InitializeOperatorAtts(AttributeSubject *atts,
                                        const ViewerPlot *plot);

    static void InitializeGlobalObjects();
  private:
    static RevolveAttributes *defaultAtts;
    static RevolveAttributes *clientAtts;
};

class RevolveEnginePluginInfo : public virtual EngineOperatorPluginInfo, public virtual RevolveCommonPluginInfo
{
  public:
    virtual avtPluginFilter *AllocAvtPluginFilter();
};

class RevolveScriptingPluginInfo : public virtual ScriptingOperatorPluginInfo, public virtual RevolveCommonPluginInfo
{
  public:
    virtual void InitializePlugin(AttributeSubject *subj, FILE *log);
    virtual void *GetMethodTable(int *nMethods);
    virtual bool TypesMatch(void *pyobject);
    virtual void SetLogging(bool val);
    virtual void SetDefaults(const AttributeSubject *atts);
};

#endif
