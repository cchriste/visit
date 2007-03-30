// ************************************************************************* //
//  File: SphereSlicePluginInfo.h
// ************************************************************************* //

#ifndef SPHERESLICE_PLUGIN_INFO_H
#define SPHERESLICE_PLUGIN_INFO_H
#include <OperatorPluginInfo.h>
#include <operator_plugin_exports.h>

class SphereSliceAttributes;

// ****************************************************************************
//  Class: SphereSlicePluginInfo
//
//  Purpose:
//    Five classes that provide all the information about an SphereSlice operator
//
//  Programmer: whitlocb -- generated by xml2info
//  Creation:   Thu Mar 13 15:31:30 PST 2003
//
//  Modifications:
//
// ****************************************************************************

class SphereSliceGeneralPluginInfo : public virtual GeneralOperatorPluginInfo
{
  public:
    virtual char *GetName() const;
    virtual char *GetVersion() const;
    virtual char *GetID() const;
};

class SphereSliceCommonPluginInfo : public virtual CommonOperatorPluginInfo, public virtual SphereSliceGeneralPluginInfo
{
  public:
    virtual AttributeSubject *AllocAttributes();
    virtual void CopyAttributes(AttributeSubject *to, AttributeSubject *from);
};

class SphereSliceGUIPluginInfo : public virtual GUIOperatorPluginInfo, public virtual SphereSliceCommonPluginInfo
{
  public:
    virtual const char *GetMenuName() const;
    virtual QvisPostableWindowObserver *CreatePluginWindow(int type,
        AttributeSubject *attr, QvisNotepadArea *notepad);
    virtual const char **XPMIconData() const;
};

class SphereSliceViewerPluginInfo : public virtual ViewerOperatorPluginInfo, public virtual SphereSliceCommonPluginInfo
{
  public:
    virtual AttributeSubject *GetClientAtts();
    virtual AttributeSubject *GetDefaultAtts();
    virtual void SetClientAtts(AttributeSubject *atts);
    virtual void GetClientAtts(AttributeSubject *atts);

    virtual void InitializeOperatorAtts(AttributeSubject *atts,
                                        const ViewerPlot *plot);
    virtual const char **XPMIconData() const;

    static void InitializeGlobalObjects();
  private:
    static SphereSliceAttributes *defaultAtts;
    static SphereSliceAttributes *clientAtts;
};

class SphereSliceEnginePluginInfo : public virtual EngineOperatorPluginInfo, public virtual SphereSliceCommonPluginInfo
{
  public:
    virtual avtPluginFilter *AllocAvtPluginFilter();
};

class SphereSliceScriptingPluginInfo : public virtual ScriptingOperatorPluginInfo, public virtual SphereSliceCommonPluginInfo
{
  public:
    virtual void InitializePlugin(AttributeSubject *subj, FILE *log);
    virtual void *GetMethodTable(int *nMethods);
    virtual bool TypesMatch(void *pyobject);
    virtual void SetLogging(bool val);
    virtual void SetDefaults(const AttributeSubject *atts);
};

#endif
