// ************************************************************************* //
//  File: ThreeSlicePluginInfo.h
// ************************************************************************* //

#ifndef THREESLICE_PLUGIN_INFO_H
#define THREESLICE_PLUGIN_INFO_H
#include <OperatorPluginInfo.h>
#include <operator_plugin_exports.h>

class ThreeSliceAttributes;

// ****************************************************************************
//  Class: ThreeSlicePluginInfo
//
//  Purpose:
//    Five classes that provide all the information about an ThreeSlice operator
//
//  Programmer: haddox1 -- generated by xml2info
//  Creation:   Mon Jun 9 13:18:31 PST 2003
//
//  Modifications:
//
// ****************************************************************************

class ThreeSliceGeneralPluginInfo : public virtual GeneralOperatorPluginInfo
{
  public:
    virtual char *GetName() const;
    virtual char *GetVersion() const;
    virtual char *GetID() const;
};

class ThreeSliceCommonPluginInfo : public virtual CommonOperatorPluginInfo, public virtual ThreeSliceGeneralPluginInfo
{
  public:
    virtual AttributeSubject *AllocAttributes();
    virtual void CopyAttributes(AttributeSubject *to, AttributeSubject *from);
};

class ThreeSliceGUIPluginInfo : public virtual GUIOperatorPluginInfo, public virtual ThreeSliceCommonPluginInfo
{
  public:
    virtual const char *GetMenuName() const;
    virtual QvisPostableWindowObserver *CreatePluginWindow(int type,
        AttributeSubject *attr, QvisNotepadArea *notepad);
    virtual const char **XPMIconData() const;
};

class ThreeSliceViewerPluginInfo : public virtual ViewerOperatorPluginInfo, public virtual ThreeSliceCommonPluginInfo
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
    static ThreeSliceAttributes *defaultAtts;
    static ThreeSliceAttributes *clientAtts;
};

class ThreeSliceEnginePluginInfo : public virtual EngineOperatorPluginInfo, public virtual ThreeSliceCommonPluginInfo
{
  public:
    virtual avtPluginFilter *AllocAvtPluginFilter();
};

class ThreeSliceScriptingPluginInfo : public virtual ScriptingOperatorPluginInfo, public virtual ThreeSliceCommonPluginInfo
{
  public:
    virtual void InitializePlugin(AttributeSubject *subj, FILE *log);
    virtual void *GetMethodTable(int *nMethods);
    virtual bool TypesMatch(void *pyobject);
    virtual void SetLogging(bool val);
    virtual void SetDefaults(const AttributeSubject *atts);
};

#endif
