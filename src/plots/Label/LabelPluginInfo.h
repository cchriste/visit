// ************************************************************************* //
//                               LabelPluginInfo.h                           //
// ************************************************************************* //

#ifndef LABEL_PLUGIN_INFO_H
#define LABEL_PLUGIN_INFO_H
#include <PlotPluginInfo.h>
#include <plot_plugin_exports.h>

class LabelAttributes;

// ****************************************************************************
//  Class: LabelPluginInfo
//
//  Purpose:
//    Five classes that provide all the information about a Label
//    plot plugin.  The information is broken up into five classes since
//    portions of it are only relevant to particular components within
//    visit.  There is the general information which all the components
//    are interested in, the gui information which the gui is interested in,
//    the viewer information which the viewer is interested in, the
//    engine information which the engine is interested in, and finally a.
//    scripting portion that enables the Python VisIt extension to use the
//    plugin.
//
//  Programmer: whitlocb -- generated by xml2info
//  Creation:   Mon Jan 9 09:13:15 PDT 2006
//
// ****************************************************************************

class LabelGeneralPluginInfo: public virtual GeneralPlotPluginInfo
{
  public:
    virtual char *GetName() const;
    virtual char *GetVersion() const;
    virtual char *GetID() const;
    virtual bool  EnabledByDefault() const;
};

class LabelCommonPluginInfo : public virtual CommonPlotPluginInfo, public virtual LabelGeneralPluginInfo
{
  public:
    virtual AttributeSubject *AllocAttributes();
    virtual void CopyAttributes(AttributeSubject *to, AttributeSubject *from);
};

class LabelGUIPluginInfo: public virtual GUIPlotPluginInfo, public virtual LabelCommonPluginInfo
{
  public:
    virtual const char *GetMenuName() const;
    virtual int GetVariableTypes() const;
    virtual QvisPostableWindowObserver *CreatePluginWindow(int type,
        AttributeSubject *attr, QvisNotepadArea *notepad);
    virtual const char **XPMIconData() const;
};

class LabelViewerPluginInfo: public virtual ViewerPlotPluginInfo, public virtual LabelCommonPluginInfo
{
  public:
    virtual AttributeSubject *GetClientAtts();
    virtual AttributeSubject *GetDefaultAtts();
    virtual void SetClientAtts(AttributeSubject *atts);
    virtual void GetClientAtts(AttributeSubject *atts);

    virtual avtPlot *AllocAvtPlot();

    virtual void InitializePlotAtts(AttributeSubject *atts,
        const avtDatabaseMetaData *md,
        const char *variableName);

    virtual void ReInitializePlotAtts(AttributeSubject *atts,
        const avtDatabaseMetaData *md,
        const char *variableName);

    virtual void ResetPlotAtts(AttributeSubject *atts,
        const avtDatabaseMetaData *md,
        const char *variableName);

    virtual const char **XPMIconData() const;
    virtual int GetVariableTypes() const;

    static void InitializeGlobalObjects();
  private:
    static LabelAttributes *defaultAtts;
    static LabelAttributes *clientAtts;

    virtual void PrivateSetPlotAtts(AttributeSubject *atts,
        const avtDatabaseMetaData *md,
        const char *variableName);
};

class LabelEnginePluginInfo: public virtual EnginePlotPluginInfo, public virtual LabelCommonPluginInfo
{
  public:
    virtual avtPlot *AllocAvtPlot();
};

class LabelScriptingPluginInfo : public virtual ScriptingPlotPluginInfo, public virtual LabelCommonPluginInfo
{
  public:
    virtual void InitializePlugin(AttributeSubject *subj, void *data);
    virtual void *GetMethodTable(int *nMethods);
    virtual bool TypesMatch(void *pyobject);
    virtual char *GetLogString();
    virtual void SetDefaults(const AttributeSubject *atts);
};

#endif
