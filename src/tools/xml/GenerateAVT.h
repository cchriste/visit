#ifndef GENERATE_AVT_H
#define GENERATE_AVT_H

#include "Field.h"
#include <stdlib.h>
#include <sys/types.h>
#include <time.h>

// ****************************************************************************
//  File:  GenerateAVT
//
//  Purpose:
//    Contains a set of classes which override the default implementation
//    to create AVT files (avtPlot, avtFilter) for the plugin.
//
//  Note: This file overrides --
//    Plugin
//
//  Programmer:  Jeremy Meredith
//  Creation:    August 29, 2001
//
//  Modifications:
//    Jeremy Meredith, Wed Feb  6 17:26:42 PST 2002
//    Changed references to plugin->name+"Attributes" to instead refer
//    directly to atts->name.
//
//    Jeremy Meredith, Tue Aug 27 14:32:50 PDT 2002
//    Added mfiles, dbtype, and libs.  Allowed NULL atts.
//
//    Hank Childs, Thu Sep 12 19:34:31 PDT 2002
//    Add string argument to ExecuteData.
//
//    Jeremy Meredith, Thu Oct 17 15:58:29 PDT 2002
//    Added some enhancements for the XML editor.
//
//    Kathleen Bonnell, Wed Oct 23 18:10:26 PDT 2002  
//    Added new plot method ApplyRenderingTransformation. 
//
// ****************************************************************************

// ----------------------------------------------------------------------------
//                             Utility Functions
// ----------------------------------------------------------------------------

QString
CurrentTime()
{
    char *tstr[] = {"PDT", "PST"};
    char s1[10], s2[10], s3[10], tmpbuf[200];
    time_t t;
    char *c = NULL;
    int h,m,s,y;
    t = time(NULL);
    c = asctime(localtime(&t));
    // Read the hour.
    sscanf(c, "%s %s %s %d:%d:%d %d", s1, s2, s3, &h, &m, &s, &y);
    // Reformat the string a little.
    sprintf(tmpbuf, "%s %s %s %02d:%02d:%02d %s %d",
            s1, s2, s3, h, m, s, tstr[h > 12], y);

    return QString(tmpbuf);
}

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
class AVTGeneratorPlugin
{
  public:
    QString name;
    QString type;
    QString label;
    QString version;
    QString vartype;
    QString dbtype;

    vector<QString> cxxflags;
    vector<QString> ldflags;
    vector<QString> libs;
    vector<QString> extensions; // for DB plugins
    bool customgfiles;
    vector<QString> gfiles;     // gui
    bool customsfiles;
    vector<QString> sfiles;     // scripting
    bool customvfiles;
    vector<QString> vfiles;     // viewer
    bool custommfiles;
    vector<QString> mfiles;     // mdserver
    bool customefiles;
    vector<QString> efiles;     // engine
    bool customwfiles;
    vector<QString> wfiles;     // widget

    Attribute *atts;
  public:
    AVTGeneratorPlugin(const QString &n,const QString &l,const QString &t,const QString &vt,const QString &dt,const QString &v,const QString &, const QString &)
        : name(n), type(t), label(l), version(v), vartype(vt), dbtype(dt), atts(NULL)
    {
    }
    void Print(ostream &out)
    {
        out << "Plugin: "<<name<<" (\""<<label<<"\", type="<<type<<") -- version "<<version<< endl;
        if (atts)
            atts->Print(cout);
    }
    void WritePlotHeader(ostream &h)
    {
        if (type!="plot")
        {
            cerr << "Must be of type plot!" << endl;
            return;
        }

        h << "// ************************************************************************* //" << endl;
        h << "//                                 avt"<<name<<"Plot.h                             //" << endl;
        h << "// ************************************************************************* //" << endl;
        h << endl;
        h << "#ifndef AVT_"<<name<<"_PLOT_H" << endl;
        h << "#define AVT_"<<name<<"_PLOT_H" << endl;
        h << endl;
        h << endl;
        h << "#include <avtLegend.h>" << endl;
        h << "#include <avtPlot.h>" << endl;
        h << endl;
        h << "#include <"<<atts->name<<".h>" << endl;
        h << endl;
        h << "class     avt"<<name<<"Filter;" << endl;
        h << endl;
        h << endl;
        h << "// ****************************************************************************" << endl;
        h << "//  Class:  avt"<<name<<"Plot" << endl;
        h << "//" << endl;
        h << "//  Purpose:" << endl;
        h << "//      A concrete type of avtPlot, this is the "<<name<<" plot." << endl;
        h << "//" << endl;
        h << "//  Programmer: "<<getenv("USER")<<" -- generated by xml2info" << endl;
        h << "//  Creation:   "<<CurrentTime()<< endl;
        h << "//" << endl;
        h << "// ****************************************************************************" << endl;
        h << endl;
        h << "YOU SHOULD INHERIT FROM ONE OF THE FOLLOWING PLOT TYPES:" << endl;
        h << "avtImageDataPlot" << endl;
        h << "avtPointDataPlot" << endl;
        h << "avtLineDataPlot" << endl;
        h << "avtSurfaceDataPlot" << endl;
        h << "avtVolumeDataPlot" << endl;
        h << endl;
        h << "class avt"<<name<<"Plot : public avtPointDataPlot" << endl;
        h << "{" << endl;
        h << "  public:" << endl;
        h << "                                avt"<<name<<"Plot();" << endl;
        h << "    virtual                    ~avt"<<name<<"Plot();" << endl;
        h << endl;
        h << "    virtual const char         *GetName(void) { return \""<<name<<"Plot\"; };" << endl;
        h << endl;
        h << "    static avtPlot             *Create();" << endl;
        h << endl;
        h << "    virtual void                SetAtts(const AttributeGroup*);" << endl;
        h << endl;
        h << "  protected:" << endl;
        h << "    "<<atts->name<<"              atts;" << endl;
        h << endl;
        h << "    YOU MUST HAVE SOME SORT OF MAPPER FOR THE PLOT." << endl;
        h << "    avt...Mapper               *myMapper;" << endl;
        h << "    avt"<<name<<"Filter              *"<<name<<"Filter;" << endl;
        h << endl;
        h << "    virtual avtMapper          *GetMapper(void);" << endl;
        h << "    virtual avtDataObject_p     ApplyOperators(avtDataObject_p);" << endl;
        h << "    virtual avtDataObject_p     ApplyRenderingTransformation(avtDataObject_p);"<< endl;
        h << "    virtual void                CustomizeBehavior(void);" << endl;
        h << "    virtual void                CustomizeMapper(avtDataObjectInformation &);" << endl;
        h << endl;
        h << "    virtual avtLegend_p         GetLegend(void) { return NULL; };" << endl;
        h << "};" << endl;
        h << endl;
        h << endl;
        h << "#endif" << endl;
    }
    void WritePlotSource(ostream &c)
    {
        if (type!="plot")
        {
            cerr << "Must be of type plot!" << endl;
            return;
        }
        c << "// ************************************************************************* //" << endl;
        c << "//                             avt"<<name<<"Plot.C                                 //" << endl;
        c << "// ************************************************************************* //" << endl;
        c << endl;
        c << "#include <avt"<<name<<"Plot.h>" << endl;
        c << endl;
        c << "#include <avt"<<name<<"Filter.h>" << endl;
        c << endl;
        c << endl;
        c << "// ****************************************************************************" << endl;
        c << "//  Method: avt"<<name<<"Plot constructor" << endl;
        c << "//" << endl;
        c << "//  Programmer: "<<getenv("USER")<<" -- generated by xml2info" << endl;
        c << "//  Creation:   "<<CurrentTime()<< endl;
        c << "//" << endl;
        c << "// ****************************************************************************" << endl;
        c << endl;
        c << "avt"<<name<<"Plot::avt"<<name<<"Plot()" << endl;
        c << "{" << endl;
        c << "    "<<name<<"Filter = new avt"<<name<<"Filter(ARGS FOR FILTER);" << endl;
        c << "    myMapper   = ....;" << endl;
        c << "}" << endl;
        c << endl;
        c << endl;
        c << "// ****************************************************************************" << endl;
        c << "//  Method: avt"<<name<<"Plot destructor" << endl;
        c << "//" << endl;
        c << "//  Programmer: "<<getenv("USER")<<" -- generated by xml2info" << endl;
        c << "//  Creation:   "<<CurrentTime()<< endl;
        c << "//" << endl;
        c << "// ****************************************************************************" << endl;
        c << endl;
        c << "avt"<<name<<"Plot::~avt"<<name<<"Plot()" << endl;
        c << "{" << endl;
        c << "    if (myMapper != NULL)" << endl;
        c << "    {" << endl;
        c << "        delete myMapper;" << endl;
        c << "        myMapper = NULL;" << endl;
        c << "    }" << endl;
        c << "    if ("<<name<<"Filter != NULL)" << endl;
        c << "    {" << endl;
        c << "        delete "<<name<<"Filter;" << endl;
        c << "        "<<name<<"Filter = NULL;" << endl;
        c << "    }" << endl;
        c << "}" << endl;
        c << endl;
        c << endl;
        c << "// ****************************************************************************" << endl;
        c << "//  Method:  avt"<<name<<"Plot::Create" << endl;
        c << "//" << endl;
        c << "//  Purpose:" << endl;
        c << "//    Call the constructor." << endl;
        c << "//" << endl;
        c << "//  Programmer: "<<getenv("USER")<<" -- generated by xml2info" << endl;
        c << "//  Creation:   "<<CurrentTime()<< endl;
        c << "//" << endl;
        c << "// ****************************************************************************" << endl;
        c << endl;
        c << "avtPlot*" << endl;
        c << "avt"<<name<<"Plot::Create()" << endl;
        c << "{" << endl;
        c << "    return new avt"<<name<<"Plot;" << endl;
        c << "}" << endl;
        c << endl;
        c << endl;
        c << "// ****************************************************************************" << endl;
        c << "//  Method: avt"<<name<<"Plot::GetMapper" << endl;
        c << "//" << endl;
        c << "//  Purpose:" << endl;
        c << "//      Gets a mapper for this plot, it is actually a variable mapper." << endl;
        c << "//" << endl;
        c << "//  Returns:    The variable mapper typed as its base class mapper." << endl;
        c << "//" << endl;
        c << "//  Programmer: "<<getenv("USER")<<" -- generated by xml2info" << endl;
        c << "//  Creation:   "<<CurrentTime()<< endl;
        c << "//" << endl;
        c << "// ****************************************************************************" << endl;
        c << endl;
        c << "avtMapper *" << endl;
        c << "avt"<<name<<"Plot::GetMapper(void)" << endl;
        c << "{" << endl;
        c << "    return myMapper;" << endl;
        c << "}" << endl;
        c << endl;
        c << endl;
        c << "// ****************************************************************************" << endl;
        c << "//  Method: avt"<<name<<"Plot::ApplyOperators" << endl;
        c << "//" << endl;
        c << "//  Purpose:" << endl;
        c << "//      Applies the operators associated with a "<<name<<" plot.  " << endl;
        c << "//      The output from this method is a query-able object." << endl;
        c << "//" << endl;
        c << "//  Arguments:" << endl;
        c << "//      input   The input data object." << endl;
        c << "//" << endl;
        c << "//  Returns:    The data object after the "<<name<<" plot has been applied." << endl;
        c << "//" << endl;
        c << "//  Programmer: "<<getenv("USER")<<" -- generated by xml2info" << endl;
        c << "//  Creation:   "<<CurrentTime()<< endl;
        c << "//" << endl;
        c << "// ****************************************************************************" << endl;
        c << endl;
        c << "avtDataObject_p" << endl;
        c << "avt"<<name<<"Plot::ApplyOperators(avtDataObject_p input)" << endl;
        c << "{" << endl;
        c << "    "<<name<<"Filter->SetInput(input);" << endl;
        c << "    return "<<name<<"Filter->GetOutput();" << endl;
        c << "}" << endl;
        c << endl;
        c << endl;
        c << "// ****************************************************************************" << endl;
        c << "//  Method: avt"<<name<<"Plot::ApplyRenderingTransformation" << endl;
        c << "//" << endl;
        c << "//  Purpose:" << endl;
        c << "//      Applies the rendering transformation associated with a "<<name<<" plot.  " << endl;
        c << "//" << endl;
        c << "//  Arguments:" << endl;
        c << "//      input   The input data object." << endl;
        c << "//" << endl;
        c << "//  Returns:    The data object after the "<<name<<" plot has been applied." << endl;
        c << "//" << endl;
        c << "//  Programmer: "<<getenv("USER")<<" -- generated by xml2info" << endl;
        c << "//  Creation:   "<<CurrentTime()<< endl;
        c << "//" << endl;
        c << "// ****************************************************************************" << endl;
        c << endl;
        c << "avtDataObject_p" << endl;
        c << "avt"<<name<<"Plot::ApplyRenderingOperators(avtDataObject_p input)" << endl;
        c << "{" << endl;
        c << "    "<<name<<"Filter->SetInput(input);" << endl;
        c << "    return "<<name<<"Filter->GetOutput();" << endl;
        c << "}" << endl;
        c << endl;
        c << endl;
        c << "// ****************************************************************************" << endl;
        c << "//  Method: avt"<<name<<"Plot::CustomizeBehavior" << endl;
        c << "//" << endl;
        c << "//  Purpose:" << endl;
        c << "//      Customizes the behavior as appropriate for a "<<name<<" plot.  This includes" << endl;
        c << "//      behavior like shifting towards or away from the screen." << endl;
        c << "//" << endl;
        c << "//  Programmer: "<<getenv("USER")<<" -- generated by xml2info" << endl;
        c << "//  Creation:   "<<CurrentTime()<< endl;
        c << "//" << endl;
        c << "// ****************************************************************************" << endl;
        c << endl;
        c << "void" << endl;
        c << "avt"<<name<<"Plot::CustomizeBehavior(void)" << endl;
        c << "{" << endl;
        c << "    //behavior->SetShiftFactor(0.6);" << endl;
        c << "}" << endl;
        c << endl;
        c << endl;
        c << "// ****************************************************************************" << endl;
        c << "//  Method: avt"<<name<<"Plot::CustomizeMapper" << endl;
        c << "//" << endl;
        c << "//  Purpose:" << endl;
        c << "//      A hook from the base class that allows the plot to change its mapper" << endl;
        c << "//      based on the dataset input. " << endl;
        c << "//" << endl;
        c << "//  Arguments:" << endl;
        c << "//      doi     The data object information." << endl;
        c << "//" << endl;
        c << "//  Programmer: "<<getenv("USER")<<" -- generated by xml2info" << endl;
        c << "//  Creation:   "<<CurrentTime()<< endl;
        c << "//" << endl;
        c << "// ****************************************************************************" << endl;
        c << endl;
        c << "void" << endl;
        c << "avt"<<name<<"Plot::CustomizeMapper(avtDataObjectInformation &doi)" << endl;
        c << "{" << endl;
        c << "/* Example of usage." << endl;
        c << "    int dim = doi.GetAttributes().GetCurrentSpatialDimension();" << endl;
        c << "    if (dim == 2)" << endl;
        c << "    {" << endl;
        c << "    }" << endl;
        c << "    else" << endl;
        c << "    {" << endl;
        c << "    }" << endl;
        c << " */" << endl;
        c << "}" << endl;
        c << endl;
        c << endl;
        c << "// ****************************************************************************" << endl;
        c << "//  Method: avt"<<name<<"Plot::SetAtts" << endl;
        c << "//" << endl;
        c << "//  Purpose:" << endl;
        c << "//      Sets the atts for the "<<name<<" plot." << endl;
        c << "//" << endl;
        c << "//  Arguments:" << endl;
        c << "//      atts    The attributes for this "<<name<<" plot." << endl;
        c << "//" << endl;
        c << "//  Programmer: "<<getenv("USER")<<" -- generated by xml2info" << endl;
        c << "//  Creation:   "<<CurrentTime()<< endl;
        c << "//" << endl;
        c << "// ****************************************************************************" << endl;
        c << endl;
        c << "void" << endl;
        c << "avt"<<name<<"Plot::SetAtts(const AttributeGroup *a)" << endl;
        c << "{" << endl;
        c << "    const "<<atts->name<<" *newAtts = (const "<<atts->name<<" *)a;" << endl;
        c << endl;
        c << "    BASED ON ATTRIBUTE VALUES, CHANGE PARAMETERS IN MAPPER AND FILTER." << endl;
        c << "}" << endl;
    }
    void WriteFilterHeader(ostream &h)
    {
        if (type=="operator")
        {
            h << "// ************************************************************************* //" << endl;
            h << "//  File: avt"<<name<<"Filter.h" << endl;
            h << "// ************************************************************************* //" << endl;
            h << endl;
            h << "#ifndef AVT_"<<name<<"_FILTER_H" << endl;
            h << "#define AVT_"<<name<<"_FILTER_H" << endl;
            h << endl;
            h << endl;
            h << "#include <avtPluginStreamer.h>" << endl;
            h << "#include <"<<atts->name<<".h>" << endl;
            h << endl;
            h << endl;
            h << "class vtkDataSet;" << endl;
            h << endl;
            h << endl;
            h << "// ****************************************************************************" << endl;
            h << "//  Class: avt"<<name<<"Filter" << endl;
            h << "//" << endl;
            h << "//  Purpose:" << endl;
            h << "//      A plugin operator for "<<name<<"." << endl;
            h << "//" << endl;
            h << "//  Programmer: "<<getenv("USER")<<" -- generated by xml2info" << endl;
            h << "//  Creation:   "<<CurrentTime()<< endl;
            h << "//" << endl;
            h << "// ****************************************************************************" << endl;
            h << endl;
            h << "class avt"<<name<<"Filter : public avtPluginStreamer" << endl;
            h << "{" << endl;
            h << "  public:" << endl;
            h << "                         avt"<<name<<"Filter();" << endl;
            h << "    virtual             ~avt"<<name<<"Filter();" << endl;
            h << endl;
            h << "    static avtFilter    *Create();" << endl;
            h << endl;
            h << "    virtual const char  *GetType(void)  { return \"avt"<<name<<"Filter\"; };" << endl;
            h << "    virtual const char  *GetDescription(void)" << endl;
            h << "                             { return \""<<label<<"\"; };" << endl;
            h << endl;
            h << "    virtual void         SetAtts(const AttributeGroup*);" << endl;
            h << "    virtual bool         Equivalent(const AttributeGroup*);" << endl;
            h << endl;
            h << "  protected:" << endl;
            h << "    "<<atts->name<<"   atts;" << endl;
            h << endl;
            h << "    virtual vtkDataSet   *ExecuteData(vtkDataSet *, int, std::string);" << endl;
            h << "};" << endl;
            h << endl;
            h << endl;
            h << "#endif" << endl;
        }
        else if (type=="plot")
        {
            h << "// ************************************************************************* //" << endl;
            h << "//                              avt"<<name<<"Filter.h                              //" << endl;
            h << "// ************************************************************************* //" << endl;
            h << endl;
            h << "#ifndef AVT_"<<name<<"_FILTER_H" << endl;
            h << "#define AVT_"<<name<<"_FILTER_H" << endl;
            h << endl;
            h << endl;
            h << "#include <avtStreamer.h>" << endl;
            h << endl;
            h << endl;
            h << "// ****************************************************************************" << endl;
            h << "//  Class: avt"<<name<<"Filter" << endl;
            h << "//" << endl;
            h << "//  Purpose:" << endl;
            h << "//      This operator is the implied operator associated with an "<<name<<" plot." << endl;
            h << "//" << endl;
            h << "//  Programmer: "<<getenv("USER")<<" -- generated by xml2info" << endl;
            h << "//  Creation:   "<<CurrentTime()<< endl;
            h << "//" << endl;
            h << "// ****************************************************************************" << endl;
            h << endl;
            h << "class avt"<<name<<"Filter : public avtStreamer" << endl;
            h << "{" << endl;
            h << "  public:" << endl;
            h << "                              avt"<<name<<"Filter(YOUR INITIALIZATION ARGS);" << endl;
            h << "    virtual                  ~avt"<<name<<"Filter();" << endl;
            h << endl;
            h << "    virtual const char       *GetType(void)   { return \"avt"<<name<<"Filter\"; };" << endl;
            h << "    virtual const char       *GetDescription(void)" << endl;
            h << "                                  { return \"Performing "<<label<<"\"; };" << endl;
            h << endl;
            h << "    ADD THE SET METHODS YOU NEED HERE" << endl;
            h << endl;
            h << "  protected:" << endl;
            h << "    ADD YOUR DATA MEMBERS HERE" << endl;
            h << endl;
            h << "    virtual vtkDataSet       *ExecuteData(vtkDataSet *, int, std::string);" << endl;
            h << "    virtual void              RefashionDataObjectInfo(void);" << endl;
            h << "};" << endl;
            h << endl;
            h << endl;
            h << "#endif" << endl;
        }
    }
    void WriteFilterSource(ostream &c)
    {
        if (type=="operator")
        {
            c << "// ************************************************************************* //" << endl;
            c << "//  File: avt"<<name<<"Filter.C" << endl;
            c << "// ************************************************************************* //" << endl;
            c << endl;
            c << "#include <avt"<<name<<"Filter.h>" << endl;
            c << endl;
            c << endl;
            c << "// ****************************************************************************" << endl;
            c << "//  Method: avt"<<name<<"Filter constructor" << endl;
            c << "//" << endl;
            c << "//  Programmer: "<<getenv("USER")<<" -- generated by xml2info" << endl;
            c << "//  Creation:   "<<CurrentTime()<< endl;
            c << "//" << endl;
            c << "// ****************************************************************************" << endl;
            c << endl;
            c << "avt"<<name<<"Filter::avt"<<name<<"Filter()" << endl;
            c << "{" << endl;
            c << "}" << endl;
            c << endl;
            c << endl;
            c << "// ****************************************************************************" << endl;
            c << "//  Method: avt"<<name<<"Filter destructor" << endl;
            c << "//" << endl;
            c << "//  Programmer: "<<getenv("USER")<<" -- generated by xml2info" << endl;
            c << "//  Creation:   "<<CurrentTime()<< endl;
            c << "//" << endl;
            c << "//  Modifications:" << endl;
            c << "//" << endl;
            c << "// ****************************************************************************" << endl;
            c << endl;
            c << "avt"<<name<<"Filter::~avt"<<name<<"Filter()" << endl;
            c << "{" << endl;
            c << "}" << endl;
            c << endl;
            c << endl;
            c << "// ****************************************************************************" << endl;
            c << "//  Method:  avt"<<name<<"Filter::Create" << endl;
            c << "//" << endl;
            c << "//  Programmer: "<<getenv("USER")<<" -- generated by xml2info" << endl;
            c << "//  Creation:   "<<CurrentTime()<< endl;
            c << "//" << endl;
            c << "// ****************************************************************************" << endl;
            c << endl;
            c << "avtFilter *" << endl;
            c << "avt"<<name<<"Filter::Create()" << endl;
            c << "{" << endl;
            c << "    return new avt"<<name<<"Filter();" << endl;
            c << "}" << endl;
            c << endl;
            c << endl;
            c << "// ****************************************************************************" << endl;
            c << "//  Method:      avt"<<name<<"Filter::SetAtts" << endl;
            c << "//" << endl;
            c << "//  Purpose:" << endl;
            c << "//      Sets the state of the filter based on the attribute object." << endl;
            c << "//" << endl;
            c << "//  Arguments:" << endl;
            c << "//      a        The attributes to use." << endl;
            c << "//" << endl;
            c << "//  Programmer: "<<getenv("USER")<<" -- generated by xml2info" << endl;
            c << "//  Creation:   "<<CurrentTime()<< endl;
            c << "//" << endl;
            c << "// ****************************************************************************" << endl;
            c << endl;
            c << "void" << endl;
            c << "avt"<<name<<"Filter::SetAtts(const AttributeGroup *a)" << endl;
            c << "{" << endl;
            c << "    atts = *(const "<<atts->name<<"*)a;" << endl;
            c << "}" << endl;
            c << endl;
            c << endl;
            c << "// ****************************************************************************" << endl;
            c << "//  Method: avt"<<name<<"Filter::Equivalent" << endl;
            c << "//" << endl;
            c << "//  Purpose:" << endl;
            c << "//      Returns true if creating a new avt"<<name<<"Filter with the given" << endl;
            c << "//      parameters would result in an equivalent avt"<<name<<"Filter." << endl;
            c << "//" << endl;
            c << "//  Programmer: "<<getenv("USER")<<" -- generated by xml2info" << endl;
            c << "//  Creation:   "<<CurrentTime()<< endl;
            c << "//" << endl;
            c << "// ****************************************************************************" << endl;
            c << endl;
            c << "bool" << endl;
            c << "avt"<<name<<"Filter::Equivalent(const AttributeGroup *a)" << endl;
            c << "{" << endl;
            c << "    return (atts == *("<<atts->name<<"*)a);" << endl;
            c << "}" << endl;
            c << endl;
            c << endl;
            c << "// ****************************************************************************" << endl;
            c << "//  Method: avt"<<name<<"Filter::ExecuteData" << endl;
            c << "//" << endl;
            c << "//  Purpose:" << endl;
            c << "//      Sends the specified input and output through the "<<name<<" filter." << endl;
            c << "//" << endl;
            c << "//  Arguments:" << endl;
            c << "//      in_ds      The input dataset." << endl;
            c << "//      <unused>   The domain number." << endl;
            c << "//      <unused>   The label." << endl;
            c << "//" << endl;
            c << "//  Returns:       The output dataset." << endl;
            c << "//" << endl;
            c << "//  Programmer: "<<getenv("USER")<<" -- generated by xml2info" << endl;
            c << "//  Creation:   "<<CurrentTime()<< endl;
            c << "//" << endl;
            c << "// ****************************************************************************" << endl;
            c << endl;
            c << "vtkDataSet *" << endl;
            c << "avt"<<name<<"Filter::ExecuteData(vtkDataSet *in_ds, int, std::string)" << endl;
            c << "{" << endl;
            c << "    YOUR CODE TO MODIFY THE DATASET GOES HERE" << endl;
            c << "}" << endl;
        }
        else if (type=="plot")
        {
            c << "// ************************************************************************* //" << endl;
            c << "//                              avt"<<name<<"Filter.C                              //" << endl;
            c << "// ************************************************************************* //" << endl;
            c << endl;
            c << "#include <avt"<<name<<"Filter.h>" << endl;
            c << endl;
            c << "#include <vtkDataSet.h>" << endl;
            c << endl;
            c << endl;
            c << "// ****************************************************************************" << endl;
            c << "//  Method: avt"<<name<<"Filter constructor" << endl;
            c << "//" << endl;
            c << "//  Programmer: "<<getenv("USER")<<" -- generated by xml2info" << endl;
            c << "//  Creation:   "<<CurrentTime()<< endl;
            c << "//" << endl;
            c << "// ****************************************************************************" << endl;
            c << endl;
            c << "avt"<<name<<"Filter::avt"<<name<<"Filter(YOUR INITIALIZERS)" << endl;
            c << "{" << endl;
            c << "}" << endl;
            c << endl;
            c << endl;
            c << "// ****************************************************************************" << endl;
            c << "//  Method: avt"<<name<<"Filter destructor" << endl;
            c << "//" << endl;
            c << "//  Programmer: "<<getenv("USER")<<" -- generated by xml2info" << endl;
            c << "//  Creation:   "<<CurrentTime()<< endl;
            c << "//" << endl;
            c << "// ****************************************************************************" << endl;
            c << endl;
            c << "avt"<<name<<"Filter::~avt"<<name<<"Filter()" << endl;
            c << "{" << endl;
            c << "}" << endl;
            c << endl;
            c << endl;
            c << "YOUR ROUTINES TO SET THE PARAMETERS OF THE FILTERS" << endl;
            c << endl;
            c << endl;
            c << "// ****************************************************************************" << endl;
            c << "//  Method: avt"<<name<<"Filter::ExecuteData" << endl;
            c << "//" << endl;
            c << "//  Purpose:" << endl;
            c << "//      Does the actual VTK code to modify the dataset." << endl;
            c << "//" << endl;
            c << "//  Arguments:" << endl;
            c << "//      inDS      The input dataset." << endl;
            c << "//      <unused>  The domain number." << endl;
            c << "//      <unused>  The label." << endl;
            c << "//" << endl;
            c << "//  Returns:      The output dataset." << endl;
            c << "//" << endl;
            c << "//  Programmer: "<<getenv("USER")<<" -- generated by xml2info" << endl;
            c << "//  Creation:   "<<CurrentTime()<< endl;
            c << "//" << endl;
            c << "// ****************************************************************************" << endl;
            c << endl;
            c << "vtkDataSet *" << endl;
            c << "avt"<<name<<"Filter::ExecuteData(vtkDataSet *inDS, int, std::string)" << endl;
            c << "{" << endl;
            c << "    THIS IS THE REAL VTK CODE" << endl;
            c << "}" << endl;
            c << endl;
            c << endl;
            c << "// ****************************************************************************" << endl;
            c << "//  Method: avt"<<name<<"Filter::RefashionDataObjectInfo" << endl;
            c << "//" << endl;
            c << "//  Purpose:" << endl;
            c << "//      Allows the filter to change its output's data object information, which" << endl;
            c << "//      is a description of the data object." << endl;
            c << "//" << endl;
            c << "//  Programmer: "<<getenv("USER")<<" -- generated by xml2info" << endl;
            c << "//  Creation:   "<<CurrentTime()<< endl;
            c << "//" << endl;
            c << "// ****************************************************************************" << endl;
            c << endl;
            c << "void" << endl;
            c << "avt"<<name<<"Filter::RefashionDataObjectInfo(void)" << endl;
            c << "{" << endl;
            c << "    IF YOU SEE FUNNY THINGS WITH EXTENTS, ETC, YOU CAN CHANGE THAT HERE." << endl;
            c << "}" << endl;
        }
    }
};


// ----------------------------------------------------------------------------
//                           Override default types
// ----------------------------------------------------------------------------
#define Plugin       AVTGeneratorPlugin

#endif
