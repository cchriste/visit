// ************************************************************************* //
//                             avtMiliFileFormat.h                           //
// ************************************************************************* //

#ifndef AVT_MILI_FILE_FORMAT_H
#define AVT_MILI_FILE_FORMAT_H

#include <vector>
#include <string>
#include <fstream.h>

extern "C" {
#include <mili.h>
}

#include <avtMTMDFileFormat.h>
#include <avtTypes.h>

class avtMaterial;
class vtkDataArray;
class vtkUnstructuredGrid;

using std::vector;

// ****************************************************************************
//  Class: avtMiliFileFormat
//
//  Purpose:
//      A file format reader for Mili.
//
//  Notes:       Much of the code was taken from Doug Speck's GRIZ reader.
//      
//  Programmer:  Hank Childs
//  Creation:    April  11, 2003
//
//  Modifications:
//    Akira Haddox, Fri May 23 08:30:01 PDT 2003
//    Added in support for multiple meshes within a Mili database.
//    Changed into a MTMD file format.
//
// ****************************************************************************

class avtMiliFileFormat : public avtMTMDFileFormat
{
  public:
                          avtMiliFileFormat(const char *);
    virtual              ~avtMiliFileFormat();
    
    virtual const char   *GetType(void) { return "Mili File Format"; };
    
    virtual void          GetCycles(vector<int> &);
    virtual int           GetNTimesteps(void);
 
    virtual vtkDataSet   *GetMesh(int, int, const char *);
    virtual vtkDataArray *GetVar(int, int, const char *);
    virtual vtkDataArray *GetVectorVar(int, int, const char *);

    virtual void          PopulateDatabaseMetaData(avtDatabaseMetaData *);

    virtual void         *GetAuxiliaryData(const char *var, int, int,
                                           const char *type, void *args,
                                           DestructorFunction &);

  protected:
    char *famroot;
    char *fampath;
    
    vector<Famid>    dbid;
    int                   ntimesteps;
    int                   ndomains;
    int                   nmeshes;

    vector<bool>     validateVars;
    vector<bool>     readMesh;
    int                   dims;

    vector<vector<int> >                   nnodes;
    vector<vector<int> >                   ncells;
    vector<vector<vtkUnstructuredGrid *> > connectivity;

    vector<vector< std::string > >         sub_records;
    vector<vector< int > >                 sub_record_ids;

    vector<vector< std::string > >         element_group_name;
    vector<vector< int > >                 connectivity_offset;
    vector<vector< int > >                 group_mesh_associations;

    vector< std::string >         vars;
    vector< avtCentering >        centering;
    vector< vector< vector<bool> > >        vars_valid;
    vector< avtVarType >          vartype;
    vector< int >                 var_mesh_associations;

    vector<int>                            nmaterials;
    vector<vector< avtMaterial * > >       materials;

    void                  ReadMesh(int dom);
    void                  ValidateVariables(int dom);
    avtMaterial *         ConstructMaterials(vector< vector<int*> >&,
                                            vector< vector<int> > &);
    int                   GetVariableIndex(const char *);
    int                   GetVariableIndex(const char *, int mesh_id);
    void                  GetSizeInfoForGroup(const char *, int &, int &, int);

    void                  DecodeMultiMeshVarname(const std::string &, 
                                                 std::string &, int &);

    inline void           OpenDB(int dom);
};


#endif


