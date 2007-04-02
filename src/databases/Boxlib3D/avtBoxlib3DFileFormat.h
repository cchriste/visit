// ************************************************************************* //
//                         avtBoxlib3DFileFormat.h                           //
// ************************************************************************* //

#ifndef AVT_BOXLIB_3D_FILE_FORMAT_H
#define AVT_BOXLIB_3D_FILE_FORMAT_H

#include <avtSTMDFileFormat.h>

#include <vector>
#include <string>

class VisMF;

// ****************************************************************************
//  Class: avtBoxlib3DFileFormat
//
//  Purpose:
//      A file format reader for single-timestep, 3D AMR Boxlib files.
//
//  Programmer:  Akira Haddox
//  Creation:    July 25, 2003
//
//  Modifications:
//
//    Hank Childs, Sat Nov  8 18:38:08 PST 2003
//    Overhauled class.  Made a STMD (from MTMD) to clean up class and also
//    to overcome current (overcome-able) infrastructure limitations with
//    SILs changing in MTMD.  Also added true AMR capabilities.
//
//    Hank Childs, Sun Mar  6 16:21:15 PST 2005
//    Add support for GeoDyne material names.
//
//    Hank Childs, Thu Jun 23 14:46:22 PDT 2005
//    Broadcast Header from proc. 0
//
// ****************************************************************************

class avtBoxlib3DFileFormat : public avtSTMDFileFormat
{
  public:
                          avtBoxlib3DFileFormat(const char *);
    virtual              ~avtBoxlib3DFileFormat();
    
    virtual const char   *GetType(void) { return "Boxlib3D File Format"; };
    virtual bool          HasInvariantSIL(void) const { return false; };
    virtual bool          HasInvariantMetaData(void) const { return false; };
    
    virtual int           GetCycle(void) { return cycle; };
 
    virtual vtkDataSet   *GetMesh(int, const char *);
    virtual vtkDataArray *GetVar(int, const char *);
    virtual vtkDataArray *GetVectorVar(int, const char *);

    virtual void          PopulateDatabaseMetaData(avtDatabaseMetaData *);

    virtual void         *GetAuxiliaryData(const char *var, int,
                                           const char *type, void *args,
                                           DestructorFunction &);
    
    virtual void          FreeUpResources(void);
    virtual void          ActivateTimestep(void);

  protected:
    // This relative location of the multifab files.  It contains entries for
    // all levels. There is no predefined length.
    std::vector<std::string>                multifabFilenames;

    // For each level and variable, which multifab file the variable
    // is stored in, and at what index (component in the file it is stored as.)
    // Indexed by level, variable.
    std::vector<std::vector<int> >          fabfileIndex;
    std::vector<std::vector<int> >          componentIds;

    std::string                             rootPath;

    int                                     nLevels;
    std::vector<int>                        patchesPerLevel;
    // These entries are per patch.
    std::vector<double>                     xMin;
    std::vector<double>                     xMax;
    std::vector<double>                     yMin;
    std::vector<double>                     yMax;
    std::vector<double>                     zMin;
    std::vector<double>                     zMax;
    // These entries are per level.
    std::vector<double>                     deltaX;
    std::vector<double>                     deltaY;
    std::vector<double>                     deltaZ;
    // This entry is per level, but level 0 is omitted.
    std::vector<int>                        refinement_ratio;

    bool                                    haveReadTimeAndCycle;
    double                                  time;
    int                                     cycle;
    std::string                             timestepPath;
    bool                                    initializedReader;
    bool                                    vf_names_for_materials;
    
    // Scalar vars listed in header.
    int                                     nVars;
    std::vector<std::string>                varNames;
    std::vector<int>                        varCentering;
    std::vector<bool>                       varUsedElsewhere;

    // Vectors are combinations of scalar vars. The index to those
    // variables are stored in vectorComponents.
    int                                     nVectors;
    std::vector<std::string>                vectorNames;
    std::vector<std::vector<int> >          vectorComponents;
    std::vector<int>                        vectorCentering;

    // Variable readers.
    // Indexed by multifabFile
    std::vector<VisMF *>                    mfReaders;

    // The VisMF class is buggy when it comes to freeing up
    // memory.  The "clearlist" is used to sidestep those issues.
    std::vector<int>                        clearlist;

    // Problem range
    double                                  probLo[3];
    double                                  probHi[3];
    
    int                                     nMaterials;

    static const int                        dimension = 3;

    void                                   *GetMaterial(const char *, int,
                                                        const char *, 
                                                        DestructorFunction &);
    void                                   *GetSpatialIntervalTree(
                                                        DestructorFunction &);
    VisMF*                                  GetVisMF(int index);
    void                                    ReadHeader(void);
    void                                    InitializeReader(void);
    void                                    CalculateDomainNesting(void);

    vtkDataSet *CreateGrid(double lo[3], double hi[3], double delta[3]) const;

    int                 GetGlobalPatchNumber(int level, int patch) const;
    void                GetLevelAndLocalPatchNumber(int global_patch,
                                  int &level, int &local_patch) const;
    void                GetDimensions(int *, int level, int patch) const;
    void                GetDimensions(int *, double*, double*, double*) const;
};


#endif


