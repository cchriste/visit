// ************************************************************************* //
//                            avtEnzoFileFormat.h                           //
// ************************************************************************* //

#ifndef AVT_Enzo_FILE_FORMAT_H
#define AVT_Enzo_FILE_FORMAT_H

#include <avtSTMDFileFormat.h>

#include <vector>


// ****************************************************************************
//  Class: avtEnzoFileFormat
//
//  Purpose:
//      Reads in Enzo files as a plugin to VisIt.
//
//  Programmer: Jeremy Meredith
//  Creation:   December 3, 2004
//
//  Modifications:
//    Jeremy Meredith, Fri Feb 11 18:15:49 PST 2005
//    Added HDF5 support to the existing HDF4 support.
//
// ****************************************************************************

class avtEnzoFileFormat : public avtSTMDFileFormat
{
  public:
                       avtEnzoFileFormat(const char *);
    virtual           ~avtEnzoFileFormat();

    virtual int         GetCyle(void);

    virtual const char    *GetType(void)   { return "Enzo"; };
    virtual void           FreeUpResources(void); 

    virtual vtkDataSet    *GetMesh(int, const char *);
    virtual vtkDataArray  *GetVar(int, const char *);
    virtual vtkDataArray  *GetVectorVar(int, const char *);

    virtual void          *GetAuxiliaryData(const char *var, int,
                                            const char *type, void *args,
                                            DestructorFunction &);
    
  protected:
    enum FileType { ENZO_FT_UNKNOWN, ENZO_FT_HDF4, ENZO_FT_HDF5 };
    FileType fileType;

    // DATA MEMBERS
    struct Grid
    {
        int              ID;
        std::vector<int> childrenID;
        int              parentID;
        int              level;

        int              numberOfParticles;

        double           minSpatialExtents[3];
        double           maxSpatialExtents[3];

        int              zdims[3];
        int              ndims[3];
        int              minLogicalExtentsInParent[3];
        int              maxLogicalExtentsInParent[3];
        int              minLogicalExtentsGlobally[3];
        int              maxLogicalExtentsGlobally[3];
        double           refinementRatio[3];

      public:
        void PrintRecursive(std::vector<Grid> &grids, int level = 0);
        void Print();
        void DetermineExtentsInParent(std::vector<Grid> &grids);
        void DetermineExtentsGlobally(int numLevels,std::vector<Grid> &grids);
    };

    std::string fname_base;
    std::string fnameB;
    std::string fnameH;

    std::vector<Grid> grids;
    int numGrids;
    int numLevels;

    int curCycle;
    double curTime;

    std::vector<std::string> varNames;
    std::vector<std::string> particleVarNames;

    virtual void           PopulateDatabaseMetaData(avtDatabaseMetaData *);

    void ReadAllMetaData();
    void ReadHierachyFile();
    void ReadParameterFile();
    void DetermineVariablesFromGridFile();
    void BuildDomainNesting();
};


#endif
