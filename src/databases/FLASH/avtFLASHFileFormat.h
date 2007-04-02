// ************************************************************************* //
//                            avtFLASHFileFormat.h                           //
// ************************************************************************* //

#ifndef AVT_FLASH_FILE_FORMAT_H
#define AVT_FLASH_FILE_FORMAT_H

#include <avtSTMDFileFormat.h>

#include <vector>

#include <hdf5.h>


// ****************************************************************************
//  Class: avtFLASHFileFormat
//
//  Purpose:
//      Reads in FLASH files as a plugin to VisIt.
//
//  Programmer: Jeremy Meredith
//  Creation:   August 23, 2005
//
// ****************************************************************************

class avtFLASHFileFormat : public avtSTMDFileFormat
{
  public:
                       avtFLASHFileFormat(const char *);
    virtual           ~avtFLASHFileFormat() {;};

    virtual bool           HasInvariantMetaData(void) const { return false; };
    virtual bool           HasInvariantSIL(void) const { return false; };

    virtual const char    *GetType(void)   { return "FLASH"; };
    virtual void           FreeUpResources(void); 

    virtual void          *GetAuxiliaryData(const char *var, int,
                                            const char *type, void *args,
                                            DestructorFunction &);

    virtual vtkDataSet    *GetMesh(int, const char *);
    virtual vtkDataArray  *GetVar(int, const char *);
    virtual vtkDataArray  *GetVectorVar(int, const char *);

  protected:
    virtual void           PopulateDatabaseMetaData(avtDatabaseMetaData *);

    void ReadAllMetaData();

    void BuildDomainNesting();

    void ReadSimulationParameters();
    void ReadBlockStructure();
    void ReadBlockExtents();
    void ReadRefinementLevels();
    void ReadUnknownNames();
    void DetermineGlobalLogicalExtentsForAllBlocks();

  protected:
    struct SimParams
    {
        int total_blocks;
        double time;
        double timestep;
        double redshift;
        int nsteps;
        int nxb;
        int nyb;
        int nzb;
    };

    struct Block
    {
        int ID;
        int level;
        int parentID;
        int childrenIDs[8];
        int neighborIDs[6];
        double minSpatialExtents[3];
        double maxSpatialExtents[3];
        int    minGlobalLogicalExtents[3];
        int    maxGlobalLogicalExtents[3];
    };

  protected:
    std::string              filename;
    int                      dimension;
    int                      numBlocks;
    int                      numLevels;
    hid_t                    fileId;
    SimParams                simParams;
    std::vector<Block>       blocks;
    int                      numChildrenPerBlock;
    int                      numNeighborsPerBlock;
    int                      block_ndims[3];
    int                      block_zdims[3];
    double                   minSpatialExtents[3];
    double                   maxSpatialExtents[3];
    std::vector<std::string> varNames;
};


#endif
