#ifndef PF3DFILEFORMAT_H
#define PF3DFILEFORMAT_H
#include <PDBReader.h>
#include <avtSTMDFileFormat.h>
#include <iostream.h>
#include <vectortypes.h>

class avtFileFormatInterface;
class vtkDataSet;
class vtkDataArray;

// ****************************************************************************
// Class: PF3DFileFormat
//
// Purpose:
//   Reads the PF3D file format.
//
// Notes:      
//
// Programmer: Brad Whitlock
// Creation:   Tue Sep 16 10:28:16 PDT 2003
//
// Modifications:
//   Brad Whitlock, Tue Jul 20 10:50:54 PDT 2004
//   Rewrote for new PF3D file format.
//
// ****************************************************************************

class PF3DFileFormat : public PDBReader, public avtSTMDFileFormat
{
public:
    static avtFileFormatInterface *CreateInterface(PDBFileObject *pdb,
         const char *const *filenames, int nList, int nBlocks);

    PF3DFileFormat(const char *filename);
    virtual ~PF3DFileFormat();

    // Methods overrides for an STMD file format.
    virtual const char   *GetType();
    virtual void          ActivateTimestep(void);
    virtual int           GetCycle(void);
    virtual double        GetTime(void);

    virtual void          PopulateDatabaseMetaData(avtDatabaseMetaData *);
    virtual void         *GetAuxiliaryData(const char *var, int dom,
                                           const char *type, void *,
                                           DestructorFunction &df);
    virtual void          FreeUpResources();

    virtual vtkDataSet   *GetMesh(int, const char *var);
    virtual vtkDataArray *GetVar(int, const char *var);

protected:
    virtual void CloseFile(int);

    virtual bool IdentifyFormat();
    static  bool StaticIdentifyFormat(PDBFileObject *);

    //
    // This class contains the structured information from the
    // master file.
    //
    class MasterInformation
    {
        struct History
        {
            double time;
            int ivzsave;
            float tnowps;
            long ncyc;
            float dt;
        };

    public:
        MasterInformation();
        virtual ~MasterInformation();
        bool Read(PDBFileObject *pdb);
        int  GetNDomains() const;

        //
        // Methods to get pointers to the fields of the variable
        // sized "struct".
        //
        double        Get_time() const;
        int           Get_ivzsave() const;
        float         Get_tnowps() const;
        long          Get_ncyc() const;
        float         Get_dt() const;
        const double *Get_rhomin_vz() const;
        const double *Get_rhomax_vz() const;
        const double *Get_e0min_vz() const;
        const double *Get_e0max_vz() const;
        const long *  Get_domloc() const;
        const double *Get_xyzloc() const;
        const char *  Get_visnams() const;
        const char *  Get_visname_for_domain(int dom, int comp) const;
        const double *Get_e1min_vz() const;
        const double *Get_e1max_vz() const;
        const double *Get_iawmin_vz() const;
        const double *Get_iawmax_vz() const;
        const double *Get_e2min_vz() const;
        const double *Get_e2max_vz() const;
        const double *Get_epwmin_vz() const;
        const double *Get_epwmax_vz() const;

        static const int visnam_size;

        ostream &operator << (ostream &os);
    private:
        int CalculateBufferSize();

        int  nDomains;
        void *buffer;
    };

    class BOF
    {
    public:
        BOF();
        BOF(const BOF &);
        virtual ~BOF();
        static void Destruct(void *);
    
        int    size[3];
        float *data;
    };

private:
    void Initialize();

    std::string FilenameForDomain(int realDomain);
    PDBFileObject *GetDomainFileObject(int realDomain);

    void ReadVariableInformation();
    void SetUpDomainConnectivity();
    int  GetNDomains() const;
    void GetExtents(int dom, float *extents);
    void GetLogicalExtents(int dom, int *globalorigin, int *size);
    void GetVarMinMaxArrays(const std::string &name,
                            const double **minvals,
                            const double **maxvals);
    int  GetVariableIndex(const std::string &name) const;
    BOF *GetBOF(int realDomain, const char *varName);
    int  GetRealDomainIndex(int dom) const;
    std::string GetBOFKey(int realDomain, const char *varName) const;
    bool ReadStringVector(const char *name, stringVector &output);


    bool initialized;

    // Values that are read in when the file is opened.
    int               cycle;
    int               nx;          // .5 Number of cells in X 
    int               ny;          // .5 Number of cells in Y
    int               nz;          // Number of cells -1 in Z
    double            lenx;        // Size of the whole problem in X.
    double            leny;        // Size of the whole problem in Y.
    double            lenz;        // Size of the whole problem in Z.
    MasterInformation master;      // Contains information about all of the
                                   // domains in the problem.
    stringVector      glob_nams;   // Contains variable names.
    stringVector      int_nams;    // Contains variable names in file.
    stringVector      glob_units;  // Contains variable units.
    intVector         apply_exp;   // Contains whether we should apply exp().
    std::string       databaseComment; // Contains the database information.
    bool              compression; // Determines whether data is compressed.
};

#endif
