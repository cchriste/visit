/**
 *
 * @file        avtVsFileFormat.h
 *
 * @brief       Base class for VizSchema visit plugins
 *
 * @version $Id: avtVsFileFormat.h 27 2008-03-26 22:04:41Z sveta $
 *
 * Copyright &copy; 2007, Tech-X Corporation
 * See LICENSE file for conditions of use.
 *
 */

// ************************************************************************* //
//                            avtVsFileFormat.h                                 //
// ************************************************************************* //

#ifndef VS_FILE_FORMAT_H
#define VS_FILE_FORMAT_H

#include <vector>
using std::vector;


//#include <VsH5Reader.h>
#include <avtSTMDFileFormat.h>
#include <hdf5.h>
#include <visit-hdf5.h>

// Forward references to minimize compilation
class vtkDataSet;
class vtkDataArray;
class avtDatabaseMetaData;
class avtMeshMetaData;
class VsVariableWithMesh;
class VsUnstructuredMesh;
class VsUniformMesh;
class VsStructuredMesh;
class VsRectilinearMesh;
class VsRegistry;
class VsH5Reader;

/**
 * avtSTMDFileFormat is a base class for multi-domain, single-time
 * VisIt plugins.  It provides mesh, scalar, and vector data.  It
 * assumes an ordering of component minor in the data file.  This is
 * reversed to get Fortran ordering.
 *
 * No ability to account for tensors at the present time.
 */
class avtVsFileFormat: public avtSTMDFileFormat {
#if HDF5_VERSION_GE(1, 8, 1)
  public:

  /**
   * Construct a file reader from a data file.
   *
   * @param dfnm the name of the data file
   * @param newStride The stride to use when loading data
   */
  avtVsFileFormat(const char* dfnm, std::vector<int> settings);

  /**
   * Destructor
   */
  virtual ~avtVsFileFormat();

  /**
   * Get plugin type
   *
   * @return plugin type name
   */
  virtual const char* GetType(void) {
    return "Vs";
  };

  /**
   * Get the data selections
   *
   */
  virtual void RegisterDataSelections( const vector<avtDataSelection_p> &sels,
                                       vector<bool> *selectionsApplied );

  /**
   * Process the data selections
   *
   */
  bool ProcessDataSelections(int *mins, int *maxs, int *strides);

  /**
   * Get a mesh by name
   *
   * @param domain the number of the domain in this mesh
   * @param meshname the name of the mesh
   *
   * @return a pointer to the mesh. Caller assumes ownership (must delete when done)
   */
  virtual vtkDataSet* GetMesh(int domain, const char* meshname);

  /**
   * get a scalar variable by name
   *
   * @param domain the number of the domain in this mesh
   * @param varname the name of the variable
   *
   * @return a pointer to the variable. Caller assumes ownership (must delete when done)
   */
  virtual vtkDataArray* GetVar(int domain, const char* varname);

  /**
   * Free up any resources created by this object.
   */
  virtual void FreeUpResources(void);

  /**
   * Called to alert the database reader that VisIt is about to ask for data
   * at this timestep.  Since this reader is single-time, this is 
   * only provided for reference and does nothing.
   */
  virtual void ActivateTimestep(void);
  
  virtual void UpdateCyclesAndTimes();
  
  protected:
  /** Populate the meta data */
  virtual void PopulateDatabaseMetaData(avtDatabaseMetaData* md);

  /** The file containing the data */
  std::string dataFileName;

  /** Pointer to the reader */
  VsH5Reader* reader;

  /** Ensure data has been read **/
  void LoadData();

  private:
  /**
   * A counter to track the number of avtVsFileFormat objects in existence
   */
  static int instanceCounter;

  /**
   * A user-specified setting for the stride to use when loading data.
   * Default is 1 on all axes.
   */
  std::vector<int> stride;

  VsRegistry* registry;

    vector<avtDataSelection_p> selList;
    vector<bool>              *selsApplied;

  /**
   * Set the axis labels for a mesh.
   *
   * @param mmd a pointer to the object that needs the axis labels.
   */
  void setAxisLabels(avtMeshMetaData* mmd);

  /**
   * Create various meshes.
   */
  vtkDataSet* getUniformMesh(VsUniformMesh*);
  vtkDataSet* getUnstructuredMesh(VsUnstructuredMesh*);
  vtkDataSet* getRectilinearMesh(VsRectilinearMesh*);
  vtkDataSet* getStructuredMesh(VsStructuredMesh*);
  vtkDataSet* getPointMesh(VsVariableWithMesh*);
  vtkDataSet* getSplitPointMesh(VsUnstructuredMesh*);
  vtkDataSet* getCurve(int domain, const std::string& name);

  /**
   * Each type of object is added to the database with a separate method
   * for neatness.
   */
  void RegisterMeshes(avtDatabaseMetaData* md);
  void RegisterMdMeshes(avtDatabaseMetaData* md);
  void RegisterVarsWithMesh(avtDatabaseMetaData* md);
  void RegisterVars(avtDatabaseMetaData* md);
  void RegisterMdVars(avtDatabaseMetaData* md);
  void RegisterExpressions(avtDatabaseMetaData* md);
#else
  avtVsFileFormat(const char* dfnm) : avtSTMDFileFormat(&dfnm, 1) {;};
#endif
};

#endif

