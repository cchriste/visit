/**
 * @file VsMesh.h
 *
 *  @class VsMesh
 *  @brief Superclass for all other types of meshes.
 *  
 *  Provides a common interface for all other types of meshes, 
 *  except VsMDMeshes.
 *
 *  Created on: Apr 28, 2010
 *      Author: mdurant
 */

#ifndef VSMESH_H_
#define VSMESH_H_

#include <string>
#include <vector>
#include <map>

#include "VsRegistryObject.h"

class VsMDMesh;
class VsH5Dataset;
class VsH5Group;
class VsH5Attribute;
class VsH5Object;

class VsMesh : public VsRegistryObject {
public:
  virtual ~VsMesh();
  
  void write();
  bool isFortranOrder();
  bool isCompMinor();
  bool isCompMajor();
  
  virtual bool isUniformMesh() { return false; }
  virtual bool isRectilinearMesh() { return false; }
  virtual bool isStructuredMesh() { return false; }
  virtual bool isUnstructuredMesh() { return false; }
  
  virtual std::string getShortName();
  virtual std::string getPath();
  
  std::string getFullName();
  void getStringAttribute(std::string attName, std::string* value);
  std::string getIndexOrder();
  virtual std::string getKind() = 0;
  size_t getNumSpatialDims();
  size_t getNumTopologicalDims();

  VsH5Attribute* getAttribute(std::string name);
  static VsMesh* buildObject(VsH5Dataset* dataset);
  static VsMesh* buildObject(VsH5Group* group);
  
  std::string getAxisLabel(unsigned int axis);
  
  virtual void getMeshDataDims(std::vector<int>& dims) = 0;
  
  void setMDMesh(VsMDMesh* mdMesh, int dNumber);
  int getDomainNumber();
  VsMDMesh* getMDMesh();

  virtual bool hasTransform();
  virtual std::string getTransformName();
  virtual std::string getTransformedMeshName();
  
protected:
  VsMesh(VsH5Object* object);
  virtual bool initialize() = 0;
  bool initializeRoot();
  
  /** The spatial dimensionality */
  size_t numSpatialDims;

  /** The topological dimensionality 
   * Note: will always be <= the spatial dimensionality
   * (i.e. can have lower-dimension object in a higher dimension
   * but not vice versa)
   */
  size_t numTopologicalDims;

  /** Index order (Fortran vs C style) */
  std::string indexOrder;

  VsH5Object* h5Object;
  
  /** Stuff for meshes that are subordinate blocks of MD Meshes**/
  int domainNumber;
  VsMDMesh* mdMesh;
};

#endif /* VSMESH_H_ */
