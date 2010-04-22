#include <hdf5.h>
#include <visit-hdf5.h>
#if HDF5_VERSION_GE(1, 8, 1)
/**
 * @file  VsSchema.h
 *
 * @class VsSchema
 *
 * @brief Describes how instant datasets and meshes are found in the
 * file compliant with the schema.
 *
 * Copyright &copy; 2008 by Tech-X Corporation
 */

#ifndef VS_SCHEMA
#define VS_SCHEMA

#include <string>

struct VsSchema {

  // MD elements
  static std::string mdAtt;

  // Elements of schema
  static std::string typeAtt;
  static std::string kindAtt;
  static std::string meshAtt;
  static std::string centeringAtt; // This is deprecated
  static std::string cellOffsetAtt; // Instead of offsetAtt
  static std::string indexOrderAtt; //component major/minor, index C/Fortran; compMinorC is default
  static std::string numSpatialDimsAtt;
  static std::string spatialIndicesAtt;
  static std::string labelsAtt;
  static std::string varKey;
  static std::string vsVarsKey;
  static std::string varWithMeshKey;
  static std::string meshKey;
  static std::string zonalCenteringKey;// Node center is default
  static std::string structuredMeshKey;

  // Index ordering...
  static std::string compMajorCKey;
  static std::string compMinorCKey;
  static std::string compMajorFKey;
  static std::string compMinorFKey;

  struct Uniform {
    static std::string key;
    static std::string deprecated_key;
    static std::string lowerBounds;
    static std::string startCell;
    static std::string numCells;
    static std::string upperBounds;
  };

  struct Rectilinear {
    static std::string key;

    static std::string axis0Key;
    static std::string axis0DefaultName;
    
    static std::string axis1Key;
    static std::string axis1DefaultName;
    
    static std::string axis2Key;
    static std::string axis2DefaultName;
  };

  struct Unstructured {
    static std::string key;
    static std::string defaultPolygonsName; //polygons
    static std::string defaultPolyhedraName; //polyhedra
    static std::string defaultPointsName; //points
    static std::string defaultLinesName; //Lines
    static std::string defaultTrianglesName; //Triangles
    static std::string defaultQuadrilateralsName; //Quadrilaterals
    static std::string defaultTetrahedralsName; //Tetrahedral
    static std::string defaultPyramidsName; //Pyramids
    static std::string defaultPrismsName; //Prisms
    static std::string defaultHexahedralsName; //Hexahedrals

    static std::string vsPolygons; //polygons
    static std::string vsPolyhedra; //polyhedra
    static std::string vsPoints; //points
    static std::string vsLines; //Lines
    static std::string vsTriangles; //Triangles
    static std::string vsQuadrilaterals; //Quadrilaterals
    static std::string vsTetrahedrals; //Tetrahedral
    static std::string vsPyramids; //Pyramids
    static std::string vsPrisms; //Prisms
    static std::string vsHexahedrals; //Hexahedrals
    static std::string vsPoints0; //points
    static std::string vsPoints1; //points
    static std::string vsPoints2; //points
  };

};

#endif
#endif

