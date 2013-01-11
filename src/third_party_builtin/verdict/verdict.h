
/*! \file verdict.h
  \brief Header file for verdict library that calculates metrics for finite elements.
    Also see: \ref index "Main Page" 
 *
 * verdict.h is the header file for applications/libraries to include
 *           to compute quality metrics.
 *
 * Copyright (C) 2003 Sandia National Laboratories <cubit@sandia.gov>
 *
 * This file is part of VERDICT
 *
 * This copy of VERDICT is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * VERDICT is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 * 
*/


#ifndef VERDICT_INC_LIB
#define VERDICT_INC_LIB

#define HAVE_VERDICT

#define VERDICT_VERSION 110

#define VERDICT_PI 3.1415926535897932384626

/* note:  the VERDICT_USE_FLOAT must be consistent with the build of the library */

#ifdef VERDICT_USE_FLOAT
   #define VERDICT_REAL float
   #define VERDICT_DBL_MIN 1.0E-30
   #define VERDICT_DBL_MAX 1.0E+30
#else
   #define VERDICT_REAL double
   #define VERDICT_DBL_MIN 1.0E-300
   #define VERDICT_DBL_MAX 1.0E+300
#endif

#if defined(_WIN32)
  #if defined(VERDICT_EXPORTS) || defined(visit_verdict_EXPORTS)
    #define C_FUNC_DEF __declspec(dllexport)
  #else
    #define C_FUNC_DEF __declspec(dllimport)
  #endif
#else
  #ifdef __cplusplus
    #if __GNUC__ >= 4 && defined(VERDICT_EXPORTS)
    # define C_FUNC_DEF extern "C" __attribute__ ((visibility("default")))
    #else
    # define C_FUNC_DEF extern "C"
    #endif
  #else
    #define C_FUNC_DEF
  #endif
#endif


/* typedef for the user if they want to use
 *    function pointers */
#ifdef __cplusplus
extern "C" {
#endif
  typedef VERDICT_REAL(*VerdictFunction)(int, VERDICT_REAL[][3]);
  typedef int(*ComputeNormal)(double point[3], double normal[3]); 
#ifdef __cplusplus
}
#endif


/** HexMetricVals is a <em>struct</em> used to return calculated metrics  
    when calling the function <em>v_hex_quality(...)</em> 

    HexMetricVals is used just as QuadMetricVals.
    For an example, see Detailed Description in QuadMetricVals  
*/
struct HexMetricVals
{
  /** \sa v_hex_aspect */
  VERDICT_REAL aspect ;  
  /** \sa v_hex_skew */
  VERDICT_REAL skew ;
  /** \sa v_hex_taper */
  VERDICT_REAL taper ;
  /** \sa v_hex_volume */
  VERDICT_REAL volume ;
  /** \sa v_hex_stretch */
  VERDICT_REAL stretch ;
  /** \sa v_hex_max_diagonal */
  VERDICT_REAL max_diagonal ;
  /** \sa v_hex_min_diagonal */
  VERDICT_REAL min_diagonal ;
  /** \sa v_hex_diagonal_ratio */
  VERDICT_REAL diagonal_ratio ;
  /** \sa v_hex_dimension */
  VERDICT_REAL dimension ;
  /** \sa v_hex_oddy */
  VERDICT_REAL oddy ;
  /** \sa v_hex_condition */
  VERDICT_REAL condition ;
  /** \sa v_hex_jacobian */
  VERDICT_REAL jacobian ;
  /** \sa v_hex_scaled_jacobian */
  VERDICT_REAL scaled_jacobian ;
  /** \sa v_hex_shear */
  VERDICT_REAL shear ;
  /** \sa v_hex_shape */
  VERDICT_REAL shape ;
  /** \sa v_hex_relative_size */
  VERDICT_REAL relative_size_squared;
  /** \sa v_hex_shape_and_size */
  VERDICT_REAL shape_and_size ; 
  /** \sa v_hex_shear_and_size */
  VERDICT_REAL shear_and_size ; 
  /** \sa v_hex_distortion */
  VERDICT_REAL distortion; 
};

/** EdgeMetricVals is a <em>struct</em> used to return calculated metrics  
    when calling the function <em>v_edge_quality(...)</em> 

    EdgeMetricVals is used just as QuadMetricVals.
    For an example, see Detailed Description in QuadMetricVals  
*/
struct EdgeMetricVals
{
  VERDICT_REAL length ; 
};


/** KnifeMetricVals is a <em>struct</em> used to return calculated metrics  
    when calling the function <em>v_knife_quality(...)</em> 

    KnifeMetricVals is used just as QuadMetricVals.
    For an example, see Detailed Description in QuadMetricVals  
*/
struct KnifeMetricVals
{
  VERDICT_REAL volume ; 
};


/** QuadMetricVals is a <em>struct</em> used to return calculated metrics  
    when calling the function <em>v_quad_quality(...)</em> 
 
    The following is an example of how this struct is used with Verdict.

    Example: \code
    QuadMetricVals quad_metrics = {0};
    unsigned long metrics_flag = 0;
    metrics_flag += V_QUAD_SHAPE;
    metrics_flag += V_QUAD_DISTORTION;
    metrics_flag += V_QUAD_AREA;   
    double quad_nodes[4][3];
    get_quad_nodes( quad_nodes );  //some user-defined function to load 
                                   //xyz coordinate info. into array  
    v_quad_quality( 4, quad_nodes, metrics_flag, quad_metrics );
    double my_shape      = quad_metrics.shape; 
    double my_distortion = quad_metrics.distortion; 
    double my_area       = quad_metrics.area;  \endcode
     
 */

struct QuadMetricVals
{
  /** \sa v_quad_aspect function */
  VERDICT_REAL aspect ;
  /** \sa v_quad_skew function */
  VERDICT_REAL skew ;
  /** \sa v_quad_taper function */
  VERDICT_REAL taper ;
  /** \sa v_quad_warpage function */
  VERDICT_REAL warpage ;
  /** \sa v_quad_area function */
  VERDICT_REAL area ;
  /** \sa v_quad_stretch function */
  VERDICT_REAL stretch ;
  /** \sa v_quad_smallest_angle function */
  VERDICT_REAL minimum_angle ;
  /** \sa v_quad_largest_angle function */
  VERDICT_REAL maximum_angle ;
  /** \sa v_quad_oddy function */
  VERDICT_REAL oddy ;
  /** \sa v_quad_condition function */
  VERDICT_REAL condition ;
  /** \sa v_quad_jacobian function */
  VERDICT_REAL jacobian ;
  /** \sa v_quad_scaled_jacobian function */
  VERDICT_REAL scaled_jacobian ;
  /** \sa v_quad_shear function */
  VERDICT_REAL shear ;
  /** \sa v_quad_shape function */
  VERDICT_REAL shape ;
  /** \sa v_quad_relative_size_squared function */
  VERDICT_REAL relative_size_squared ;
  /** \sa v_quad_shape_and_size function */
  VERDICT_REAL shape_and_size ; 
  /** \sa v_quad_shear_and_size function */
  VERDICT_REAL shear_and_size ;
  /** \sa v_quad_distortion function */
  VERDICT_REAL distortion; 
};

/** PyramidMetricVals is a <em>struct</em> used to return calculated metrics  
    when calling the function <em>v_pyramid_quality(...)</em> 

    PyramidMetricVals is used just as QuadMetricVals.
    For an example, see Detailed Description in QuadMetricVals  
*/
struct PyramidMetricVals
{
  VERDICT_REAL volume ; 
};
   
/** WedgeMetricVals is a <em>struct</em> used to return calculated metrics  
    when calling the function <em>v_wedge_quality(...)</em> 

    WedgeMetricVals is used just as QuadMetricVals.
    For an example, see Detailed Description in QuadMetricVals  
*/
struct WedgeMetricVals
{
  VERDICT_REAL volume ; 
};
     
/** TetMetricVals is a <em>struct</em> used to return calculated metrics  
    when calling the function <em>v_tet_quality(...)</em> 

    TetMetricVals is used just as QuadMetricVals.
    For an example, see Detailed Description in QuadMetricVals  
*/
struct TetMetricVals
{
  /** \sa v_tet_aspect*/
  VERDICT_REAL aspect_beta;
  /** \sa v_hex_aspect_gamma */
  VERDICT_REAL aspect_gamma ;
  /** \sa v_tet_volume */
  VERDICT_REAL volume ;
  /** \sa v_tet_condition */
  VERDICT_REAL condition ;
  /** \sa v_tet_jacobian */
  VERDICT_REAL jacobian ;
  /** \sa v_tet_scaled_jacobian */
  VERDICT_REAL scaled_jacobian ;
  /** \sa v_tet_shape */
  VERDICT_REAL shape ;
  /** \sa v_tet_relative_size */
  VERDICT_REAL relative_size_squared ;
  /** \sa v_tet_shape_and_size*/
  VERDICT_REAL shape_and_size ; 
  /** \sa v_tet_distortion */
  VERDICT_REAL distortion; 
};

/** TriMetricVals is a <em>struct</em> used to return calculated metrics  
    when calling the function <em>v_tri_quality(...)</em> 

    TriMetricVals is used just as QuadMetricVals.
    For an example, see Detailed Description in QuadMetricVals  
*/
struct TriMetricVals
{
  /** \sa v_tri_aspect */
  VERDICT_REAL aspect ;
  /** \sa v_tri_area*/
  VERDICT_REAL area ;
  /** \sa v_tri_minimum_angle*/
  VERDICT_REAL minimum_angle ;
  /** \sa v_tri_maximum_angle */
  VERDICT_REAL maximum_angle ;
  /** \sa v_tri_condition */
  VERDICT_REAL condition ;
  /** \sa v_tri_scaled_jacobian */
  VERDICT_REAL scaled_jacobian ;
  /** \sa v_tri_shape */
  VERDICT_REAL shape ;
  /** \sa v_tri_relative_size_squared */
  VERDICT_REAL relative_size_squared ;
  /** \sa v_tri_shape_and_size */
  VERDICT_REAL shape_and_size ; 
  /** \sa v_tri_distortion */
  VERDICT_REAL distortion; 
};



/* definition of bit fields to determine which metrics to calculate */

//! \name Hex bit fields
//! 
//@{

#define V_HEX_ASPECT                 0x000001    /*!< \hideinitializer */
#define V_HEX_SKEW                   0x000002    /*!< \hideinitializer */
#define V_HEX_TAPER                  0x000004    /*!< \hideinitializer */
#define V_HEX_VOLUME                 0x000008    /*!< \hideinitializer */
#define V_HEX_STRETCH                0x000010    /*!< \hideinitializer */
#define V_HEX_DIAGONAL_RATIO         0x000020    /*!< \hideinitializer */
#define V_HEX_DIMENSION              0x000040    /*!< \hideinitializer */
#define V_HEX_ODDY                   0x000080    /*!< \hideinitializer */
#define V_HEX_CONDITION              0x000100    /*!< \hideinitializer */
#define V_HEX_JACOBIAN               0x000200    /*!< \hideinitializer */
#define V_HEX_SCALED_JACOBIAN        0x000400    /*!< \hideinitializer */
#define V_HEX_SHEAR                  0x000800    /*!< \hideinitializer */
#define V_HEX_SHAPE                  0x001000    /*!< \hideinitializer */
#define V_HEX_RELATIVE_SIZE_SQUARED  0x002000    /*!< \hideinitializer */
#define V_HEX_SHAPE_AND_SIZE         0x004000    /*!< \hideinitializer */
#define V_HEX_SHEAR_AND_SIZE         0x008000    /*!< \hideinitializer */
#define V_HEX_DISTORTION             0x010000    /*!< \hideinitializer */
#define V_HEX_MIN_DIAGONAL           0x020000    /*!< \hideinitializer */
#define V_HEX_MAX_DIAGONAL           0x040000    /*!< \hideinitializer */
#define V_HEX_ALL                    0x07ffff    /*!< \hideinitializer */
/*!< \hideinitializer */
#define V_HEX_TRADITIONAL            V_HEX_ASPECT          + \
                                     V_HEX_SKEW            + \
                                     V_HEX_TAPER           + \
                                     V_HEX_STRETCH         + \
                                     V_HEX_DIAGONAL_RATIO  + \
                                     V_HEX_ODDY            + \
                                     V_HEX_CONDITION       + \
                                     V_HEX_JACOBIAN        + \
                                     V_HEX_SCALED_JACOBIAN + \
                                     V_HEX_DIMENSION

/*!< \hideinitializer */
#define V_HEX_DIAGNOSTIC             V_HEX_VOLUME
/*!< \hideinitializer */
#define V_HEX_ALGEBRAIC              V_HEX_SHAPE                  + \
                                     V_HEX_SHEAR                  + \
                                     V_HEX_RELATIVE_SIZE_SQUARED  + \
                                     V_HEX_SHAPE_AND_SIZE         + \
                                     V_HEX_SHEAR_AND_SIZE
/*!< \hideinitializer */
#define V_HEX_ROBINSON               V_HEX_SKEW + \
                                     V_HEX_TAPER    
//@}
                                     
//! \name Tet bit fields
//! 
//@{
#define V_TET_ASPECT_BETA            1   /*!< \hideinitializer */
#define V_TET_ASPECT_GAMMA           2   /*!< \hideinitializer */
#define V_TET_VOLUME                 4   /*!< \hideinitializer */
#define V_TET_CONDITION              8   /*!< \hideinitializer */
#define V_TET_JACOBIAN               16   /*!< \hideinitializer */
#define V_TET_SCALED_JACOBIAN        32   /*!< \hideinitializer */
#define V_TET_SHAPE                  64   /*!< \hideinitializer */
#define V_TET_RELATIVE_SIZE_SQUARED  128   /*!< \hideinitializer */
#define V_TET_SHAPE_AND_SIZE         256   /*!< \hideinitializer */
#define V_TET_DISTORTION             512   /*!< \hideinitializer */
#define V_TET_ALL                    1023   /*!< \hideinitializer */
/*!< \hideinitializer */
#define V_TET_TRADITIONAL            V_TET_ASPECT_BETA + \
                                     V_TET_ASPECT_GAMMA + \
                                     V_TET_CONDITION + \
                                     V_TET_JACOBIAN + \
                                     V_TET_SCALED_JACOBIAN  
/*!< \hideinitializer */
#define V_TET_DIAGNOSTIC             V_TET_VOLUME
/*!< \hideinitializer */
#define V_TET_ALGEBRAIC              V_TET_SHAPE                  + \
                                     V_TET_RELATIVE_SIZE_SQUARED  + \
                                     V_TET_SHAPE_AND_SIZE
//@}

 
//! \name Pyramid bit fields
//! 
//@{
#define V_PYRAMID_VOLUME             1   /*!< \hideinitializer */
//@}

//! \name Wedge bit fields
//! 
//@{
#define V_WEDGE_VOLUME               1   /*!< \hideinitializer */
//@}
 
//! \name Knife bit fields
//! 
//@{
#define V_KNIFE_VOLUME               1   /*!< \hideinitializer */
//@}
 
//! \name Quad bit fields
//! 
//@{
#define V_QUAD_ASPECT                1   /*!< \hideinitializer */
#define V_QUAD_SKEW                  2   /*!< \hideinitializer */
#define V_QUAD_TAPER                 4   /*!< \hideinitializer */
#define V_QUAD_WARPAGE               8   /*!< \hideinitializer */
#define V_QUAD_AREA                  16   /*!< \hideinitializer */
#define V_QUAD_STRETCH               32   /*!< \hideinitializer */
#define V_QUAD_MINIMUM_ANGLE         64   /*!< \hideinitializer */
#define V_QUAD_MAXIMUM_ANGLE         128   /*!< \hideinitializer */
#define V_QUAD_ODDY                  256   /*!< \hideinitializer */
#define V_QUAD_CONDITION             512   /*!< \hideinitializer */
#define V_QUAD_JACOBIAN              1024   /*!< \hideinitializer */
#define V_QUAD_SCALED_JACOBIAN       2048   /*!< \hideinitializer */
#define V_QUAD_SHEAR                 4096   /*!< \hideinitializer */
#define V_QUAD_SHAPE                 8192   /*!< \hideinitializer */
#define V_QUAD_RELATIVE_SIZE_SQUARED 16384   /*!< \hideinitializer */
#define V_QUAD_SHAPE_AND_SIZE        32768   /*!< \hideinitializer */
#define V_QUAD_SHEAR_AND_SIZE        65536   /*!< \hideinitializer */
#define V_QUAD_DISTORTION            131072   /*!< \hideinitializer */
#define V_QUAD_ALL                   262143   /*!< \hideinitializer */
/*!< \hideinitializer */
#define V_QUAD_TRADITIONAL           V_QUAD_ASPECT         + \
                                     V_QUAD_SKEW           + \
                                     V_QUAD_TAPER          + \
                                     V_QUAD_WARPAGE        + \
                                     V_QUAD_STRETCH        + \
                                     V_QUAD_MINIMUM_ANGLE + \
                                     V_QUAD_MAXIMUM_ANGLE  + \
                                     V_QUAD_ODDY           + \
                                     V_QUAD_CONDITION      + \
                                     V_QUAD_JACOBIAN       + \
                                     V_QUAD_SCALED_JACOBIAN 
/*!< \hideinitializer */
#define V_QUAD_DIAGNOSTIC            V_QUAD_AREA
/*!< \hideinitializer */
#define V_QUAD_ALGEBRAIC             V_QUAD_SHEAR                 + \
                                     V_QUAD_SHAPE                 + \
                                     V_QUAD_RELATIVE_SIZE_SQUARED + \
                                     V_QUAD_SHAPE_AND_SIZE     
/*!< \hideinitializer */
#define V_QUAD_ROBINSON              V_QUAD_ASPECT + \
                                     V_QUAD_SKEW   + \
                                     V_QUAD_TAPER
//@}


//! \name Tri bit fields
//! 
//@{
#define V_TRI_ASPECT                 1   /*!< \hideinitializer */
#define V_TRI_AREA                   2   /*!< \hideinitializer */
#define V_TRI_MINIMUM_ANGLE          4   /*!< \hideinitializer */
#define V_TRI_MAXIMUM_ANGLE          8   /*!< \hideinitializer */
#define V_TRI_CONDITION              16   /*!< \hideinitializer */
#define V_TRI_SCALED_JACOBIAN        32   /*!< \hideinitializer */
#define V_TRI_SHAPE                  64   /*!< \hideinitializer */
#define V_TRI_RELATIVE_SIZE_SQUARED  128   /*!< \hideinitializer */
#define V_TRI_SHAPE_AND_SIZE         256   /*!< \hideinitializer */
#define V_TRI_DISTORTION             512   /*!< \hideinitializer */
#define V_TRI_ALL                    1023   /*!< \hideinitializer */
/*!< \hideinitializer */
#define V_TRI_TRADITIONAL            V_TRI_ASPECT + \
                                     V_TRI_MINIMUM_ANGLE + \
                                     V_TRI_MAXIMUM_ANGLE + \
                                     V_TRI_CONDITION + \
                                     V_TRI_SCALED_JACOBIAN 
/*!< \hideinitializer */
#define V_TRI_DIAGNOSTIC             V_TRI_AREA
/*!< \hideinitializer */
#define V_TRI_ALGEBRAIC              V_TRI_SHAPE + \
                                     V_TRI_SHAPE_AND_SIZE + \
                                     V_TRI_RELATIVE_SIZE_SQUARED
 
#define V_EDGE_LENGTH                1   /*!< \hideinitializer */
//@}

                                     
/*! \mainpage
  Verdict is a library used to calculate metrics on the following type of elements:

    \li Hexahedra
    \li Tetrahedra
    \li Pryamid 
    \li Wedge 
    \li Knife 
    \li Quadrilateral
    \li Triangle
    \li Edge 

  Verdict calculates individual or multiple metrics on a single elment.  
  The v_*_quality(...) functions allow for efficient calculations of 
  multiple metrics on a single element.  Individual metrics may be 
  calculated on a single element as well. 

  \section jack Using Verdict

  The v_*_quality functions take the following parameters: 

  \param num_nodes Number of nodes in the element. 
  \param coordinates 2D array containing x,y,z coordinate data of the nodes.
  \param metrics_request_flag Bitfield to define which metrics to calculate.
  \param metric_vals Struct to hold the metric calculated values. 

  All other functions take these parameters below and return the calculated
  metric value.

  \param num_nodes Number of nodes in the element. 
  \param coordinates 2D array containing x,y,z coordinate data of the nodes.

  
  \par Setting the metrics_request_flag: 
  In order to use v_*_quality functions you must know how to set the bitfield argument 
  correctly.  To calculate aspect ratio, condition number, shape and shear metrics on a triangle, set
  the "metrics_request_flag" like so: 

  \code
  unsigned int metrics_request_flag = 0;
  metrics_request_flag += V_TRI_ASPECT;  
  metrics_request_flag += V_CONDITION;
  metrics_request_flag += V_SHAPE;
  metrics_request_flag += V_SHEAR;
  \endcode 

  The bitwise field can also be set for many metrics at once using #deinfed numbers.  V_HEX_ALL,
  V_HEX_DIAGNOSTIC, V_TRI_ALGEBRAIC are examples.

  Below is an example of how use Verdict's functions:

    
    Example: \code
    QuadMetricVals quad_metrics = {0};
    unsigned long metrics_flag = 0;
    metrics_flag += V_QUAD_SHAPE;
    metrics_flag += V_QUAD_DISTORTION;
    metrics_flag += V_QUAD_AREA;
    double quad_nodes[4][3];
 
    //node 1
    quad_node[0][0] = 0;  //x
    quad_node[0][1] = 0;  //y
    quad_node[0][2] = 0;  //z

    //node 2
    quad_node[1][0] = 1; 
    quad_node[1][1] = 0.1;
    quad_node[1][2] = 0.1;

    //node 3
    quad_node[2][0] = 0.9; 
    quad_node[2][1] = 0.9;
    quad_node[2][2] = -0.1;

    //node 4
    quad_node[3][0] = -0.05; 
    quad_node[3][1] = 1;
    quad_node[3][2] = 0;

    //calculate multiple metrics with one call
    v_quad_quality( 4, quad_nodes, metrics_flag, quad_metrics );
    double my_shape      = quad_metrics.shape; 
    double my_distortion = quad_metrics.distortion; 
    double my_area       = quad_metrics.area;  
  
    //calculate an individual metric 
    double my_relative_size = v_quad_relative_size( 4, quad_nodes );
    \endcode

*/

    //! Calculates quality metrics for hexahedral elements.
    C_FUNC_DEF void v_hex_quality( int num_nodes, VERDICT_REAL coordinates[][3], 
        unsigned int metrics_request_flag, struct HexMetricVals *metric_vals ); 
    
    //! Calculates quality metrics for tetrahedral elements.
    C_FUNC_DEF void v_tet_quality( int num_nodes, VERDICT_REAL coordinates[][3], 
        unsigned int metrics_request_flag, struct TetMetricVals *metric_vals ); 
    
    //! Calculates quality metrics for pyramid elements.
    C_FUNC_DEF void v_pyramid_quality( int num_nodes, VERDICT_REAL coordinates[][3], 
        unsigned int metrics_request_flag, struct PyramidMetricVals *metric_vals ); 

    //! Calculates quality metrics for wedge elements.
    C_FUNC_DEF void v_wedge_quality( int num_nodes, VERDICT_REAL coordinates[][3], 
        unsigned int metrics_request_flag, struct WedgeMetricVals *metric_vals ); 

    //! Calculates quality metrics for knife elements.
    C_FUNC_DEF void v_knife_quality( int num_nodes, VERDICT_REAL coordinates[][3], 
        unsigned int metrics_request_flag, struct KnifeMetricVals *metric_vals ); 

    //! Calculates quality metrics for quadralateral elements.
    C_FUNC_DEF void v_quad_quality( int num_nodes, VERDICT_REAL coordinates[][3], 
        unsigned int metrics_request_flag, struct QuadMetricVals *metric_vals ); 

    //! Calculates quality metrics for triangle elements.
    C_FUNC_DEF void v_tri_quality( int num_nodes, VERDICT_REAL coordinates[][3], 
        unsigned int metrics_request_flag, struct TriMetricVals *metric_vals );

    //! Calculates quality metrics for edge elements.
    C_FUNC_DEF void v_edge_quality( int num_nodes, VERDICT_REAL coordinates[][3], 
        unsigned int metrics_request_flag, struct EdgeMetricVals *metric_vals ); 



/* individual quality functions for hex elements */

    //! Sets average size (volume) of hex, needed for v_hex_relative_size(...)
    C_FUNC_DEF void v_set_hex_size( VERDICT_REAL size );

    //! Calculates hex aspect ratio 
    /**Maximum edge length ratios at hex center.
      Reference --- L.M. Taylor, and D.P. Flanagan, Pronto3D - A Three Dimensional Transient
         Solid Dynamics Program, SAND87-1912, Sandia National Laboratories, 1989. */
    C_FUNC_DEF VERDICT_REAL v_hex_aspect( int num_nodes, VERDICT_REAL coordinates[][3] );

    //! Calculates hex skew metric. 
    /** Maximum |cos A| where A is the angle between edges at hex center.   
      Reference --- L.M. Taylor, and D.P. Flanagan, Pronto3D - A Three Dimensional Transient
         Solid Dynamics Program, SAND87-1912, Sandia National Laboratories, 1989. */
    C_FUNC_DEF VERDICT_REAL v_hex_skew( int num_nodes, VERDICT_REAL coordinates[][3] ); 

    //! Calculates hex taper metric 
    /**  Maximum ratio of lengths derived from opposite edges. 
      Reference --- L.M. Taylor, and D.P. Flanagan, Pronto3D - A Three Dimensional Transient
         Solid Dynamics Program, SAND87-1912, Sandia National Laboratories, 1989. */
    C_FUNC_DEF VERDICT_REAL v_hex_taper( int num_nodes, VERDICT_REAL coordinates[][3] ); 

    //! Calculates hex volume 
    /**  Jacobian at hex center. 
      Reference --- L.M. Taylor, and D.P. Flanagan, Pronto3D - A Three Dimensional Transient
         Solid Dynamics Program, SAND87-1912, Sandia National Laboratories, 1989. */
    C_FUNC_DEF VERDICT_REAL v_hex_volume( int num_nodes, VERDICT_REAL coordinates[][3] ); 

    //! Calculates hex stretch metric   
    /**  Sqrt(3) * minimum edge length / maximum diagonal length. 
      Reference --- FIMESH code */ 
    C_FUNC_DEF VERDICT_REAL v_hex_stretch( int num_nodes, VERDICT_REAL coordinates[][3] ); 

    //! Calculates hex diagonal metric   
    /** Minimum diagonal length / maximum diagonal length. 
      Reference --- Unknown */ 
    C_FUNC_DEF VERDICT_REAL v_hex_diagonal_ratio( int num_nodes, VERDICT_REAL coordinates[][3] ); 

    //! Calculates max hex diagonal metric   
    C_FUNC_DEF VERDICT_REAL v_hex_max_diagonal( int num_nodes, VERDICT_REAL coordinates[][3] ); 

    //! Calculates min hex diagonal metric   
    C_FUNC_DEF VERDICT_REAL v_hex_min_diagonal( int num_nodes, VERDICT_REAL coordinates[][3] ); 

    //! Calculates hex dimension metric   
    /** Pronto-specific characteristic length for stable time step calculation.  
        Char_length = Volume / 2 grad Volume. 
      Reference --- L.M. Taylor, and D.P. Flanagan, Pronto3D - A Three Dimensional Transient
         Solid Dynamics Program, SAND87-1912, Sandia National Laboratories, 1989. */
    C_FUNC_DEF VERDICT_REAL v_hex_dimension( int num_nodes, VERDICT_REAL coordinates[][3] ); 

    //! Calculates hex oddy metric   
    C_FUNC_DEF VERDICT_REAL v_hex_oddy( int num_nodes, VERDICT_REAL coordinates[][3] ); 

    //! Calculates hex condition metric   
    /** Maximum condition number of the Jacobian matrix at 8 corners.
       Reference --- P. Knupp, Achieving Finite Element Mesh Quality via 
       Optimization of the Jacobian Matrix Norm and Associated Quantities, 
       Intl. J. Numer. Meth. Engng. 2000, 48:1165-1185. */ 
    C_FUNC_DEF VERDICT_REAL v_hex_condition( int num_nodes, VERDICT_REAL coordinates[][3] ); 

    //! Calculates hex jacobian metric   
    /** Minimum pointwise volume of local map at 8 corners & center of hex. 
       Reference --- P. Knupp, Achieving Finite Element Mesh Quality via 
       Optimization of the Jacobian Matrix Norm and Associated Quantities, 
       Intl. J. Numer. Meth. Engng. 2000, 48:1165-1185. */ 
    C_FUNC_DEF VERDICT_REAL v_hex_jacobian( int num_nodes, VERDICT_REAL coordinates[][3] ); 
    
    //! Calculates hex scaled jacobian metric   
    /** Minimum Jacobian divided by the lengths of the 3 edge vectors. 
       Reference --- P. Knupp, Achieving Finite Element Mesh Quality via 
       Optimization of the Jacobian Matrix Norm and Associated Quantities, 
       Intl. J. Numer. Meth. Engng. 2000, 48:1165-1185. */ 
    C_FUNC_DEF VERDICT_REAL v_hex_scaled_jacobian( int num_nodes, VERDICT_REAL coordinates[][3] ); 

    //! Calculates hex shear metric   
    /** 3/Mean Ratio of Jacobian Skew matrix.
       Reference --- P. Knupp, Algebraic Mesh Quality Metrics for
       Unstructured Initial Meshes, submitted for publication.  */
    C_FUNC_DEF VERDICT_REAL v_hex_shear( int num_nodes, VERDICT_REAL coordinates[][3] ); 

    //! Calculates hex shape metric.
    /** 3/Mean Ratio of weighted Jacobian matrix. 
       Reference --- P. Knupp, Algebraic Mesh Quality Metrics for
       Unstructured Initial Meshes, submitted for publication.  */
    C_FUNC_DEF VERDICT_REAL v_hex_shape( int num_nodes, VERDICT_REAL coordinates[][3] ); 

    //! Calculates hex relative size metric. 
    /** 3/Mean Ratio of weighted Jacobian matrix.
       Reference --- P. Knupp, Algebraic Mesh Quality Metrics for
       Unstructured Initial Meshes, submitted for publication.  */
    C_FUNC_DEF VERDICT_REAL v_hex_relative_size_squared( int num_nodes, VERDICT_REAL coordinates[][3] ); 

    //! Calculates hex shape-size metric.
    /** Product of Shape and Relative Size.
       Reference --- P. Knupp, Algebraic Mesh Quality Metrics for
       Unstructured Initial Meshes, submitted for publication.  */
    C_FUNC_DEF VERDICT_REAL v_hex_shape_and_size( int num_nodes, VERDICT_REAL coordinates[][3] ); 

    //! Calculates hex shear-size metric   
    /** Product of Shear and Relative Size.
       Reference --- P. Knupp, Algebraic Mesh Quality Metrics for
       Unstructured Initial Meshes, submitted for publication.  */
    C_FUNC_DEF VERDICT_REAL v_hex_shear_and_size( int num_nodes, VERDICT_REAL coordinates[][3] );

    //! Calculates hex distortion metric   
    /** {min(|J|)/actual volume}*parent volume, parent volume = 8 for hex.
       Reference --- SDRC/IDEAS Simulation: Finite Element Modeling--User's Guide */
    C_FUNC_DEF VERDICT_REAL v_hex_distortion( int num_nodes, VERDICT_REAL coordinates[][3] );

/* individual quality functions for tet elements */

    //! Sets average size (volume) of tet, needed for v_tet_relative_size(...)
    C_FUNC_DEF void v_set_tet_size( VERDICT_REAL size );

    //! Calculates tet aspect ratio metric.
    /** CR / (3.0 * IR)  where CR = circumsphere radius, IR = inscribed sphere radius.
       Reference ---  V. N. Parthasarathy et al, A comparison of tetrahedron 
       quality measures, Finite Elem. Anal. Des., Vol 15(1993), 255-261. */ 
    C_FUNC_DEF VERDICT_REAL v_tet_aspect_beta( int num_nodes, VERDICT_REAL coordinates[][3] ); 

    //! Calculates tet aspect ratio gamma metric.
    /**  Srms**3 / (8.479670*V) where Srms = sqrt(Sum(Si**2)/6), Si = edge length. 
       Reference ---  V. N. Parthasarathy et al, A comparison of tetrahedron 
       quality measures, Finite Elem. Anal. Des., Vol 15(1993), 255-261. */ 
    C_FUNC_DEF VERDICT_REAL v_tet_aspect_gamma( int num_nodes, VERDICT_REAL coordinates[][3] ); 

    //! Calculates tet volume.
    /** (1/6) * Jacobian at corner node.
       Reference ---  V. N. Parthasarathy et al, A comparison of tetrahedron 
       quality measures, Finite Elem. Anal. Des., Vol 15(1993), 255-261. */ 
    C_FUNC_DEF VERDICT_REAL v_tet_volume( int num_nodes, VERDICT_REAL coordinates[][3] ); 

    //! Calculates tet condition metric.
    /** Condition number of the Jacobian matrix at any corner. 
       Reference --- P. Knupp, Achieving Finite Element Mesh Quality via 
       Optimization of the Jacobian Matrix Norm and Associated Quantities,
       Intl. J. Numer. Meth. Engng. 2000, 48:1165-1185. */
    C_FUNC_DEF VERDICT_REAL v_tet_condition( int num_nodes, VERDICT_REAL coordinates[][3] ); 

    //! Calculates tet jacobian. 
    /** Minimum pointwise volume at any corner. 
       Reference --- P. Knupp, Achieving Finite Element Mesh Quality via 
       Optimization of the Jacobian Matrix Norm and Associated Quantities,
       Intl. J. Numer. Meth. Engng. 2000, 48:1165-1185. */
    C_FUNC_DEF VERDICT_REAL v_tet_jacobian( int num_nodes, VERDICT_REAL coordinates[][3] ); 

    //! Calculates tet scaled jacobian. 
    /** Minimum Jacobian divided by the lengths of 3 edge vectors 
       Reference --- P. Knupp, Achieving Finite Element Mesh Quality via 
       Optimization of the Jacobian Matrix Norm and Associated Quantities,
       Intl. J. Numer. Meth. Engng. 2000, 48:1165-1185. */
    C_FUNC_DEF VERDICT_REAL v_tet_scaled_jacobian( int num_nodes, VERDICT_REAL coordinates[][3] ); 

    //! Calculates tet shape metric.
    /** 3/Mean Ratio of weighted Jacobian matrix.
       Reference --- P. Knupp, Algebraic Mesh Quality Metrics for
       Unstructured Initial Meshes, submitted for publication. */ 
    C_FUNC_DEF VERDICT_REAL v_tet_shape( int num_nodes, VERDICT_REAL coordinates[][3] ); 

    //! Calculates tet relative size metric.
    /** Min( J, 1/J ), where J is determinant of weighted Jacobian matrix.
       Reference --- P. Knupp, Algebraic Mesh Quality Metrics for
       Unstructured Initial Meshes, submitted for publication. */ 
    C_FUNC_DEF VERDICT_REAL v_tet_relative_size_squared( int num_nodes, VERDICT_REAL coordinates[][3] ); 

    //! Calculates tet shape-size metric.
    /** Product of Shape and Relative Size. 
       Reference --- P. Knupp, Algebraic Mesh Quality Metrics for
       Unstructured Initial Meshes, submitted for publication. */ 
    C_FUNC_DEF VERDICT_REAL v_tet_shape_and_size( int num_nodes, VERDICT_REAL coordinates[][3] ); 

    //! Calculates tet distortion metric.
    /** {min(|J|)/actual volume}*parent volume, parent volume = 1/6 for tet.
       Reference --- SDRC/IDEAS Simulation: Finite Element Modeling--User's Guide */ 
    C_FUNC_DEF VERDICT_REAL v_tet_distortion( int num_nodes, VERDICT_REAL coordinates[][3] ); 
    
/* individual quality functions for pyramid elements */ 

    //! Calculates pyramid volume.
    C_FUNC_DEF VERDICT_REAL v_pyramid_volume( int num_nodes, VERDICT_REAL coordinates[][3] ); 


/* individual quality functions for wedge elements */

    //! Calculates wedge volume.
    C_FUNC_DEF VERDICT_REAL v_wedge_volume( int num_nodes, VERDICT_REAL coordinates[][3] ); 

   
/* individual quality functions for knife elements */

    //! Calculates knife volume.
    C_FUNC_DEF VERDICT_REAL v_knife_volume( int num_nodes, VERDICT_REAL coordinates[][3] ); 

    
/* individual quality functions for edge elements */

    //! Calculates edge length. 
    C_FUNC_DEF VERDICT_REAL v_edge_length( int num_nodes, VERDICT_REAL coordinates[][3] ); 

    
/* individual quality functions for quad elements */
    //! Sets average size (area) of quad, needed for v_quad_relative_size(...)
    C_FUNC_DEF void v_set_quad_size( VERDICT_REAL size );

    //! Calculates quad aspect ratio metric.
    /** Maximum edge length ratios at quad center.
       Reference --- J. Robinson, CRE Method of element testing and the 
       Jacobian shape parameters, Eng. Comput., Vol 4, 1987. */ 
    C_FUNC_DEF VERDICT_REAL v_quad_aspect( int num_nodes, VERDICT_REAL coordinates[][3] ); 

    //! Calculates quad skew metric.
    /** Maximum |cos A| where A is the angle between edges at quad center. 
       Reference --- J. Robinson, CRE Method of element testing and the 
       Jacobian shape parameters, Eng. Comput., Vol 4, 1987. */ 
    C_FUNC_DEF VERDICT_REAL v_quad_skew( int num_nodes, VERDICT_REAL coordinates[][3] ); 

    //! Calculates quad taper metric.
    /** Maximum ratio of lengths derived from opposite edges. 
       Reference --- J. Robinson, CRE Method of element testing and the 
       Jacobian shape parameters, Eng. Comput., Vol 4, 1987. */ 
    C_FUNC_DEF VERDICT_REAL v_quad_taper( int num_nodes, VERDICT_REAL coordinates[][3] ); 

    //! Calculates quad warpage metric.
    /** Cosine of Minimum Dihedral Angle formed by Planes Intersecting in Diagonals. 
       Reference --- J. Robinson, CRE Method of element testing and the 
       Jacobian shape parameters, Eng. Comput., Vol 4, 1987. */ 
    C_FUNC_DEF VERDICT_REAL v_quad_warpage( int num_nodes, VERDICT_REAL coordinates[][3] ); 

    //! Calculates quad area.
    /** Jacobian at quad center.
       Reference --- J. Robinson, CRE Method of element testing and the 
       Jacobian shape parameters, Eng. Comput., Vol 4, 1987. */ 
    C_FUNC_DEF VERDICT_REAL v_quad_area( int num_nodes, VERDICT_REAL coordinates[][3] ); 

    //! Calculates quad strech metric.
    /** Sqrt(2) * minimum edge length / maximum diagonal length.
       Reference --- FIMESH code. */
    C_FUNC_DEF VERDICT_REAL v_quad_stretch( int num_nodes, VERDICT_REAL coordinates[][3] ); 

    //! Calculates quad's smallest angle.
    /** Smallest included quad angle (degrees).
       Reference --- Unknown. */ 
    C_FUNC_DEF VERDICT_REAL v_quad_minimum_angle( int num_nodes, VERDICT_REAL coordinates[][3] ); 

    //! Calculates quad's largest angle.
    /** Largest included quad angle (degrees). 
       Reference --- Unknown. */ 
    C_FUNC_DEF VERDICT_REAL v_quad_maximum_angle( int num_nodes, VERDICT_REAL coordinates[][3] ); 

    //! Calculates quad oddy metric.
    C_FUNC_DEF VERDICT_REAL v_quad_oddy( int num_nodes, VERDICT_REAL coordinates[][3] ); 

    //! Calculates quad condition number metric.
    /** Maximum condition number of the Jacobian matrix at 4 corners.
       Reference --- P. Knupp, Achieving Finite Element Mesh Quality via 
       Optimization of the Jacobian Matrix Norm and Associated Quantities,
       Intl. J. Numer. Meth. Engng. 2000, 48:1165-1185. */
    C_FUNC_DEF VERDICT_REAL v_quad_condition( int num_nodes, VERDICT_REAL coordinates[][3] ); 

    //! Calculates quad jacobian.
    /** Minimum pointwise volume of local map at 4 corners & center of quad. 
       Reference --- P. Knupp, Achieving Finite Element Mesh Quality via 
       Optimization of the Jacobian Matrix Norm and Associated Quantities,
       Intl. J. Numer. Meth. Engng. 2000, 48:1165-1185. */
    C_FUNC_DEF VERDICT_REAL v_quad_jacobian( int num_nodes, VERDICT_REAL coordinates[][3] ); 

    //! Calculates quad scaled jacobian.
    /** Minimum Jacobian divided by the lengths of the 2 edge vectors. 
       Reference --- P. Knupp, Achieving Finite Element Mesh Quality via 
       Optimization of the Jacobian Matrix Norm and Associated Quantities,
       Intl. J. Numer. Meth. Engng. 2000, 48:1165-1185. */
    C_FUNC_DEF VERDICT_REAL v_quad_scaled_jacobian( int num_nodes, VERDICT_REAL coordinates[][3] ); 

    //! Calculates quad shear metric.
    /** 2/Condition number of Jacobian Skew matrix.
       Reference --- P. Knupp, Algebraic Mesh Quality Metrics for
       Unstructured Initial Meshes, submitted for publication. */
    C_FUNC_DEF VERDICT_REAL v_quad_shear( int num_nodes, VERDICT_REAL coordinates[][3] ); 

    //! Calculates quad shape metric.
    /** 2/Condition number of weighted Jacobian matrix. 
       Reference --- P. Knupp, Algebraic Mesh Quality Metrics for
       Unstructured Initial Meshes, submitted for publication. */
    C_FUNC_DEF VERDICT_REAL v_quad_shape( int num_nodes, VERDICT_REAL coordinates[][3] ); 

    //! Calculates quad relative size metric.
    /** Min( J, 1/J ), where J is determinant of weighted Jacobian matrix. 
       Reference --- P. Knupp, Algebraic Mesh Quality Metrics for
       Unstructured Initial Meshes, submitted for publication. */
    C_FUNC_DEF VERDICT_REAL v_quad_relative_size_squared( int num_nodes, VERDICT_REAL coordinates[][3] ); 

    //! Calculates quad shape-size metric.
    /** Product of Shape and Relative Size. 
       Reference --- P. Knupp, Algebraic Mesh Quality Metrics for
       Unstructured Initial Meshes, submitted for publication. */
    C_FUNC_DEF VERDICT_REAL v_quad_shape_and_size( int num_nodes, VERDICT_REAL coordinates[][3] ); 

    //! Calculates quad shear-size metric.
    /** Product of Shear and Relative Size. 
       Reference --- P. Knupp, Algebraic Mesh Quality Metrics for
       Unstructured Initial Meshes, submitted for publication. */
    C_FUNC_DEF VERDICT_REAL v_quad_shear_and_size( int num_nodes, VERDICT_REAL coordinates[][3] ); 

    //! Calculates quad distortion metric.
    /** {min(|J|)/actual area}*parent area, parent area = 4 for quad.
       Reference --- SDRC/IDEAS Simulation: Finite Element Modeling--User's Guide */
    C_FUNC_DEF VERDICT_REAL v_quad_distortion( int num_nodes, VERDICT_REAL coordinates[][3] ); 


/* individual quality functions for tri elements */

    //! Sets average size (area) of tri, needed for v_tri_relative_size(...) 
    C_FUNC_DEF void  v_set_tri_size( VERDICT_REAL size );

    //! Calculates tri metric.
    /**  */
    C_FUNC_DEF VERDICT_REAL v_tri_aspect( int num_nodes, VERDICT_REAL coordinates[][3] ); 

    //! Calculates tri metric.
    /** Maximum included angle in triangle */
    C_FUNC_DEF VERDICT_REAL v_tri_area( int num_nodes, VERDICT_REAL coordinates[][3] ); 

    //! Calculates tri metric.
    /** Minimum included angle in triangle */
    C_FUNC_DEF VERDICT_REAL v_tri_minimum_angle( int num_nodes, VERDICT_REAL coordinates[][3] ); 

    //! Calculates tri metric.
    /** Maximum included angle in triangle */
    C_FUNC_DEF VERDICT_REAL v_tri_maximum_angle( int num_nodes, VERDICT_REAL coordinates[][3] ); 

    //! Calculates tri metric.
    /** Condition number of the Jacobian matrix.
       Reference --- P. Knupp, Achieving Finite Element Mesh Quality via 
       Optimization of the Jacobian Matrix Norm and Associated Quantities,
       Intl. J. Numer. Meth. Engng. 2000, 48:1165-1185. */ 
    C_FUNC_DEF VERDICT_REAL v_tri_condition( int num_nodes, VERDICT_REAL coordinates[][3] ); 

    //! Calculates tri metric.
    /** Minimum Jacobian divided by the lengths of 2 edge vectors. 
       Reference --- P. Knupp, Achieving Finite Element Mesh Quality via 
       Optimization of the Jacobian Matrix Norm and Associated Quantities,
       Intl. J. Numer. Meth. Engng. 2000, 48:1165-1185. */ 
    C_FUNC_DEF VERDICT_REAL v_tri_scaled_jacobian( int num_nodes, VERDICT_REAL coordinates[][3] );

    //! Calculates tri metric.
    /**  */
    C_FUNC_DEF VERDICT_REAL v_tri_shear( int num_nodes, VERDICT_REAL coordinates[][3] ); 
    
    //! Calculates tri metric.
    /** Min( J, 1/J ), where J is determinant of weighted Jacobian matrix. 
       Reference ---  P. Knupp, Algebraic Mesh Quality Metrics for
       Unstructured Initial Meshes, submitted for publication. */
    C_FUNC_DEF VERDICT_REAL v_tri_relative_size_squared( int num_nodes, VERDICT_REAL coordinates[][3] ); 

    //! Calculates tri metric.
    /** 2/Condition number of weighted Jacobian matrix. 
       Reference ---  P. Knupp, Algebraic Mesh Quality Metrics for
       Unstructured Initial Meshes, submitted for publication. */
    C_FUNC_DEF VERDICT_REAL v_tri_shape( int num_nodes, VERDICT_REAL coordinates[][3] ); 

    //! Calculates tri metric.
    /**  Product of Shape and Relative Size. 
       Reference ---  P. Knupp, Algebraic Mesh Quality Metrics for
       Unstructured Initial Meshes, submitted for publication. */
    C_FUNC_DEF VERDICT_REAL v_tri_shape_and_size( int num_nodes, VERDICT_REAL coordinates[][3] );

    //! Calculates tri metric.
    /** {min(|J|)/actual area}*parent area, parent area = 1/2 for triangular element. 
       Reference --- SDRC/IDEAS Simulation: Finite Element Modeling--User's Guide */
    C_FUNC_DEF VERDICT_REAL v_tri_distortion( int num_nodes, VERDICT_REAL coordinates[][3] );



#endif  /* VERDICT_INC_LIB */



