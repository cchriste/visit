/*
   For more information, please see: http://software.sci.utah.edu

   The MIT License

   Copyright (c) 2009 Scientific Computing and Imaging Institute,
   University of Utah.

   License for the specific language governing rights and limitations under
   Permission is hereby granted, free of charge, to any person obtaining a
   copy of this software and associated documentation files (the "Software"),
   to deal in the Software without restriction, including without limitation
   the rights to use, copy, modify, merge, publish, distribute, sublicense,
   and/or sell copies of the Software, and to permit persons to whom the
   Software is furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included
   in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
   OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
   THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
   DEALINGS IN THE SOFTWARE.
*/


/*
 *  StreamlineAnalyzer.h:
 *
 *  Written by:
 *   Allen Sanderson
 *   SCI Institute
 *   University of Utah
 *   September 2006
 *
 */

#if !defined(StreamlineAnalyzerLib_h)
#define StreamlineAnalyzerLib_h

//#define STRAIGHTLINE_SKELETON 1

//#include <Core/Geometry/Point.h>
//#include <Core/Geometry/Vector.h>

#include <avtVector.h>
#include <DebugStream.h>

typedef avtVector Point ;
typedef avtVector Vector;

#include <vector>

using namespace std;


class FieldlineInfo {

public:

enum FieldlineType { UNKNOWN  = 0,

                     PERIODIC = 10,
                     RATIONAL = 11,
                     O_POINT  = 12,
                     X_POINT  = 13,
                     
                     QUASI_PERIODIC = 20,
                     FLUX_SURFACE   = 21,
                     ISLAND_CHAIN   = 22,
                     
                     CHAOTIC = 30 };
  
enum analysisStatus { IN_PROCESS = 0x0001,
                      COMPLETED  = 0x0002,
                      TERMINATED = 0x0004,

                      DELETE     = 0x0008,

                      ADD          = 0x00F0,
                      ADD_O_POINTS = 0x0010,
                      ADD_X_POINTS = 0x0020 };

public:

  FieldlineType type;

  int analysisStatus;

  unsigned int toroidalWinding;
  unsigned int poloidalWinding;

  unsigned int toroidalPeriod;
  unsigned int poloidalPeriod;

  unsigned int windingGroupOffset;
  unsigned int islands;

  float nnodes;

  float confidence;
  float ridgelineVariance;


  unsigned int nPuncturesNeeded;

  std::vector< Point > OPoints;
  bool seedOPoints;

  std::vector< int > parentIds;
  std::vector< int > childIds;
};

namespace FusionPSE {

//using namespace SCIRun;

class FieldlineLib
{
public:

  Point interpert( Point lastPt, Point currPt, double t );

  int ccw( Vector v0, Vector v1 );

  int intersect( Point l0_p0, Point l0_p1,
                 Point l1_p0, Point l1_p1 );

  void convexHull( vector< pair< Point, unsigned int > > &hullPts,
                   unsigned int &m,
                   unsigned int toroidalWinding,
                   int dir );

  bool hullCheck( vector< Point > &points,
                  int &direction);

  unsigned int GCD( unsigned int a, unsigned int b );

  Point circle(Point &pt1, Point &pt2, Point &pt3);

  bool IsPerpendicular(Point &pt1, Point &pt2, Point &pt3);

  Point CalcCircle(Point &pt1, Point &pt2, Point &pt3);

  bool
  IntersectCheck( vector< Point >& points, unsigned int nbins );

  unsigned int Blankinship( unsigned int toroidalWinding,
                            unsigned int poloidalWinding,
                            unsigned int offset = 1 );

  void poloidalWindingCheck( vector< unsigned int > &poloidalWindingset,
                             vector< pair < pair<unsigned int, unsigned int >,
                                            double > > &windingSetList );

  unsigned int
  periodicityStats( vector< Point >& points,
                    vector< pair< unsigned int, double > >& stats );


  unsigned int
  periodicityChecks( vector< Point >& points,
                     vector< pair< unsigned int, double > >& stats,
                     double &consistency,
                     bool useBest );

  double
  calculateSumOfSquares( vector< Point >& poloidalWinding_points,
                         unsigned int poloidalWinding,
                         int checkType );

  bool
  rationalCheck( vector< Point >& points,
                 unsigned int toroidalWinding,
                 unsigned int &nnodes,
                 float delta=0.01);

  bool
  islandChecks( vector< Point >& points,
                unsigned int toroidalWinding,
                unsigned int &islands,
                unsigned int &nnodes,
                bool &complete );

  FieldlineInfo 
  fieldlineProperties( vector< Point > &ptList,
                       unsigned int override,
                       unsigned int maxToroidalWinding,
                       double windingPairConfidence,
                       double periodicityConsistency,
                       bool findIslandCenters );

  unsigned int
  islandProperties( vector< Point > &points,
                    Vector &baseCentroid,
                    unsigned int &startIndex,
                    unsigned int &middleIndex,
                    unsigned int &stopIndex,
                    unsigned int &nodes );

  unsigned int
  surfaceOverlapCheck( vector< vector< Point > > &bins,
                    unsigned int toroidalWinding,
                    unsigned int skip,
                    unsigned int &nnodes );

  unsigned int
  surfaceGroupCheck( vector< vector< Point > > &bins,
                     unsigned int i,
                     unsigned int j,
                     unsigned int nnodes );

  unsigned int
  removeOverlap( vector< vector < Point > > &bins,
                 unsigned int &nnodes,
                 unsigned int toroidalWinding,
                 unsigned int poloidalWinding,
                 unsigned int skip,
                 unsigned int island );

  unsigned int
  smoothCurve( vector< vector < Point > > &bins,
               unsigned int &nnodes,
               unsigned int toroidalWinding,
               unsigned int poloidalWinding,
               unsigned int skip,
               unsigned int island );

  unsigned int
  mergeOverlap( vector< vector < Point > > &bins,
                unsigned int &nnodes,
                unsigned int toroidalWinding,
                unsigned int poloidalWinding,
                unsigned int skip,
                unsigned int island );

  void
  findIslandCenter( vector< Point > &points,
                    unsigned int nnodes,
                    unsigned int toroidalWinding,
                    unsigned int poloidalWinding,
                    vector< Point > &centers );

  bool verboseFlag;
};

} // End namespace FusionPSE

#endif // StreamlineAnalyzerLib_h
