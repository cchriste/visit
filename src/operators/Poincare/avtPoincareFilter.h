/*****************************************************************************
*
* Copyright (c) 2000 - 2016, Lawrence Livermore National Security, LLC
* Produced at the Lawrence Livermore National Laboratory
* LLNL-CODE-442911
* All rights reserved.
*
* This file is  part of VisIt. For  details, see https://visit.llnl.gov/.  The
* full copyright notice is contained in the file COPYRIGHT located at the root
* of the VisIt distribution or at http://www.llnl.gov/visit/copyright.html.
*
* Redistribution  and  use  in  source  and  binary  forms,  with  or  without
* modification, are permitted provided that the following conditions are met:
*
*  - Redistributions of  source code must  retain the above  copyright notice,
*    this list of conditions and the disclaimer below.
*  - Redistributions in binary form must reproduce the above copyright notice,
*    this  list of  conditions  and  the  disclaimer (as noted below)  in  the
*    documentation and/or other materials provided with the distribution.
*  - Neither the name of  the LLNS/LLNL nor the names of  its contributors may
*    be used to endorse or promote products derived from this software without
*    specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT  HOLDERS AND CONTRIBUTORS "AS IS"
* AND ANY EXPRESS OR  IMPLIED WARRANTIES, INCLUDING,  BUT NOT  LIMITED TO, THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND  FITNESS FOR A PARTICULAR  PURPOSE
* ARE  DISCLAIMED. IN  NO EVENT  SHALL LAWRENCE  LIVERMORE NATIONAL  SECURITY,
* LLC, THE  U.S.  DEPARTMENT OF  ENERGY  OR  CONTRIBUTORS BE  LIABLE  FOR  ANY
* DIRECT,  INDIRECT,   INCIDENTAL,   SPECIAL,   EXEMPLARY,  OR   CONSEQUENTIAL
* DAMAGES (INCLUDING, BUT NOT  LIMITED TO, PROCUREMENT OF  SUBSTITUTE GOODS OR
* SERVICES; LOSS OF  USE, DATA, OR PROFITS; OR  BUSINESS INTERRUPTION) HOWEVER
* CAUSED  AND  ON  ANY  THEORY  OF  LIABILITY,  WHETHER  IN  CONTRACT,  STRICT
* LIABILITY, OR TORT  (INCLUDING NEGLIGENCE OR OTHERWISE)  ARISING IN ANY  WAY
* OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
* DAMAGE.
*
*****************************************************************************/

// ************************************************************************* //
//                          avtPoincareFilter.h                              //
// ************************************************************************* //

#ifndef AVT_POINCARE_FILTER_H
#define AVT_POINCARE_FILTER_H

/** header file for plugin development */
#include <avtPluginFilter.h>
/** header file for parallel integral curve system via the streamline filter */
#include <avtPICSFilter.h>

/** included attributes for Poincare */
#include <PoincareAttributes.h>

#include <vtkPlane.h>

#include <vector>

class avtPoincareIC;

// ****************************************************************************
//  Class: avtPoincareFilter
//
//  Purpose:
//      This operator is the implied operator associated with an Poincare plot.
//
//  Programmer: Dave Pugmire -- generated by xml2avt
//  Creation:   Tue Oct 7 09:02:52 PDT 2008
//
//  Modifications:
//
// ****************************************************************************

class avtPoincareFilter : public virtual avtPluginFilter,
                          public virtual avtPICSFilter
{
  public:
    // default constructor
    avtPoincareFilter();
    // virtual destructor
    virtual            ~avtPoincareFilter();
    // create filter
    static avtFilter   *Create();

    virtual const char       *GetType(void) { return "avtPoincareFilter"; };
    virtual const char       *GetDescription(void) {
      return "Performing Poincare"; };

    void SetPointSource(const double *p);
    void SetPointListSource(const std::vector<double> &ptList);
    void SetLineSource(const double *p0, const double *p1,
                       int den, bool rand, int seed, int numPts);

    void SetPuncturePlane( unsigned int val ) { puncturePlane = val; }
    void SetPuncturePlotType( unsigned int val ) { puncturePlotType = val; }
    void SetPuncturePeriod( double val ) { puncturePeriod = val; }
    void SetPuncturePeriodTolerance( double val ) { puncturePeriodTolerance = val; }
    void SetAnalysis( unsigned int val ) { analysis = val; }

    void SetMaxPunctures( double punctures ) { maxPunctures = punctures; }

    void SetMaximumToroidalWinding( unsigned int value ) {
      maximumToroidalWinding = value; };
    void SetOverrideToroidalWinding( unsigned int value) {
      overrideToroidalWinding = value; }
    void SetOverridePoloidalWinding( unsigned int value) {
      overridePoloidalWinding = value; }

    void SetWindingPairConfidence( double val ) { windingPairConfidence = val; }
    void SetRationalSurfaceFactor( double val ) { rationalSurfaceFactor = val; }

    void SetOverlaps( unsigned int val ) { overlaps = val; }


    void SetShowCurves( unsigned int val ) { is_curvemesh = val; }

    void SetClipPlanes( std::vector< double > planeAngles ) {
      planes = planeAngles; }

    void SetDataValue( unsigned int val ) { dataValue = val; }

    void SetShowRationalSurfaces( bool val ) { showRationalSurfaces = val; }
    void SetRationalSurfaceMaxIterations( int val ) { rationalSurfaceMaxIterations = val; }

    void SetShowOPoints( bool val ) { showOPoints = val; }
    void SetOPointMaxIterations( int val ) { OPointMaxIterations = val; }

    void SetShowXPoints( bool val ) { showXPoints = val; }
    void SetXPointMaxIterations( int val ) { XPointMaxIterations = val; }

    void SetPerformOLineAnalysis( bool val ) { performOLineAnalysis = val; }
    void SetOLineToroidalWinding( int val ) { OLineToroidalWinding = val; }
    void SetOLineAxisFileName( std::string val ) { OLineAxisFileName = val; }

    void SetShowChaotic( bool val ) { showChaotic = val; }
    void SetShowIslands( bool val ) { showIslands = val; }
    void SetShowLines( bool val )   { showLines = val; }
    void SetShowPoints( bool val )  { showPoints = val; }
    void SetPointScale(int scale)   { pointScale = scale; }
    void SetSummaryFlag( bool val ) { summaryFlag = val; }
    void SetVerboseFlag( bool val ) { verboseFlag = val; }
    void SetShow1DPlots( bool val )   { show1DPlots = val; }

    // Methods to set the filter's attributes.
    void                 SetIntersectionCriteria();

    virtual void         SetAtts(const AttributeGroup*);
    virtual bool         Equivalent(const AttributeGroup*);

  protected:
    // Fieldline overides.
    virtual void ExamineContract(avtContract_p in_contract);
    virtual avtContract_p ModifyContract(avtContract_p in_spec);

    virtual void Execute(void);
    virtual bool ContinueExecute();
    virtual void PreExecute(void);
    virtual void PostExecute(void);

    virtual void UpdateDataObjectInfo(void);
    virtual void ReportWarnings(std::vector<avtIntegralCurve *> &ics);

    virtual void GetIntegralCurvePoints(std::vector<avtIntegralCurve *> &ic);
    /** virtual functions for avtPICSFilter */
    /** construct default integral curve */
    virtual avtIntegralCurve *CreateIntegralCurve();
    /** construct integral curve given parameters */
    virtual avtIntegralCurve *CreateIntegralCurve( const avtIVPSolver* model,
                                                   avtIntegralCurve::Direction,
                                                   const double& t_start,
                                                   const avtVector &p_start,
                                                   const avtVector &v_start,
                                                   long ID );

    virtual bool             GetAllSeedsSentToAllProcs() { return true; };

    /** Construct the initial locations to emanate integral curves */
    virtual std::vector<avtVector> GetInitialLocations();
    virtual std::vector<avtVector> GetInitialVelocities();

    void GenerateSeedPointsFromPoint(std::vector<avtVector> &pts);
    void GenerateSeedPointsFromPointList(std::vector<avtVector> &pts);
    void GenerateSeedPointsFromLine(std::vector<avtVector> &pts);

    void CreateIntegralCurveOutput(std::vector<avtIntegralCurve*,
                                   std::allocator<avtIntegralCurve*> >&) {};

    /** type of communication between processors
      options include: RestoreSequence,LeaveOnCurrentProcessor, and 
             ReturnToOriginatingProcessor(NOT IMPLEMENTED) */
//    virtual CommunicationPattern GetCommunicationPattern()
//                { return LeaveOnCurrentProcessor; }
    virtual CommunicationPattern GetCommunicationPattern()
                { return ReturnToOriginatingProcessor; }

    virtual void drawPoints( avtDataTree *dt,
                             std::vector < avtVector > &nodes );

    virtual void drawRationalCurve( avtDataTree *dt,
                                    std::vector< std::vector < std::vector < avtVector > > > &nodes,
                                    unsigned int islands,
                                    unsigned int skip,
                                    unsigned int color,
                                    double color_value );
  
    virtual void drawIrrationalCurve( avtDataTree *dt,
                                      std::vector< std::vector < std::vector < avtVector > > > &nodes,
                                      unsigned int nnodes,
                                      unsigned int islands,
                                      unsigned int skip,
                                      unsigned int color,
                                      double color_value,
                                      bool connect = false,
                                      bool modulo = false);
  
    virtual void drawSurface( avtDataTree *dt,
                              std::vector< std::vector < std::vector < avtVector > > > &nodes,
                              unsigned int islands,
                              unsigned int skip,
                              unsigned int color,
                              double color_value);
  
    virtual void drawPeriodicity( avtDataTree *dt,
                                  std::vector < avtVector > &nodes,
                                  unsigned int period,
                                  unsigned int nnodes,
                                  unsigned int islands,
                                  unsigned int poloidalWindings,
                                  unsigned int color,
                                  double color_value,
                                  bool ptFlag );

    // Poincare filter methods.
    bool ClassifyFieldlines(std::vector<avtIntegralCurve *> &ic);
    void SetupPlaneOrdering( avtPoincareIC *poincare_ic );
    void CreatePoincareOutput(avtDataTree *dt,
                              std::vector<avtIntegralCurve *> &ic);

    bool ClassifyRationals(std::vector<avtIntegralCurve *> &ic);
    void CreateRationalOutput(avtDataTree *dt,
                              std::vector<avtIntegralCurve *> &ic); /////////RATAIONAL

    void IssueWarningForMaxStepsTermination(bool v) 
                 { issueWarningForMaxStepsTermination = v; };
    void IssueWarningForStepsize(bool v) 
                 { issueWarningForStepsize = v; };
    void IssueWarningForStiffness(bool v) 
                 { issueWarningForStiffness = v; };
    void IssueWarningForCriticalPoints(bool v, double speed) 
                 { issueWarningForCriticalPoints = v;
                   criticalPointThreshold = speed; };

    std::string outVarName;

    PoincareAttributes atts;
    bool               needsRecalculation;

    int    sourceType;   

    int      maxSteps;
    bool     doTime;
    double   maxTime;

    // Various starting locations for streamlines.
    std::vector< avtVector > points;
    int       numSamplePoints;
    int       sampleDensity[3];
    bool      randomSamples;
    int       randomSeed;


    unsigned int puncturePlane;
    unsigned int puncturePlotType;
    double puncturePeriod;
    double puncturePeriodTolerance;
    unsigned int analysis;

    double maxPunctures;

    unsigned int maximumToroidalWinding;
    unsigned int overrideToroidalWinding;
    unsigned int overridePoloidalWinding;

    double windingPairConfidence;
    double rationalSurfaceFactor;

    std::vector< double > planes;

    unsigned int overlaps;
    bool is_curvemesh;

    unsigned int dataValue;
  
    vtkPlane *intPlane; 
    int maxIntersections;

    bool showRationalSurfaces, showOPoints, showXPoints, showIslands,
      showChaotic;
    bool showLines, showPoints, summaryFlag, verboseFlag, show1DPlots;
    int  pointScale;

    unsigned int rationalSurfaceMaxIterations;
    unsigned int OPointMaxIterations;
    unsigned int XPointMaxIterations;

    bool performOLineAnalysis;
    int  OLineToroidalWinding;
    std::string OLineAxisFileName;

    avtVector seedVelocity;

    //input seed points..
    std::vector<avtVector> seedPoints;

    bool      issueWarningForMaxStepsTermination;
    bool      issueWarningForStepsize;
    bool      issueWarningForStiffness;
    bool      issueWarningForCriticalPoints;
    double    criticalPointThreshold;
};
#endif
