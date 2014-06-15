/*****************************************************************************
*
* Copyright (c) 2000 - 2014, Lawrence Livermore National Security, LLC
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
//                       avtConnComponentsExpression.h                       //
// ************************************************************************* //

#ifndef AVT_CONNECTED_COMPONENTS_EXPRESSION_H
#define AVT_CONNECTED_COMPONENTS_EXPRESSION_H

#include <avtExpressionFilter.h>

#include <vector>

class     vtkDataArray;
class     vtkIntArray;
class     vtkDataSet;
class     vtkCell;
class     vtkUnstructuredGrid;
class     avtIntervalTree;

// ****************************************************************************
//  Class: avtConnComponentsExpression
//
//  Purpose:
//      An expression that finds and labels the connected components
//      of an unstructured mesh.
//
//  Programmer: Cyrus Harrison
//  Creation:   January 22, 2007
//
//  Modifications:
//    Cyrus Harrison, Fri Mar 16 15:52:47 PDT 2007
//    Added variables to track progress. 
//
//    Cyrus Harrison, Sat Aug 11 14:41:01 PDT 2007
//    Added LabelGhostNeighbors and ModifyContract
//
//    Cyrus Harrison, Thu Aug 23 08:35:12 PDT 2007
//    Added ProcessArguments and option to explicitly disable ghost neighbors
//
//    Cyrus Harrison, Fri Oct 17 11:35:48 PDT 2008
//    Added CheckForProperGhostZones, a new sanity check for proper ghost info
//    before using the ghost zone neighbors enhancement in the parallel case.
//
//    Cyrus Harrison, Mon Mar 30 12:06:50 PDT 2009
//    Added variable name argument to BoundarySet::RelocateUsingPartition.
//    This allows us to ony send the label variable and reduce communication.
//
// ****************************************************************************

class EXPRESSION_API avtConnComponentsExpression : public avtExpressionFilter
{
  public:
                              avtConnComponentsExpression();
    virtual                  ~avtConnComponentsExpression();

    virtual const char       *GetType(void)
                                     { return "avtConnComponentsExpression"; };
    virtual const char       *GetDescription(void)
                                      {return "Finding connected components";};
    virtual int               NumVariableArguments() { return 1; };
    virtual int               GetVariableDimension() { return 1; };
    virtual bool              IsPointVariable(void) { return false; }
    virtual int               GetNumberOfComponents();
    virtual void              ProcessArguments(ArgsExpr*, ExprPipelineState *);

    // ************************************************************************
    //  Class: avtConnComponentsExpression::UnionFind
    //
    //  Purpose:
    //      Helps efficiently solve the disjoint set union problem.
    //
    //  Programmer: Cyrus Harrison
    //  Creation:   January 23, 2007
    //
    //  Modifications:
    //   Use pointers for faster access.
    //
    // ************************************************************************
    class UnionFind
    {
      public:
                      UnionFind(int,bool);
        virtual      ~UnionFind();

        int           Find(int);
        void          Union(int, int);

        bool          IsValid(int);
        void          SetValid(int,bool);

        int           FinalizeLabels();
        int           GetFinalLabel(int);
        int           GetNumberOfFinalLabels();


      protected:
        std::vector<int>      ranks;
        std::vector<int>      parents;
        std::vector<int>      valid;
        std::vector<int>      finalLabels;
        int                   nFinalLabels;

        int  *ranksPtr;
        int  *parentsPtr;
        int  *validPtr;
        int  *finalLabelsPtr;

    };

    class SpatialPartition;

    // ************************************************************************
    //  Class: avtConnComponentsExpression::BoundarySet
    //
    //  Purpose:
    //      Provides efficient access to cell intersections across sets.
    //
    //  Programmer: Cyrus Harrison
    //  Creation:   January 25, 2007
    //
    // ************************************************************************
    class BoundarySet
    {
      public:
                              BoundarySet();
        virtual              ~BoundarySet();

        void                  AddMesh(vtkDataSet *);
        void                  Finalize();
        void                  Clear();

        std::vector<vtkDataSet *>  GetMeshes() const;

        void                  GetBounds(double *) const;


        void                  GetIntersectionSet(int,
                                                 int,
                                                 std::vector<int> &,
                                                 std::vector<int> &) const;

        bool                  GetBoundsIntersection(double *,
                                                    double *,
                                                    double *) const;

        void                  RelocateUsingPartition(const SpatialPartition &,
                                                     const char *label_var_name); 


      protected:
        std::vector<vtkDataSet *>       meshes;
        std::vector<avtIntervalTree *>  itrees;

        bool                       empty;
        double                     bounds[6];

    };

    // ************************************************************************
    //  Class: SpatialPartition
    //
    //  Purpose:
    //      Finds an efficient way to partition meshes across available 
    //      processors.
    //
    //  Notes: Adapted from Hank Child's avtPosCMFEAlgorithm::SpatialPartition
    //
    //
    //  Programmer: Cyrus Harrison
    //  Creation:   February 2, 2007
    //
    // ************************************************************************
    class SpatialPartition
    {
      public:
                              SpatialPartition();
        virtual              ~SpatialPartition();

        void                  CreatePartition(const BoundarySet &, 
                                              double *);

        void                  GetProcessorList(vtkCell *,
                                               std::vector<int> &) const;

      protected:
        avtIntervalTree      *itree;
    };

    // ************************************************************************
    //  Class: avtConnComponentsExpression::RcbOcTree
    //
    //  Purpose:
    //      Creates an efficient data structure used when creating a spatial
    //      partition.
    //
    //  Programmer: Ryan Bleile
    //  Creation:   October 23, 2013
    //
    // ************************************************************************
    class RcbOcTree
    {
      public:
                      RcbOcTree();
                      RcbOcTree(double *min, double *max);
                      ~RcbOcTree();

        void          AddPoint( double *point , bool checkforduplicates = false );
        void          SetBB( double* min, double *max );
        void          SplitNode( double *point );

        int           GetNumPoints();
        int           NumPointsInRange( double* min, double* max );
        int           WhichChild( double* point );

      private:
        int           numPoints;
        double        BBMax[3];
        double        BBMin[3];
        double        min_extent[3];
        double        max_extent[3];

        double        *X;
        double        *Y;
        double        *Z;
        RcbOcTree     *child;
    };


  protected:
    int                       nFinalComps;
    int                       currentProgress;
    int                       totalSteps;

    int                       enableGhostNeighbors;

    virtual void              Execute(void);

    virtual avtContract_p
                              ModifyContract(avtContract_p);
    
    virtual bool              CheckForProperGhostZones(vtkDataSet **sets,int nsets);
    virtual void              LabelGhostNeighbors(vtkDataSet *);
    virtual void              LabelBoundaryNeighbors(vtkDataSet *);
    
    virtual vtkIntArray      *SingleSetLabel(vtkDataSet *, int &);

    virtual int               MultiSetResolve(int,
                                              const BoundarySet &,
                                              const std::vector<vtkDataSet*> &,
                                              const std::vector<vtkIntArray*> &);

    virtual void              MultiSetList(int,
                                           const BoundarySet &,
                                           const std::vector<vtkDataSet*> &,
                                           const std::vector<vtkIntArray*> &,
                                           std::vector<int> &,
                                           std::vector<int> &);

    virtual int               GlobalLabelShift(int,
                                               const std::vector<vtkIntArray*> &);

    virtual int               GlobalResolve(int,
                                            BoundarySet &,
                                            const std::vector<vtkDataSet*> &,
                                            const std::vector<vtkIntArray*> &);

    virtual int               GlobalUnion(int,
                                          const std::vector<int>&,
                                          const std::vector<int>&,
                                          const std::vector<vtkIntArray*> &);

    virtual void              ShiftLabels(vtkIntArray *, int);
};


#endif

