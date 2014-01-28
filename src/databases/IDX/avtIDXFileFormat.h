/*****************************************************************************
*
* Copyright (c) 2000 - 2013, Lawrence Livermore National Security, LLC
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
//                            avtIDXFileFormat.h                           //
// ************************************************************************* //

#ifndef AVT_IDX_FILE_FORMAT_H
#define AVT_IDX_FILE_FORMAT_H

#include <avtMTMDFileFormat.h>

#include <vector>

#include <visuscpp/appkit/visus_appkit.h>

///////////////////////////////////////////////////////////////////////////////
class avtIDXFileFormat;
class avtIDXQueryNode : public Visus::QueryNode
{
public:
    
    //constructor
    avtIDXQueryNode(avtIDXFileFormat *node_):node(node_) {VisusAssert(node);}
    
    //destructor
    ~avtIDXQueryNode() {}

    //publish
    virtual bool publish(const Visus::DictObject& evt);

private:

    avtIDXFileFormat *node;

};

// ****************************************************************************
//  Class: avtIDXFileFormat
//
//  Purpose:
//      Reads in IDX files as a plugin to VisIt.
//
//  Programmer: camc -- generated by xml2avt
//  Creation:   Tue Jan 7 12:20:07 MST 2014
//
// ****************************************************************************

class DummyNode;
class avtView3D;
class avtIDXFileFormat : public avtMTMDFileFormat, public Visus::Object
{
    friend class avtIDXQueryNode;

  public:
                       avtIDXFileFormat(const char *);
    virtual           ~avtIDXFileFormat();

    //
    // This is used to return unconvention data -- ranging from material
    // information to information about block connectivity.
    //
    // virtual void      *GetAuxiliaryData(const char *var, int timestep, 
    //                                     int domain, const char *type, void *args, 
    //                                     DestructorFunction &);
    //

    //
    // If you know the times and cycle numbers, overload this function.
    // Otherwise, VisIt will make up some reasonable ones for you.
    //
    // virtual void        GetCycles(std::vector<int> &);
    // virtual void        GetTimes(std::vector<double> &);
    //

    virtual int            GetNTimesteps(void);

    virtual const char    *GetType(void)   { return "IDX"; };

    virtual bool           CanCacheVariable(const char *);

    //virtual bool          HasInvariantMetaData(void) const { return false; }

    virtual void           RegisterDataSelections(
                             const std::vector<avtDataSelection_p>&,
                             std::vector<bool>* applied);

    virtual void           PopulateDatabaseMetaData(avtDatabaseMetaData *, int);

    virtual vtkDataSet    *GetMesh(int, int, const char *);
    virtual vtkDataArray  *GetVar(int, int, const char *);
    virtual vtkDataArray  *GetVectorVar(int, int, const char *);

    virtual void           FreeUpResources(void); 

  protected:

    std::string            filename;

    // int                    meshNx, meshNy;
    // double                 meshXmin, meshXmax, meshYmin, meshYmax;
    // int                    coarseNx, coarseNy;

    std::vector<avtDataSelection_p>  selectionsList;
    std::vector<bool>               *selectionsApplied;
    
    void                   CalculateMesh(/*double &, double &,
                                           double &, double &, */int timestate);

    //
    //<ctc> old ViSUS plugin stuff, may not all be necessary
    //

    double                 frustum2d[6];
    double                 frustum3d[6];

    int                    resolution;

    int                    nblocks;
    int                    rank;

    bool                    haveData;
    int                     dim;         //2d or 3d
    int                     bounds[3]; 
    double                  extents[6];
    double                  fullextents[6];

    static int              num_instances;

    UniquePtr<Visus::Application> app;
    UniquePtr<Visus::Dataflow>    dataflow;
    void                          onDataflowInput(Visus::DataflowNode* dnode);
    // bool connect(Visus::DataflowPortPtr oport,Visus::DataflowPortPtr iport);
    // bool connect(Visus::DataflowNodePtr from,Visus::String oport,Visus::String iport,Visus::DataflowNodePtr to);

    SharedPtr<Visus::Dataset>       dataset;
    SharedPtr<Visus::DatasetNode>   dataflow_dataset;

    avtIDXQueryNode         *query; //<ctc> memory?

    SharedPtr<Visus::Array>  data;

    // Visus::KDataflowNodePtr  K                      (Visus::ObjectPtr value,Visus::String oport="value") {return dataflow->K(value,oport);}

    VISUS_DECLARE_BINDABLE(avtIDXFileFormat);
};


#endif
