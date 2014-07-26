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
//                       avtConnComponentsSummaryQuery.h                     //
// ************************************************************************* //

#ifndef AVT_CONN_COMPONENTS_SUMMARY_QUERY_H
#define AVT_CONN_COMPONENTS_SUMMARY_QUERY_H
#include <query_exports.h>

#include <avtDatasetQuery.h>
#include <avtConnComponentsQuery.h>

#include <string>
#include <vector>

class avtEdgeLength;
class avtRevolvedVolume;
class avtVMetricArea;
class avtVMetricVolume;

class vtkDataSet;


// ****************************************************************************
//  Class: avtConnComponentsSummaryQuery
//
//  Purpose:
//      Finds the connected components of a mesh and creates a point dataset 
//      with per component quantities. This dataset is saved in as an Xmdv
//      file.  
// 
//      Component Fields:
//       x,y,z              Component centroid location
//       comp_length        Component length (if 1D topo data)
//       comp_area          Component area   (if 2D topo data)
//       comp_volume        Component volume (if 3D topo or revolved volume data)
//       comp_sum           Component variable sum 
//       comp_weighted_sum  Component weighted variable sum
//       comp_bb_{x,y,z}_{min,max} Component bounding box values.
//
//
//  Programmer: Cyrus Harrison
//  Creation:   March 1, 2007 
//
//  Modifications:
//    Cyrus Harrison, Mon Mar 30 11:34:27 PDT 2009
//    Added support for per component bounding boxes.
//
//    Cyrus Harrison, Wed Dec 22 14:53:33 PST 2010
//    Added output for the number of processors each component spans.
//
//    Cyrus Harrison, Mon Jun  6 17:02:12 PDT 2011
//    Added lengthFilter.
//
//    Kathleen Biagas, Tue Jun 21 10:42:07 PDT 2011
//    Added SetInputParams.
//
//    Kathleen Biagas, Fri Jul 15 16:05:35 PDT 2011
//    Added GetDefaultInputParams.
//
//    Cyrus Harrison, Wed Jul 16 13:54:00 PDT 2014
//    Added PrepareMapNodeResult.
//
// ****************************************************************************

class QUERY_API avtConnComponentsSummaryQuery 
: public avtConnComponentsQuery
{
  public:
                               avtConnComponentsSummaryQuery();
    virtual                   ~avtConnComponentsSummaryQuery();

    virtual const char        *GetType(void)
                         { return "avtConnComponentsSummaryQuery"; };
    virtual const char        *GetDescription(void)
                         { return "Generating connected components summary.";};

    virtual void               SetInputParams(const MapNode &);
    static  void               GetDefaultInputParams(MapNode &);

    void                       SetOutputFileName(const std::string &fname)
                                 {outputFileName = fname; }
  protected:

    std::string                outputFileName;

    avtEdgeLength             *lengthFilter;
    avtRevolvedVolume         *revolvedVolumeFilter;
    avtVMetricArea            *areaFilter;
    avtVMetricVolume          *volumeFilter;

    std::string                variableName;

    std::vector<int>           nCellsPerComp;
    std::vector<int>           nProcsPerComp;

    std::vector<double>        xCentroidPerComp;
    std::vector<double>        yCentroidPerComp;
    std::vector<double>        zCentroidPerComp;

    // for computing per component bounds
    std::vector<double>        xMinPerComp;
    std::vector<double>        xMaxPerComp;
    std::vector<double>        yMinPerComp;
    std::vector<double>        yMaxPerComp;
    std::vector<double>        zMinPerComp;
    std::vector<double>        zMaxPerComp;

    std::vector<double>        lengthPerComp;
    std::vector<double>        areaPerComp;
    std::vector<double>        volPerComp;

    std::vector<double>        sumPerComp;
    std::vector<double>        wsumPerComp;

    bool                       findLength;
    bool                       findArea;
    bool                       findVolume;

    virtual void               Execute(vtkDataSet *, const int);
    virtual void               PreExecute(void);
    virtual void               PostExecute(void);
    virtual avtDataObject_p    ApplyFilters(avtDataObject_p);
    virtual void               VerifyInput(void);

    void                       SaveComponentResults(std::string fname);
    void                       PrepareComponentResults(std::vector<double> &);
    void                       PrepareMapNodeResult(MapNode &result_node);
};


#endif



