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
//                             avtMinMaxQuery.h                              //
// ************************************************************************* //

#ifndef AVT_MINMAX_QUERY_H
#define AVT_MINMAX_QUERY_H
#include <query_exports.h>

#include <avtDatasetQuery.h>

#include <MinMaxInfo.h>

#include <string>

class avtMatrix;
class vtkDataArray;
class vtkDataSet;


// ****************************************************************************
//  Class: avtMinMaxQuery
//
//  Purpose:
//    A query that retrieves min and max information about a variable. 
//
//  Programmer: Kathleen Bonnell
//  Creation:   July 23, 2003
//
//  Modifications:
//    Kathleen Bonnell, Tue Feb  3 17:54:19 PST 2004
//    Renamed from avtPlotMinMaxQuery. Made into parent class. 
//
//    Kathleen Bonnell, Wed Mar 31 16:07:50 PST 2004 
//    Added args to constructor. 
//
//    Kathleen Bonnell, Tue Jul  6 16:59:26 PDT 2004 
//    Encapsulated elNum, vals, domain, coords in class MinMaxInfo.
//    Removed CreateMinMessage, CreateMaxMessage.  
//    Added InfoToString, CreateMessage, FindElement, FinalizeZoneCoord,
//    FinalizeNodeCoord.
//
//    Mark C. Miller, Tue Mar 27 08:39:55 PDT 2007
//    Added support for node origin
//
//    Hank Childs, Sun Dec 26 12:13:19 PST 2010
//    Add support for time parallelization.
//
//    Kathleen Biagas, Tue Jul 26 13:48:11 PDT 2011
//    Add GetDefaultInputParams.
//
// ****************************************************************************

class QUERY_API avtMinMaxQuery : virtual public avtDatasetQuery
{
  public:
                            avtMinMaxQuery(bool mn = true, bool mx = true);
    virtual                ~avtMinMaxQuery();


    virtual const char     *GetType(void)   
                                { return "avtMinMaxQuery"; };
    virtual const char     *GetDescription(void)
                                { return "Calculating min/max."; };

    static  void            GetDefaultInputParams(MapNode &);

  protected:
    virtual void            Execute(vtkDataSet *, const int);
    virtual void            PreExecute(void);
    virtual void            PostExecute(void);
            void            StandardPostExecute(void);
            void            TimeVaryingPostExecute(void);
    virtual void            VerifyInput(void);   
            void            Preparation(avtDataObject_p);   

  private:

    int                     dimension;
    int                     topoDim;
    int                     blockOrigin;
    int                     cellOrigin;
    int                     nodeOrigin;
    bool                    singleDomain;
    bool                    scalarCurve;
    bool                    nodeCentered;
    std::string             minMsg;
    std::string             maxMsg;
    std::string             elementName;

    MinMaxInfo              minInfo1;
    MinMaxInfo              minInfo2;
    MinMaxInfo              maxInfo1;
    MinMaxInfo              maxInfo2;

    std::string             nodeMsg1;
    std::string             nodeMsg2;
    std::string             zoneMsg1;
    std::string             zoneMsg2;

    bool                    doMin;
    bool                    doMax;

    avtQueryableSource     *src;

    const avtMatrix        *invTransform;

    void                    GetNodeCoord(vtkDataSet *ds, const int id, 
                                        double coord[3]);
    void                    GetCellCoord(vtkDataSet *ds, const int id, 
                                        double coord[3]);

    void                    CreateResultMessage(const int);

    std::string             InfoToString(const MinMaxInfo &);
    void                    CreateMessage(const int, const MinMaxInfo &, 
                                          const MinMaxInfo &, std::string &,
                                          doubleVector &);

    void                    FindElement(MinMaxInfo &);
    void                    FinalizeZoneCoord(vtkDataSet *, 
                                              vtkDataArray *, 
                                              MinMaxInfo &, bool);
    void                    FinalizeNodeCoord(vtkDataSet *, MinMaxInfo &);
};


#endif
