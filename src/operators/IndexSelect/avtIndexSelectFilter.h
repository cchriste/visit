/*****************************************************************************
*
* Copyright (c) 2000 - 2008, Lawrence Livermore National Security, LLC
* Produced at the Lawrence Livermore National Laboratory
* LLNL-CODE-400142
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
//  File: avtIndexSelectFilter.h
// ************************************************************************* //

#ifndef AVT_IndexSelect_FILTER_H
#define AVT_IndexSelect_FILTER_H


#include <avtPluginStreamer.h>
#include <IndexSelectAttributes.h>


class vtkDataSet;
class vtkMaskPoints;
class vtkVisItExtractGrid;
class vtkVisItExtractRectilinearGrid;


// ****************************************************************************
//  Class: avtIndexSelectFilter
//
//  Purpose:
//      A plugin operator for IndexSelect.
//
//  Programmer: childs -- generated by xml2info
//  Creation:   Wed Jun 5 09:09:10 PDT 2002
//
//  Modifications:
//
//    Hank Childs, Sat Jun 29 16:22:48 PDT 2002
//    Added support for groups.
//
//    Mark C. Miller, Tue Sep 28 19:57:42 PDT 2004
//    Added data selection id
//
//    Kathleen Bonnell, Tue May 10 11:19:24 PDT 2005 
//    Use VisIt versions of vtkExtractGrid and vtkExtractRectilinearGrid, 
//    they have been modified to correctly handle cell data when VOI is
//    along max boundary. 
//
//    Kathleen Bonnell, Thu Aug  4 15:47:59 PDT 2005 
//    Added UpdateDataObjectInfo.
//
//    Kathleen Bonnell,  Mon Jan 30 15:10:26 PST 2006 
//    Add vtkMaskPoints for a points filter. 
//
//    Jeremy Meredith, Wed Jan 17 11:41:51 EST 2007
//    Added support for transformed rectilinear grids.
//
//    Kathleen Bonnell, Thu Jun 21 16:31:59 PDT 2007 
//    Added amrLevel, amrMesh, int* arg to PrepareFilters.
//
// ****************************************************************************

class avtIndexSelectFilter : public avtPluginStreamer
{
  public:
                         avtIndexSelectFilter();
    virtual             ~avtIndexSelectFilter();

    static avtFilter    *Create();

    virtual const char  *GetType(void)  { return "avtIndexSelectFilter"; };
    virtual const char  *GetDescription(void)
                             { return "Index selecting"; };
    virtual void         ReleaseData(void);

    virtual void         SetAtts(const AttributeGroup*);
    virtual bool         Equivalent(const AttributeGroup*);

  protected:
    IndexSelectAttributes       atts;
    bool                        haveIssuedWarning;
    bool                        successfullyExecuted;
    int                         selID;
    bool                        groupCategory;
    bool                        amrMesh;
    int                         amrLevel;

    vtkVisItExtractGrid                  *curvilinearFilter;
    vtkVisItExtractRectilinearGrid       *rectilinearFilter;
    vtkMaskPoints                        *pointsFilter;

    void                        PrepareFilters(int [3], int *);

    virtual vtkDataSet         *ExecuteData(vtkDataSet *, int, std::string);
    virtual void                PreExecute(void);
    virtual void                PostExecute(void);
    virtual void                UpdateDataObjectInfo(void);
    virtual void                VerifyInput(void);

    virtual avtContract_p
                                ModifyContract(avtContract_p);
    virtual bool                FilterUnderstandsTransformedRectMesh();
};


#endif
