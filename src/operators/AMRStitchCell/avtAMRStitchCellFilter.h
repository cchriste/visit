/*****************************************************************************
*
* Copyright (c) 2000 - 2011, Lawrence Livermore National Security, LLC
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
//  File: avtAMRStitchCellFilter.h
// ************************************************************************* //

#ifndef AVT_AMRStitchCell_FILTER_H
#define AVT_AMRStitchCell_FILTER_H


#include <avtPluginDataTreeIterator.h>
#include <AMRStitchCellAttributes.h>

#include <vtkType.h>

#include <string>
#include <vector>

class vtkDataSet;
class vtkRectilinearGrid;

// ****************************************************************************
//  Class: avtAMRStitchCellFilter
//
//  Purpose:
//      A plugin operator for AMRStitchCell.
//
//  Programmer: ghweber -- generated by xml2avt
//  Creation:   Thu Jul 8 15:14:01 PST 2010
//
// ****************************************************************************

class avtAMRStitchCellFilter :  public virtual avtSIMODataTreeIterator,
                                public virtual avtPluginFilter
{
  public:
                          avtAMRStitchCellFilter();
    virtual              ~avtAMRStitchCellFilter();

    static avtFilter     *Create();

    virtual const char   *GetType(void)  { return "avtAMRStitchCellFilter"; };
    virtual const char   *GetDescription(void)
                              { return "Creating dual grid and stitch cells for AMR data"; };

    virtual void          SetAtts(const AttributeGroup*);
    virtual bool          Equivalent(const AttributeGroup*);

  protected:
    AMRStitchCellAttributes   atts;

    virtual void          PreExecute(void);
    virtual avtContract_p ModifyContract(avtContract_p);
    virtual avtDataTree_p ExecuteDataTree(vtkDataSet *, int, std::string);
    vtkDataSet           *CreateStitchCells(vtkRectilinearGrid*, int, int,
                              int[3], int[3], const std::vector<int>&, const vtkIdType*);
    vtkDataSet           *CreateDualGrid(vtkRectilinearGrid*, int, int, int[3],
                              const vtkIdType*);
    vtkDataArray         *ContractDual(vtkDataArray *coords);

  private:
    int                  topologicalDimension;
    double               domainOrigin[3];
    double               domainBoundingBox[6];
};

#endif
