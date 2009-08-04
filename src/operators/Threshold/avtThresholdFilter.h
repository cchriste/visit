/*****************************************************************************
*
* Copyright (c) 2000 - 2009, Lawrence Livermore National Security, LLC
* Produced at the Lawrence Livermore National Laboratory
* LLNL-CODE-400124
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
//  File: avtThresholdFilter.h
// ************************************************************************* //

#ifndef AVT_Threshold_FILTER_H
#define AVT_Threshold_FILTER_H

#include <map>
#include <string>

#include <avtPluginStructuredChunkDataTreeIterator.h>
#include <ThresholdAttributes.h>

#include <avtGhostData.h>

class     vtkDataSet;


// ****************************************************************************
//  Class: avtThresholdFilter
//
//  Purpose:
//      A plugin operator for Threshold.
//
//  Programmer: childs -- generated by xml2info
//  Creation:   Tue Oct 23 16:38:18 PST 2001
//
//  Modifications:
//
//    Hank Childs, Sat Mar 19 10:18:52 PST 2005
//    Add support for structured chunking.
//
//    Hank Childs, Sun Mar 27 11:34:04 PST 2005
//    Inherit from new base type that supports structured chunking.
//
//    Hank Childs, Tue Sep 13 09:07:05 PDT 2005
//    Add support for PointsOnly mode.
//
//    Mark Blair, Tue Mar  7 13:25:00 PST 2006
//    Reworked to support multi-variable thresholding.
//
//    Mark Blair, Tue Aug  8 17:47:00 PDT 2006
//    Now accommodates an empty list of threshold variables; does pass-through.
//
//    Mark Blair, Wed Oct  4 17:45:48 PDT 2006
//    Keeps track of the "active variable".
//
//    Sean Ahern (Markus Glatter), Tue Sep 25 10:14:26 EDT 2007
//    Keep around a mapping of data selections and variable names.
//
//    Hank Childs, Mon Feb 23 11:54:22 PST 2009
//    Added CreateNamedSelection.
//
// ****************************************************************************

class avtThresholdFilter : public avtPluginStructuredChunkDataTreeIterator
{
  public:
                          avtThresholdFilter();
    virtual              ~avtThresholdFilter() {;};

    static avtFilter     *Create();

    virtual const char   *GetType(void)  { return "avtThresholdFilter"; };
    virtual const char   *GetDescription(void) { return "Thresholding"; };

    virtual void          SetAtts(const AttributeGroup*);
    virtual bool          Equivalent(const AttributeGroup*);

    virtual avtNamedSelection *
                          CreateNamedSelection(avtContract_p c,
                                               const std::string &);

  protected:
    ThresholdAttributes   atts;
    std::string           activeVarName;
    std::map<std::string,int> selIDs;

    virtual avtContract_p ModifyContract(avtContract_p);

    virtual vtkDataSet   *ProcessOneChunk(vtkDataSet *, int, std::string,bool);
    virtual void          GetAssignments(vtkDataSet *, const int *,
                      std::vector<avtStructuredMeshChunker::ZoneDesignation>&);
    vtkDataSet           *ThresholdToPointMesh(vtkDataSet *in_ds);

    virtual void          UpdateDataObjectInfo(void);
    virtual void          PreExecute(void);
};


#endif


