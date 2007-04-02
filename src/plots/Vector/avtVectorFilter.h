/*****************************************************************************
*
* Copyright (c) 2000 - 2006, The Regents of the University of California
* Produced at the Lawrence Livermore National Laboratory
* All rights reserved.
*
* This file is part of VisIt. For details, see http://www.llnl.gov/visit/. The
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
*    documentation and/or materials provided with the distribution.
*  - Neither the name of the UC/LLNL nor  the names of its contributors may be
*    used to  endorse or  promote products derived from  this software without
*    specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT  HOLDERS AND CONTRIBUTORS "AS IS"
* AND ANY EXPRESS OR  IMPLIED WARRANTIES, INCLUDING,  BUT NOT  LIMITED TO, THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND  FITNESS FOR A PARTICULAR  PURPOSE
* ARE  DISCLAIMED.  IN  NO  EVENT  SHALL  THE  REGENTS  OF  THE  UNIVERSITY OF
* CALIFORNIA, THE U.S.  DEPARTMENT  OF  ENERGY OR CONTRIBUTORS BE  LIABLE  FOR
* ANY  DIRECT,  INDIRECT,  INCIDENTAL,  SPECIAL,  EXEMPLARY,  OR CONSEQUENTIAL
* DAMAGES (INCLUDING, BUT NOT  LIMITED TO, PROCUREMENT OF  SUBSTITUTE GOODS OR
* SERVICES; LOSS OF  USE, DATA, OR PROFITS; OR  BUSINESS INTERRUPTION) HOWEVER
* CAUSED  AND  ON  ANY  THEORY  OF  LIABILITY,  WHETHER  IN  CONTRACT,  STRICT
* LIABILITY, OR TORT  (INCLUDING NEGLIGENCE OR OTHERWISE)  ARISING IN ANY  WAY
* OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
* DAMAGE.
*
*****************************************************************************/

// ************************************************************************* //
//                             avtVectorFilter.h                             //
// ************************************************************************* //

#ifndef AVT_VECTOR_FILTER_H
#define AVT_VECTOR_FILTER_H


#include <avtStreamer.h>

class  vtkVectorReduceFilter;
class  vtkVertexFilter;


// ****************************************************************************
//  Class: avtVectorFilter
//
//  Purpose:
//      A filter that takes in vector data and creates vector glyphs as poly
//      data.
//
//  Programmer: Hank Childs
//  Creation:   March 21, 2001
//
//  Modifications:
//
//    Kathleen Bonnell, Tue Apr 10 11:51:18 PDT 2001
//    Renamed ExecuteDomain as ExecuteData.
//
//    Hank Childs, Thu Aug 30 17:30:48 PDT 2001
//    Added the vertex filter.
//
//    Kathleen Bonnell, Mon Aug  9 14:27:08 PDT 2004 
//    Added magVarName, SetMagVarName and PerformRestriction.
//
//    Kathleen Bonnell, Tue Oct 12 16:18:37 PDT 2004 
//    Added keepNodeZone.
//
//    Hank Childs, Fri Mar 11 15:00:05 PST 2005
//    Instantiate VTK filters on the fly.
//
// ****************************************************************************

class avtVectorFilter : public avtStreamer
{
  public:
                              avtVectorFilter(bool, int);
    virtual                  ~avtVectorFilter();

    virtual const char       *GetType(void)   { return "avtVectorFilter"; };
    virtual const char       *GetDescription(void)
                                  { return "Creating vectors"; };

    bool                      Equivalent(bool, int);

    void                      SetStride(int);
    void                      SetNVectors(int);
    void                      SetMagVarName(const std::string &);

  protected:
    bool                      useStride;
    int                       stride;
    int                       nVectors;
    std::string               magVarName;
    bool                      keepNodeZone;

    virtual vtkDataSet       *ExecuteData(vtkDataSet *, int, std::string);
    virtual void              RefashionDataObjectInfo(void);
    virtual avtPipelineSpecification_p
                              PerformRestriction(avtPipelineSpecification_p);
};


#endif


