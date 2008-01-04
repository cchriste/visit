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
//  File: avtIsovolumeFilter.h
// ************************************************************************* //

#ifndef AVT_Isovolume_FILTER_H
#define AVT_Isovolume_FILTER_H


#include <avtPluginStreamer.h>
#include <IsovolumeAttributes.h>


class vtkDataSet;

// ****************************************************************************
//  Class: avtIsovolumeFilter
//
//  Purpose:
//      A plugin operator for Isovolume.
//
//  Programmer: meredith -- generated by xml2info
//  Creation:   Fri Jan 30 14:50:21 PST 2004
//
//  Modifications:
//    Jeremy Meredith, Mon Feb 16 19:12:11 PST 2004
//    Added RefashionDataObjectInfo.  This was needed for correct support
//    on various mesh types.
//
//    Jeremy Meredith, Thu May  6 11:37:47 PDT 2004
//    Split some code from ExecuteData into a new function to avoid
//    code duplication.
//
//    Hank Childs, Wed Aug 11 09:23:17 PDT 2004
//    Added PerformRestriction.
//
// ****************************************************************************

class avtIsovolumeFilter : public avtPluginStreamer
{
  public:
                         avtIsovolumeFilter();
    virtual             ~avtIsovolumeFilter();

    static avtFilter    *Create();

    virtual const char  *GetType(void)  { return "avtIsovolumeFilter"; };
    virtual const char  *GetDescription(void)
                             { return "Isovolume"; };

    virtual void         SetAtts(const AttributeGroup*);
    virtual bool         Equivalent(const AttributeGroup*);

  protected:
    IsovolumeAttributes   atts;

    virtual vtkDataSet   *ExecuteData(vtkDataSet *, int, std::string);
    virtual void          RefashionDataObjectInfo(void);

    virtual avtPipelineSpecification_p
                          PerformRestriction(avtPipelineSpecification_p);

  private:
    virtual vtkDataSet   *ExecuteSingleClip(vtkDataSet *, float, bool);
};


#endif
