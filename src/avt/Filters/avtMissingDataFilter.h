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
#ifndef AVT_MISSING_DATA_FILTER_H
#define AVT_MISSING_DATA_FILTER_H

#include <filters_exports.h>

#include <avtDatabaseMetaData.h>
#include <avtDataTreeIterator.h>


// ****************************************************************************
// Class: avtMissingDataFilter
//
// Purpose:
//   This filter generates/removes missing data from its inputs.
//
// Notes:      
//
// Programmer: Brad Whitlock
// Creation:   Tue Jan 10 09:36:35 PST 2012
//
// Modifications:
//
//   Dave Pugmire, Thu Mar 22 13:06:30 EDT 2012
//   Added canDoCollectiveCommunication flag to detect and handle when we
//   are streaming.
//   
//   Eric Brugger, Mon Jul 21 14:40:46 PDT 2014
//   Modified the class to work with avtDataRepresentation.
//
// ****************************************************************************

class AVTFILTERS_API avtMissingDataFilter : public avtDataTreeIterator
{
public:
    avtMissingDataFilter();
    virtual ~avtMissingDataFilter();

    virtual const char                 *GetType(void);
    virtual const char                 *GetDescription(void);

    void SetMetaData(const avtDatabaseMetaData *md);

    void SetGenerateMode(bool);
    void SetRemoveMode(bool);

protected:
    virtual void          PreExecute(void);
    virtual avtDataRepresentation *ExecuteData(avtDataRepresentation *);
    virtual void          PostExecute(void);

    virtual avtContract_p ModifyContract(avtContract_p);
    virtual bool          FilterUnderstandsTransformedRectMesh();

    stringVector          MissingDataVariables(avtDataRequest_p spec, 
                              const avtDatabaseMetaData *md) const;
    avtCentering          MissingDataCentering(const stringVector &vars) const;
    vtkDataArray         *MissingDataBuildMask(vtkDataSet *in_ds,
                              avtDataRequest_p spec, 
                              const avtDatabaseMetaData *md, 
                              bool &missing, avtCentering &cent) const;
    bool                  TagMissingData(vtkDataSet *in_ds, 
                              vtkDataArray *missingData, 
                              const stringVector &varsMissingData, 
                              avtCentering centering) const;

    bool                removedData;
    bool                generateMode;
    bool                removeMode;
    bool                canDoCollectiveCommunication;
    avtContract_p       contract;
    avtDatabaseMetaData metadata;
};

#endif
