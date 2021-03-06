/*****************************************************************************
*
* Copyright (c) 2000 - 2015, Lawrence Livermore National Security, LLC
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
//                    avtSamplePointsToSamplePointsFilter.C                  //
// ************************************************************************* //

#include <avtSamplePointsToSamplePointsFilter.h>

avtSamplePointsToSamplePointsFilter::avtSamplePointsToSamplePointsFilter()
{
}

avtSamplePointsToSamplePointsFilter::~avtSamplePointsToSamplePointsFilter()
{
}

// ****************************************************************************
//  Method: avtSamplePointsToSamplePointsFilter::PreExecute
//
//  Purpose:
//      Executes before the avtFilter calls Execute.  The sample points need
//      to know how many variables they will be sampling over before they get
//      going.  This is a chance to set that up.
//
//  Programmer: Hank Childs
//  Creation:   November 28, 2001
//
//  Modifications:
//    Jeremy Meredith, Thu Feb 15 11:55:03 EST 2007
//    Call inherited PreExecute before everything else.
//
//    Hank Childs, Thu May 31 13:54:52 PDT 2007
//    Added support for non-scalar variables.
//
//    Hank Childs, Thu Jun 21 11:40:04 PDT 2007
//    Fix stupid bug where variables weren't being set up correctly.
//
// ****************************************************************************

void
avtSamplePointsToSamplePointsFilter::PreExecute(void)
{
    avtSamplePointsToDataObjectFilter::PreExecute();

    avtSamplePoints_p in  = GetTypedInput();
    avtSamplePoints_p out = GetTypedOutput();
    int nvar = in->GetNumberOfRealVariables();
    std::vector<std::string> varnames;
    std::vector<int>         varsizes;
    for (int i = 0 ; i < nvar ; i++)
    {
         varnames.push_back(in->GetVariableName(i));
         varsizes.push_back(in->GetVariableSize(i));
    }

    out->SetNumberOfVariables(varsizes, varnames);
}


