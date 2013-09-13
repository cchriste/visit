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
//                        avtDataRangeSelection.h                           //
// ************************************************************************* //

#ifndef AVT_DATA_RANGE_SELECTION_H
#define AVT_DATA_RANGE_SELECTION_H
#include <float.h>
#include <string>
 
#include <pipeline_exports.h>

#include <ref_ptr.h>

#include <avtDataSelection.h>

// ****************************************************************************
//  Class: avtDataRangeSelection
//
//  Purpose: Specify a data selection by a scalar range for a named variable.
// 
//  The default is a range from -FLT_MAX to +FLT_MAX for variable "default"
//
//  Programmer: Markus Glatter
//  Creation:   July 27, 2007
//
//  Modifications:
//
//    Hank Childs, Tue Dec 18 10:04:43 PST 2007
//    Define private copy constructor and assignment operator to prevent
//    accidental use of default, bitwise copy implementations.
//
//    Hank Childs, Tue Dec 20 14:43:08 PST 2011
//    Add method DescriptionString.
//
// ****************************************************************************

class PIPELINE_API avtDataRangeSelection : public avtDataSelection 
{
  public:
                            avtDataRangeSelection();
                            avtDataRangeSelection(const std::string _var, 
                                                  const double _min, 
                                                  const double _max);
    virtual                ~avtDataRangeSelection();

    virtual const char *    GetType() const
                                { return "Data Range Selection"; }; 
    virtual std::string     DescriptionString(void);

    void                    SetVariable(const std::string _var)
                                { var = _var; };
    void                    SetMin(const double _min)
                                { min = _min; };
    void                    SetMax(const double _max)
                                { max = _max; };

    std::string             GetVariable() const
                                { return var; };
    double                  GetMin() const
                                { return min; };
    double                  GetMax() const
                                { return max; };

    bool                    operator==(const avtDataRangeSelection &s) const;

  private:
    std::string var;
    double min;
    double max;

    // These methods are defined to prevent accidental use of bitwise copy
    // implementations.  If you want to re-define them to do something
    // meaningful, that's fine.
                            avtDataRangeSelection(const avtDataRangeSelection&)
                                                                           {;};
    avtDataRangeSelection  &operator=(const avtDataRangeSelection &) 
                                                             { return *this; };
};

typedef ref_ptr<avtDataRangeSelection> avtDataRangeSelection_p;


#endif
