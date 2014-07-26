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
//                       avtIdentifierSelection.h                            //
// ************************************************************************* //

#ifndef AVT_IDENTIFIER_SELECTION_H
#define AVT_IDENTIFIER_SELECTION_H 

#include <pipeline_exports.h>

#include <ref_ptr.h>

#include <avtDataSelection.h>

#include <vector>
#include <string>

// ****************************************************************************
//  Class: avtIdentifierSelection
//
//  Purpose: 
//      Specify data selection using identifiers.
//
//  Programmer: Hank Childs
//  Creation:   March 5, 2008
//
//  Modifications:
//
//    Hank Childs, Thu Mar  6 09:05:59 PST 2008
//    Add Destruct method.
//
//    Hank Childs, Tue Dec 20 14:43:08 PST 2011
//    Add method DescriptionString.
//
//    Brad Whitlock, Thu Mar 15 14:13:59 PDT 2012
//    Added idVar.
//
// ****************************************************************************

class PIPELINE_API avtIdentifierSelection : public avtDataSelection 
{
  public:
                            avtIdentifierSelection();
    virtual                ~avtIdentifierSelection();

    static void             Destruct(void *);

    virtual const char *    GetType() const
                                { return "Identifier Data Selection"; }; 
    virtual std::string     DescriptionString(void);

    void                    SetIdentifiers(const std::vector<double> &a)
                                { ids = a; };
    const std::vector<double> &GetIdentifiers(void) { return ids; };

    void                    SetIdVariable(const std::string &id) {idVar = id; }
    const std::string      &GetIdVariable() const {return idVar; }

    bool                    operator==(const avtIdentifierSelection &) const;

  private:
    std::vector<double>     ids;
    std::string             idVar;

    // These methods are defined to prevent accidental use of bitwise copy
    // implementations.  If you want to re-define them to do something
    // meaningful, that's fine.
                    avtIdentifierSelection(const avtIdentifierSelection &) {;};
    avtIdentifierSelection &operator=(const avtIdentifierSelection &) 
                                                            { return *this; };
};

typedef ref_ptr<avtIdentifierSelection> avtIdentifierSelection_p;


#endif


