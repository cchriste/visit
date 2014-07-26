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
//                              avtSILGenerator.h                            //
// ************************************************************************* //

#ifndef AVT_SIL_GENERATOR_H
#define AVT_SIL_GENERATOR_H

#include <database_exports.h>

#include <string>
#include <vector>

#include <avtDatabaseMetaData.h>
#include <avtSILCollection.h>


class     avtSIL;
class     NameschemeAttributes;


// ****************************************************************************
//  Class: avtSILGenerator
//
//  Purpose:
//      Takes avtDatabaseMetaData and constructs a reasonable SIL from it.
//
//  Notes:      Created from pre-existing routines in the avtGenericDatabase.
//
//  Programmer: Hank Childs
//  Creation:   September 6, 2002
//
//  Modifications:
//    Jeremy Meredith, Wed Aug 24 10:34:13 PDT 2005
//    Added origin for the groups.
//
//    Jeremy Meredith, Tue Aug 29 16:52:36 EDT 2006
//    Added enumerated scalar subset support.
//
//    Brad Whitlock, Thu Mar 8 09:46:17 PDT 2007
//    Changed AddSpecies method.
//
//    Dave Bremer, Tue Apr  1 16:43:16 PDT 2008
//    Added a flag to AddSubsets to control the use of SIL Arrays
//
//    Mark C. Miller, Wed Aug 26 11:08:47 PDT 2009
//    Removed custom SIL stuff.
//
//    Hank Childs, Mon Dec  7 14:05:00 PST 2009
//    Added various optimizations that enable the SIL object to get smaller.
//
//    Kathleen Biagas, Thu Aug 22 09:58:50 PDT 2013
//    Add groupNames argment to AddGroups.
//
// ****************************************************************************

class DATABASE_API avtSILGenerator
{
  public:
    void                CreateSIL(avtDatabaseMetaData *, avtSIL *);

  protected:
    void                AddSubsets(avtSIL *, int, int, int, std::vector<int> &,
                                   const std::string &, const std::string &,
                                   const std::vector<std::string> &,
                                   const NameschemeAttributes &,
                                   bool useSILArrays,
                                   SILCategoryRole cat = SIL_DOMAIN,
                                   bool onlyCreateSets = false);
    std::vector<int>    AddGroups(avtSIL *, int, int, int, const std::string &,
                                  const std::string &,
                                  const std::vector< std::string > &);
    void                AddGroupCollections(avtSIL *, int, int,
                                  const std::vector<int> &,
                                  const std::vector<int>&, 
                                  const std::vector<int>&, 
                                  const std::string &, const std::vector<int> &);
    void                AddMaterials(avtSIL *, int, const std::string &,
                                     const std::vector<std::string> &,
                                     std::vector<int> &, int);
    void                AddSpecies(avtSIL *, int,
                                   const std::vector<std::string>&,
                                   const std::string &,
                                   const avtSpeciesMetaData *,
                                   int);
    void                AddMaterialSubsets(avtSIL *, const std::vector<int> &,
                                           int, int , const std::vector<int> &,
                                           const std::string &,
                                           const std::vector<std::string> &,
                                           const std::vector<std::string> &);
    void                AddEnumScalars(avtSIL *, int,
                                       const avtScalarMetaData *);
};


#endif
 

