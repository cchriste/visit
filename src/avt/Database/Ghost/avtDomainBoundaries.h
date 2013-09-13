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
//                            avtDomainBoundaries.h                          //
// ************************************************************************* //

#ifndef AVT_DOMAIN_BOUNDARIES_H
#define AVT_DOMAIN_BOUNDARIES_H

#include <database_exports.h>

#include <avtGhostData.h>

#include <vector>

class vtkDataSet;
class vtkDataArray;
class avtMixedVariable;
class avtMaterial;

// ****************************************************************************
//  Class:  avtDomainBoundaries
//
//  Purpose:
//    Encapsulate domain boundary information.
//
//  Programmer:  Jeremy Meredith
//  Creation:    October 25, 2001
//
//  Modifications:
//    Jeremy Meredith, Thu Dec 13 11:47:06 PST 2001
//    Added mats to the exchange mixvars call.
//
//    Kathleen Bonnell, Fri Feb  8 11:03:49 PST 2002
//    vtkScalars and vtkVectors have been deprecated in VTK 4.0, 
//    use vtkDataArray instead.
//
//    Kathleen Bonnell, Mon May 20 13:40:17 PDT 2002
//    Made ExhangeVector into two methods to handle different underlying
//    data types (int, float).
//
//    Hank Childs, Sat Aug 14 06:41:00 PDT 2004
//    Added ghost nodes.
//
//    Hank Childs, Sun Feb 27 12:00:12 PST 2005
//    Added pure virtual methods RequiresCommunication.  Added "allDomains"
//    argument to CreateGhostNodes.
//
//    Hank Childs, Mon Jun 27 16:28:22 PDT 2005
//    Added virtual method ResetCachedMembers.
//
//    Hank Childs, Thu Jan 26 10:04:34 PST 2006
//    Add virtual method "CreatesRobustGhostNodes".
//
//    Hank Childs, Thu Feb 14 17:12:38 PST 2008
//    Add virtual method "CanOnlyCreateGhostNodes".
//
//    Brad Whitlock, Sun Apr 22 10:32:28 PDT 2012
//    Added ExchangeDoubleVector.
//
// ****************************************************************************

class DATABASE_API avtDomainBoundaries
{
  public:
                 avtDomainBoundaries();
    virtual      ~avtDomainBoundaries();

    virtual std::vector<vtkDataSet*>       ExchangeMesh(std::vector<int>       domainNum,
                                               std::vector<vtkDataSet*>   meshes)  =0;

    virtual std::vector<vtkDataArray*>     ExchangeScalar(std::vector<int>     domainNum,
                                               bool                  isPointData,
                                               std::vector<vtkDataArray*> scalars) =0;

    virtual std::vector<vtkDataArray*>     ExchangeFloatVector(std::vector<int> domainNum,
                                               bool                   isPointData,
                                               std::vector<vtkDataArray*>  vectors) =0;

    virtual std::vector<vtkDataArray*>     ExchangeDoubleVector(std::vector<int> domainNum,
                                               bool                   isPointData,
                                               std::vector<vtkDataArray*>  vectors) =0;

    virtual std::vector<vtkDataArray*>     ExchangeIntVector(std::vector<int>  domainNum,
                                               bool                  isPointData,
                                               std::vector<vtkDataArray*> vectors) =0;

    virtual std::vector<avtMaterial*>      ExchangeMaterial(std::vector<int>   domainNum,
                                              std::vector<avtMaterial*>   mats)    =0;

    virtual std::vector<avtMixedVariable*> ExchangeMixVar(std::vector<int>     domainNum,
                                        const std::vector<avtMaterial*>   mats,
                                        std::vector<avtMixedVariable*>    mixvars) =0;
    virtual void                      CreateGhostNodes(std::vector<int>   domainNum,
                                               std::vector<vtkDataSet*>   meshes,
                                               std::vector<int> &)  =0;
    virtual bool                      CreatesRobustGhostNodes(void) 
                                                              { return true; };
    virtual bool                      CanOnlyCreateGhostNodes(void) 
                                                              { return false; };
    virtual bool                      RequiresCommunication(avtGhostDataType) = 0;
    virtual bool                      ConfirmMesh(std::vector<int>      domainNum,
                                               std::vector<vtkDataSet*> meshes)  =0;
    virtual void                      ResetCachedMembers(void) {;};
};

#endif
