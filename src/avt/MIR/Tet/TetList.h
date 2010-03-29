/*****************************************************************************
*
* Copyright (c) 2000 - 2010, Lawrence Livermore National Security, LLC
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

#ifndef TETLIST_H
#define TETLIST_H
#include <mir_exports.h>

#include "Tet.h"
#include "VisItArray.h"
#include "mat3d_tet.h"

// ****************************************************************************
// ****************************************************************************
//                             class TetList
// ****************************************************************************
// ****************************************************************************

// ****************************************************************************
//  Class:  TetList
//
//  Purpose:
//    Encapsulation of a list of Tets.
//
//  Programmer:  Jeremy Meredith
//  Creation:    February 14, 2001
//
//  Modifications:
//    Jeremy Meredith, Thu Feb 15 14:09:13 PST 2001
//    Made it use a vector, since there is no current need to pop off the
//    front of the list.
//
//    Jeremy Meredith, Thu May 31 17:13:39 PDT 2001
//    Added the Clear method.
//
//    Jeremy Meredith, Mon Sep 17 15:33:35 PDT 2001
//    Made it use an Array instead of a vector, for faster push_backs.
//
//    Jeremy Meredith, Wed Dec 11 10:10:31 PST 2002
//    Added a "forced material" where if it is >=0, any added tet will have
//    the forced material instead of the normal requested material.
//
//    Mark C. Miller, Thu Feb  9 21:06:10 PST 2006
//    Renamed Array class to VisItArray to avoid name collisions with
//    third-party libs
//
// ****************************************************************************
class MIR_API TetList
{
  private:
    VisItArray<Tet> list;
    void operator=(const TetList &rhs) { };
  public:
    void       Clear()                 {list.clear();}
    bool       Empty()                 {return list.empty();}
    const int &Size()            const {return list.size();}
    Tet       &operator[](const int i) {return list[i];}
    void Add(const Tet&, int);
    void AddTet(int, int,   const Tet::Node&,const Tet::Node&,const Tet::Node&,const Tet::Node&, int);
    static void Swap(TetList &a, TetList &b) {VisItArray<Tet>::Swap(a.list, b.list);}
};

// ****************************************************************************
//  Method:  TetList::Add
//
//  Purpose:
//    Add a tet directly to the list.
//
//  Arguments:
//    t      : the Tet
//
//  Programmer:  Jeremy Meredith
//  Creation:    February 14, 2001
//
// ****************************************************************************
void
TetList::Add(const Tet &t, int forcedMat)
{
    list.push_back(t);

    if (forcedMat >= 0)
        list[list.size()-1].mat = forcedMat;
}

// ****************************************************************************
//  Method:  TetList::AddTet
//
//  Purpose:
//    Construct a tet from its nodes and add it to the list.
//
//  Arguments:
//    id     : the cell id
//    mat    : the mat number
//    n0..n3 : the nodes of the tet
//
//  Programmer:  Jeremy Meredith
//  Creation:    February 14, 2001
//
// ****************************************************************************
void
TetList::AddTet(int id,int mat,
                const Tet::Node &n0, const Tet::Node &n1, const Tet::Node &n2, const Tet::Node &n3,
                int forcedMat)
{
    if (forcedMat >= 0)
        mat = forcedMat;

    list.push_back(Tet(id, n0, n1, n2, n3, mat));
}

#endif
