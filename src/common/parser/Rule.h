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

#ifndef RULE_H
#define RULE_H
#include <parser_exports.h>

#include "Symbol.h"
#include "Sequence.h"

// ****************************************************************************
//  Class:  Rule
//
//  Purpose:
//    A production rule in a grammar.  It has a nonterminal on the left
//    and a sequence on the right.
//
//  Programmer:  Jeremy Meredith
//  Creation:    April  5, 2002
//
// ****************************************************************************
class PARSER_API Rule
{
  public:
    Rule();
    Rule(int, const Symbol&);

    Rule &operator>>(const Sequence&);

    bool  operator==(const Rule &r) const {return lhs == r.lhs && rhs == r.rhs; }
    void  operator=(const Rule &r) { index = r.index; lhs = r.lhs; rhs = r.rhs; }

    void  Print(ostream &o, int pos = -1) const;
    void  PrintNoColor(ostream &o, int pos = -1) const;
    bool  IsNullable(const std::vector<const Rule*> &rules) const;

    const Symbol   *GetLHS()        const { return lhs; }
    const Sequence &GetRHS()        const { return rhs; }
    int             GetID()         const { return id; }
    int             GetPrec()       const { return prec; }
    int             GetIndex()      const { return index; }
    void            SetPrec(int p)        { prec = p; }
    void            SetIndex(int i)       { index = i; }
  private:
    const Symbol  *lhs;
    Sequence       rhs;
    int            id;     // client-specified identifier, used by same client
    int            prec;   // precision of rule
    int            index;  // index among all rules
};

PARSER_API ostream &operator<<(ostream &o, const Rule &r);


#endif
