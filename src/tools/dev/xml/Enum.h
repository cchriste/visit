/*****************************************************************************
*
* Copyright (c) 2000 - 2018, Lawrence Livermore National Security, LLC
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

#ifndef ENUM_H
#define ENUM_H

#include <QTextStream>

#include <vector>

// ****************************************************************************
//  Class:  EnumType
//
//  Purpose:
//    Abstraction for an enumerated type.
//
//  Programmer:  Jeremy Meredith
//  Creation:    August 28, 2001
//
//  Modifications:
//    Jeremy Meredith, Thu Oct 17 15:57:44 PDT 2002
//    Moved the static data to a .C file.
//
//    Eric Brugger, Mon Jul 26 15:00:00 PDT 2004
//    I changed cout to out references in the Print method.
//
//    Brad Whitlock, Thu May  8 11:31:42 PDT 2008
//    Qt 4. Use QTextStream.
//
//    Brad Whitlock, Tue Sep 26 12:12:34 PDT 2017
//    Support adding enum int values.
//
// ****************************************************************************
class EnumType
{
  public:
    static std::vector<EnumType*> enums;
    static EnumType *FindEnum(const QString &s)
    {
        EnumType *e = NULL;
        for (size_t i=0; i<enums.size(); i++)
        {
            if (enums[i]->type == s)
            {
                e = enums[i];
            }
        }
        if (!e)
            throw QString("unknown enum subtype '%1'").arg(s);
        return e;
    }
  public:
    QString         type;
    std::vector<QString> values;
    std::vector<int>     ivalues;
  public:
    EnumType(const QString &s) : type(s), values(), ivalues()
    { 
    }
    void AddValue(const QString &s)
    {
        int n;
        if((n = s.indexOf("=")) != -1)
        {
            values.push_back(s.left(n).simplified());
            bool ok = false;
            int ival = s.mid(n+1).simplified().toInt(&ok);
            ivalues.push_back(ok ? ival : -1);
        }
        else
        {
            values.push_back(s);
            ivalues.push_back(-1);
        }
    }
    const QString& GetValue(size_t index)
    {
        if (index >= values.size())
            throw QString("tried to access out-of-bounds enum type %1").arg(index);
        return values[index];
    }
    const int GetIValue(size_t index)
    {
        if (index >= values.size())
            throw QString("tried to access out-of-bounds enum type %1").arg(index);
        return ivalues[index];
    }
    void Print(QTextStream &out)
    {
        out << "Enum: " << type << endl;
        for (size_t i=0; i<values.size(); i++)
        {
            out << "    " << values[i];
            if(ivalues[i] >= 0)
                out << " = " << ivalues[i];
            out << endl;
        }
    }
};

#endif
