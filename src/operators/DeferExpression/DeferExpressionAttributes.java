// ***************************************************************************
//
// Copyright (c) 2000 - 2017, Lawrence Livermore National Security, LLC
// Produced at the Lawrence Livermore National Laboratory
// LLNL-CODE-442911
// All rights reserved.
//
// This file is  part of VisIt. For  details, see https://visit.llnl.gov/.  The
// full copyright notice is contained in the file COPYRIGHT located at the root
// of the VisIt distribution or at http://www.llnl.gov/visit/copyright.html.
//
// Redistribution  and  use  in  source  and  binary  forms,  with  or  without
// modification, are permitted provided that the following conditions are met:
//
//  - Redistributions of  source code must  retain the above  copyright notice,
//    this list of conditions and the disclaimer below.
//  - Redistributions in binary form must reproduce the above copyright notice,
//    this  list of  conditions  and  the  disclaimer (as noted below)  in  the
//    documentation and/or other materials provided with the distribution.
//  - Neither the name of  the LLNS/LLNL nor the names of  its contributors may
//    be used to endorse or promote products derived from this software without
//    specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT  HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR  IMPLIED WARRANTIES, INCLUDING,  BUT NOT  LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND  FITNESS FOR A PARTICULAR  PURPOSE
// ARE  DISCLAIMED. IN  NO EVENT  SHALL LAWRENCE  LIVERMORE NATIONAL  SECURITY,
// LLC, THE  U.S.  DEPARTMENT OF  ENERGY  OR  CONTRIBUTORS BE  LIABLE  FOR  ANY
// DIRECT,  INDIRECT,   INCIDENTAL,   SPECIAL,   EXEMPLARY,  OR   CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT  LIMITED TO, PROCUREMENT OF  SUBSTITUTE GOODS OR
// SERVICES; LOSS OF  USE, DATA, OR PROFITS; OR  BUSINESS INTERRUPTION) HOWEVER
// CAUSED  AND  ON  ANY  THEORY  OF  LIABILITY,  WHETHER  IN  CONTRACT,  STRICT
// LIABILITY, OR TORT  (INCLUDING NEGLIGENCE OR OTHERWISE)  ARISING IN ANY  WAY
// OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ***************************************************************************

package llnl.visit.operators;

import llnl.visit.AttributeSubject;
import llnl.visit.CommunicationBuffer;
import llnl.visit.Plugin;
import java.util.Vector;

// ****************************************************************************
// Class: DeferExpressionAttributes
//
// Purpose:
//    Attributes for the DeferExpression operator
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

public class DeferExpressionAttributes extends AttributeSubject implements Plugin
{
    private static int DeferExpressionAttributes_numAdditionalAtts = 1;

    public DeferExpressionAttributes()
    {
        super(DeferExpressionAttributes_numAdditionalAtts);

        exprs = new Vector();
    }

    public DeferExpressionAttributes(int nMoreFields)
    {
        super(DeferExpressionAttributes_numAdditionalAtts + nMoreFields);

        exprs = new Vector();
    }

    public DeferExpressionAttributes(DeferExpressionAttributes obj)
    {
        super(obj);

        int i;

        exprs = new Vector(obj.exprs.size());
        for(i = 0; i < obj.exprs.size(); ++i)
            exprs.addElement(new String((String)obj.exprs.elementAt(i)));


        SelectAll();
    }

    public int Offset()
    {
        return super.Offset() + super.GetNumAdditionalAttributes();
    }

    public int GetNumAdditionalAttributes()
    {
        return DeferExpressionAttributes_numAdditionalAtts;
    }

    public boolean equals(DeferExpressionAttributes obj)
    {
        int i;

        // Compare the elements in the exprs vector.
        boolean exprs_equal = (obj.exprs.size() == exprs.size());
        for(i = 0; (i < exprs.size()) && exprs_equal; ++i)
        {
            // Make references to String from Object.
            String exprs1 = (String)exprs.elementAt(i);
            String exprs2 = (String)obj.exprs.elementAt(i);
            exprs_equal = exprs1.equals(exprs2);
        }
        // Create the return value
        return (exprs_equal);
    }

    public String GetName() { return "DeferExpression"; }
    public String GetVersion() { return "1.0"; }

    // Property setting methods
    public void SetExprs(Vector exprs_)
    {
        exprs = exprs_;
        Select(0);
    }

    // Property getting methods
    public Vector GetExprs() { return exprs; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteStringVector(exprs);
    }

    public void ReadAtts(int index, CommunicationBuffer buf)
    {
        SetExprs(buf.ReadStringVector());
    }

    public String toString(String indent)
    {
        String str = new String();
        str = str + stringVectorToString("exprs", exprs, indent) + "\n";
        return str;
    }


    // Attributes
    private Vector exprs; // vector of String objects
}

