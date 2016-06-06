// ***************************************************************************
//
// Copyright (c) 2000 - 2015, Lawrence Livermore National Security, LLC
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
import java.lang.Integer;
import java.util.Vector;

// ****************************************************************************
// Class: AxisAlignedSlice4DAttributes
//
// Purpose:
//    Attributes for AxisAlignedSlice4D
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

public class AxisAlignedSlice4DAttributes extends AttributeSubject implements Plugin
{
    private static int AxisAlignedSlice4DAttributes_numAdditionalAtts = 4;

    public AxisAlignedSlice4DAttributes()
    {
        super(AxisAlignedSlice4DAttributes_numAdditionalAtts);

        I = new Vector();
        J = new Vector();
        K = new Vector();
        L = new Vector();
    }

    public AxisAlignedSlice4DAttributes(int nMoreFields)
    {
        super(AxisAlignedSlice4DAttributes_numAdditionalAtts + nMoreFields);

        I = new Vector();
        J = new Vector();
        K = new Vector();
        L = new Vector();
    }

    public AxisAlignedSlice4DAttributes(AxisAlignedSlice4DAttributes obj)
    {
        super(AxisAlignedSlice4DAttributes_numAdditionalAtts);

        int i;

        I = new Vector();
        for(i = 0; i < obj.I.size(); ++i)
        {
            Integer iv = (Integer)obj.I.elementAt(i);
            I.addElement(new Integer(iv.intValue()));
        }
        J = new Vector();
        for(i = 0; i < obj.J.size(); ++i)
        {
            Integer iv = (Integer)obj.J.elementAt(i);
            J.addElement(new Integer(iv.intValue()));
        }
        K = new Vector();
        for(i = 0; i < obj.K.size(); ++i)
        {
            Integer iv = (Integer)obj.K.elementAt(i);
            K.addElement(new Integer(iv.intValue()));
        }
        L = new Vector();
        for(i = 0; i < obj.L.size(); ++i)
        {
            Integer iv = (Integer)obj.L.elementAt(i);
            L.addElement(new Integer(iv.intValue()));
        }

        SelectAll();
    }

    public int Offset()
    {
        return super.Offset() + super.GetNumAdditionalAttributes();
    }

    public int GetNumAdditionalAttributes()
    {
        return AxisAlignedSlice4DAttributes_numAdditionalAtts;
    }

    public boolean equals(AxisAlignedSlice4DAttributes obj)
    {
        int i;

        // Compare the elements in the I vector.
        boolean I_equal = (obj.I.size() == I.size());
        for(i = 0; (i < I.size()) && I_equal; ++i)
        {
            // Make references to Integer from Object.
            Integer I1 = (Integer)I.elementAt(i);
            Integer I2 = (Integer)obj.I.elementAt(i);
            I_equal = I1.equals(I2);
        }
        // Compare the elements in the J vector.
        boolean J_equal = (obj.J.size() == J.size());
        for(i = 0; (i < J.size()) && J_equal; ++i)
        {
            // Make references to Integer from Object.
            Integer J1 = (Integer)J.elementAt(i);
            Integer J2 = (Integer)obj.J.elementAt(i);
            J_equal = J1.equals(J2);
        }
        // Compare the elements in the K vector.
        boolean K_equal = (obj.K.size() == K.size());
        for(i = 0; (i < K.size()) && K_equal; ++i)
        {
            // Make references to Integer from Object.
            Integer K1 = (Integer)K.elementAt(i);
            Integer K2 = (Integer)obj.K.elementAt(i);
            K_equal = K1.equals(K2);
        }
        // Compare the elements in the L vector.
        boolean L_equal = (obj.L.size() == L.size());
        for(i = 0; (i < L.size()) && L_equal; ++i)
        {
            // Make references to Integer from Object.
            Integer L1 = (Integer)L.elementAt(i);
            Integer L2 = (Integer)obj.L.elementAt(i);
            L_equal = L1.equals(L2);
        }
        // Create the return value
        return (I_equal &&
                J_equal &&
                K_equal &&
                L_equal);
    }

    public String GetName() { return "AxisAlignedSlice4D"; }
    public String GetVersion() { return "1.0"; }

    // Property setting methods
    public void SetI(Vector I_)
    {
        I = I_;
        Select(0);
    }

    public void SetJ(Vector J_)
    {
        J = J_;
        Select(1);
    }

    public void SetK(Vector K_)
    {
        K = K_;
        Select(2);
    }

    public void SetL(Vector L_)
    {
        L = L_;
        Select(3);
    }

    // Property getting methods
    public Vector GetI() { return I; }
    public Vector GetJ() { return J; }
    public Vector GetK() { return K; }
    public Vector GetL() { return L; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteIntVector(I);
        if(WriteSelect(1, buf))
            buf.WriteIntVector(J);
        if(WriteSelect(2, buf))
            buf.WriteIntVector(K);
        if(WriteSelect(3, buf))
            buf.WriteIntVector(L);
    }

    public void ReadAtts(int index, CommunicationBuffer buf)
    {
        switch(index)
        {
        case 0:
            SetI(buf.ReadIntVector());
            break;
        case 1:
            SetJ(buf.ReadIntVector());
            break;
        case 2:
            SetK(buf.ReadIntVector());
            break;
        case 3:
            SetL(buf.ReadIntVector());
            break;
        }
    }

    public String toString(String indent)
    {
        String str = new String();
        str = str + intVectorToString("I", I, indent) + "\n";
        str = str + intVectorToString("J", J, indent) + "\n";
        str = str + intVectorToString("K", K, indent) + "\n";
        str = str + intVectorToString("L", L, indent) + "\n";
        return str;
    }


    // Attributes
    private Vector I; // vector of Integer objects
    private Vector J; // vector of Integer objects
    private Vector K; // vector of Integer objects
    private Vector L; // vector of Integer objects
}

