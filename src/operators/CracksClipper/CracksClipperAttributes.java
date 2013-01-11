// ***************************************************************************
//
// Copyright (c) 2000 - 2012, Lawrence Livermore National Security, LLC
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

// ****************************************************************************
// Class: CracksClipperAttributes
//
// Purpose:
//    Attributes for the cracks clipper operator
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

public class CracksClipperAttributes extends AttributeSubject implements Plugin
{
    private static int CracksClipperAttributes_numAdditionalAtts = 8;

    public CracksClipperAttributes()
    {
        super(CracksClipperAttributes_numAdditionalAtts);

        crack1Var = new String("crack1_dir");
        crack2Var = new String("crack2_dir");
        crack3Var = new String("crack3_dir");
        strainVar = new String("void_strain_ten");
        showCrack1 = true;
        showCrack2 = true;
        showCrack3 = true;
        inMassVar = new String("ems");
    }

    public CracksClipperAttributes(int nMoreFields)
    {
        super(CracksClipperAttributes_numAdditionalAtts + nMoreFields);

        crack1Var = new String("crack1_dir");
        crack2Var = new String("crack2_dir");
        crack3Var = new String("crack3_dir");
        strainVar = new String("void_strain_ten");
        showCrack1 = true;
        showCrack2 = true;
        showCrack3 = true;
        inMassVar = new String("ems");
    }

    public CracksClipperAttributes(CracksClipperAttributes obj)
    {
        super(CracksClipperAttributes_numAdditionalAtts);

        crack1Var = new String(obj.crack1Var);
        crack2Var = new String(obj.crack2Var);
        crack3Var = new String(obj.crack3Var);
        strainVar = new String(obj.strainVar);
        showCrack1 = obj.showCrack1;
        showCrack2 = obj.showCrack2;
        showCrack3 = obj.showCrack3;
        inMassVar = new String(obj.inMassVar);

        SelectAll();
    }

    public int Offset()
    {
        return super.Offset() + super.GetNumAdditionalAttributes();
    }

    public int GetNumAdditionalAttributes()
    {
        return CracksClipperAttributes_numAdditionalAtts;
    }

    public boolean equals(CracksClipperAttributes obj)
    {
        // Create the return value
        return ((crack1Var.equals(obj.crack1Var)) &&
                (crack2Var.equals(obj.crack2Var)) &&
                (crack3Var.equals(obj.crack3Var)) &&
                (strainVar.equals(obj.strainVar)) &&
                (showCrack1 == obj.showCrack1) &&
                (showCrack2 == obj.showCrack2) &&
                (showCrack3 == obj.showCrack3) &&
                (inMassVar.equals(obj.inMassVar)));
    }

    public String GetName() { return "CracksClipper"; }
    public String GetVersion() { return "1.0"; }

    // Property setting methods
    public void SetCrack1Var(String crack1Var_)
    {
        crack1Var = crack1Var_;
        Select(0);
    }

    public void SetCrack2Var(String crack2Var_)
    {
        crack2Var = crack2Var_;
        Select(1);
    }

    public void SetCrack3Var(String crack3Var_)
    {
        crack3Var = crack3Var_;
        Select(2);
    }

    public void SetStrainVar(String strainVar_)
    {
        strainVar = strainVar_;
        Select(3);
    }

    public void SetShowCrack1(boolean showCrack1_)
    {
        showCrack1 = showCrack1_;
        Select(4);
    }

    public void SetShowCrack2(boolean showCrack2_)
    {
        showCrack2 = showCrack2_;
        Select(5);
    }

    public void SetShowCrack3(boolean showCrack3_)
    {
        showCrack3 = showCrack3_;
        Select(6);
    }

    public void SetInMassVar(String inMassVar_)
    {
        inMassVar = inMassVar_;
        Select(7);
    }

    // Property getting methods
    public String  GetCrack1Var() { return crack1Var; }
    public String  GetCrack2Var() { return crack2Var; }
    public String  GetCrack3Var() { return crack3Var; }
    public String  GetStrainVar() { return strainVar; }
    public boolean GetShowCrack1() { return showCrack1; }
    public boolean GetShowCrack2() { return showCrack2; }
    public boolean GetShowCrack3() { return showCrack3; }
    public String  GetInMassVar() { return inMassVar; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteString(crack1Var);
        if(WriteSelect(1, buf))
            buf.WriteString(crack2Var);
        if(WriteSelect(2, buf))
            buf.WriteString(crack3Var);
        if(WriteSelect(3, buf))
            buf.WriteString(strainVar);
        if(WriteSelect(4, buf))
            buf.WriteBool(showCrack1);
        if(WriteSelect(5, buf))
            buf.WriteBool(showCrack2);
        if(WriteSelect(6, buf))
            buf.WriteBool(showCrack3);
        if(WriteSelect(7, buf))
            buf.WriteString(inMassVar);
    }

    public void ReadAtts(int index, CommunicationBuffer buf)
    {
        switch(index)
        {
        case 0:
            SetCrack1Var(buf.ReadString());
            break;
        case 1:
            SetCrack2Var(buf.ReadString());
            break;
        case 2:
            SetCrack3Var(buf.ReadString());
            break;
        case 3:
            SetStrainVar(buf.ReadString());
            break;
        case 4:
            SetShowCrack1(buf.ReadBool());
            break;
        case 5:
            SetShowCrack2(buf.ReadBool());
            break;
        case 6:
            SetShowCrack3(buf.ReadBool());
            break;
        case 7:
            SetInMassVar(buf.ReadString());
            break;
        }
    }

    public String toString(String indent)
    {
        String str = new String();
        str = str + stringToString("crack1Var", crack1Var, indent) + "\n";
        str = str + stringToString("crack2Var", crack2Var, indent) + "\n";
        str = str + stringToString("crack3Var", crack3Var, indent) + "\n";
        str = str + stringToString("strainVar", strainVar, indent) + "\n";
        str = str + boolToString("showCrack1", showCrack1, indent) + "\n";
        str = str + boolToString("showCrack2", showCrack2, indent) + "\n";
        str = str + boolToString("showCrack3", showCrack3, indent) + "\n";
        str = str + stringToString("inMassVar", inMassVar, indent) + "\n";
        return str;
    }


    // Attributes
    private String  crack1Var;
    private String  crack2Var;
    private String  crack3Var;
    private String  strainVar;
    private boolean showCrack1;
    private boolean showCrack2;
    private boolean showCrack3;
    private String  inMassVar;
}

