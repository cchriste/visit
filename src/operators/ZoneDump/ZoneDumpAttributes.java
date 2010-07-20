// ***************************************************************************
//
// Copyright (c) 2000 - 2010, Lawrence Livermore National Security, LLC
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
// Class: ZoneDumpAttributes
//
// Purpose:
//    Zone Dump Control
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

public class ZoneDumpAttributes extends AttributeSubject implements Plugin
{
    private static int numAdditionalAttributes = 5;

    public ZoneDumpAttributes()
    {
        super(numAdditionalAttributes);

        variable = new String("");
        lowerBound = -1e+37;
        upperBound = 1e+37;
        outputFile = new String("visit_zonedump.zod");
        enabled = false;
    }

    public ZoneDumpAttributes(int nMoreFields)
    {
        super(numAdditionalAttributes + nMoreFields);

        variable = new String("");
        lowerBound = -1e+37;
        upperBound = 1e+37;
        outputFile = new String("visit_zonedump.zod");
        enabled = false;
    }

    public ZoneDumpAttributes(ZoneDumpAttributes obj)
    {
        super(numAdditionalAttributes);

        variable = new String(obj.variable);
        lowerBound = obj.lowerBound;
        upperBound = obj.upperBound;
        outputFile = new String(obj.outputFile);
        enabled = obj.enabled;

        SelectAll();
    }

    public int Offset()
    {
        return super.Offset() + super.GetNumAdditionalAttributes();
    }

    public int GetNumAdditionalAttributes()
    {
        return numAdditionalAttributes;
    }

    public boolean equals(ZoneDumpAttributes obj)
    {
        // Create the return value
        return ((variable.equals(obj.variable)) &&
                (lowerBound == obj.lowerBound) &&
                (upperBound == obj.upperBound) &&
                (outputFile.equals(obj.outputFile)) &&
                (enabled == obj.enabled));
    }

    public String GetName() { return "ZoneDump"; }
    public String GetVersion() { return "1.0"; }

    // Property setting methods
    public void SetVariable(String variable_)
    {
        variable = variable_;
        Select(0);
    }

    public void SetLowerBound(double lowerBound_)
    {
        lowerBound = lowerBound_;
        Select(1);
    }

    public void SetUpperBound(double upperBound_)
    {
        upperBound = upperBound_;
        Select(2);
    }

    public void SetOutputFile(String outputFile_)
    {
        outputFile = outputFile_;
        Select(3);
    }

    public void SetEnabled(boolean enabled_)
    {
        enabled = enabled_;
        Select(4);
    }

    // Property getting methods
    public String  GetVariable() { return variable; }
    public double  GetLowerBound() { return lowerBound; }
    public double  GetUpperBound() { return upperBound; }
    public String  GetOutputFile() { return outputFile; }
    public boolean GetEnabled() { return enabled; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteString(variable);
        if(WriteSelect(1, buf))
            buf.WriteDouble(lowerBound);
        if(WriteSelect(2, buf))
            buf.WriteDouble(upperBound);
        if(WriteSelect(3, buf))
            buf.WriteString(outputFile);
        if(WriteSelect(4, buf))
            buf.WriteBool(enabled);
    }

    public void ReadAtts(int index, CommunicationBuffer buf)
    {
        switch(index)
        {
        case 0:
            SetVariable(buf.ReadString());
            break;
        case 1:
            SetLowerBound(buf.ReadDouble());
            break;
        case 2:
            SetUpperBound(buf.ReadDouble());
            break;
        case 3:
            SetOutputFile(buf.ReadString());
            break;
        case 4:
            SetEnabled(buf.ReadBool());
            break;
        }
    }

    public String toString(String indent)
    {
        String str = new String();
        str = str + stringToString("variable", variable, indent) + "\n";
        str = str + doubleToString("lowerBound", lowerBound, indent) + "\n";
        str = str + doubleToString("upperBound", upperBound, indent) + "\n";
        str = str + stringToString("outputFile", outputFile, indent) + "\n";
        str = str + boolToString("enabled", enabled, indent) + "\n";
        return str;
    }


    // Attributes
    private String  variable;
    private double  lowerBound;
    private double  upperBound;
    private String  outputFile;
    private boolean enabled;
}

