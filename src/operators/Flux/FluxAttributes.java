// ***************************************************************************
//
// Copyright (c) 2000 - 2013, Lawrence Livermore National Security, LLC
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
// Class: FluxAttributes
//
// Purpose:
//    Attributes for Flux operator
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

public class FluxAttributes extends AttributeSubject implements Plugin
{
    private static int FluxAttributes_numAdditionalAtts = 3;

    public FluxAttributes()
    {
        super(FluxAttributes_numAdditionalAtts);

        flowField = new String("");
        weight = false;
        weightField = new String("");
    }

    public FluxAttributes(int nMoreFields)
    {
        super(FluxAttributes_numAdditionalAtts + nMoreFields);

        flowField = new String("");
        weight = false;
        weightField = new String("");
    }

    public FluxAttributes(FluxAttributes obj)
    {
        super(FluxAttributes_numAdditionalAtts);

        flowField = new String(obj.flowField);
        weight = obj.weight;
        weightField = new String(obj.weightField);

        SelectAll();
    }

    public int Offset()
    {
        return super.Offset() + super.GetNumAdditionalAttributes();
    }

    public int GetNumAdditionalAttributes()
    {
        return FluxAttributes_numAdditionalAtts;
    }

    public boolean equals(FluxAttributes obj)
    {
        // Create the return value
        return ((flowField.equals(obj.flowField)) &&
                (weight == obj.weight) &&
                (weightField.equals(obj.weightField)));
    }

    public String GetName() { return "Flux"; }
    public String GetVersion() { return "1.0"; }

    // Property setting methods
    public void SetFlowField(String flowField_)
    {
        flowField = flowField_;
        Select(0);
    }

    public void SetWeight(boolean weight_)
    {
        weight = weight_;
        Select(1);
    }

    public void SetWeightField(String weightField_)
    {
        weightField = weightField_;
        Select(2);
    }

    // Property getting methods
    public String  GetFlowField() { return flowField; }
    public boolean GetWeight() { return weight; }
    public String  GetWeightField() { return weightField; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteString(flowField);
        if(WriteSelect(1, buf))
            buf.WriteBool(weight);
        if(WriteSelect(2, buf))
            buf.WriteString(weightField);
    }

    public void ReadAtts(int index, CommunicationBuffer buf)
    {
        switch(index)
        {
        case 0:
            SetFlowField(buf.ReadString());
            break;
        case 1:
            SetWeight(buf.ReadBool());
            break;
        case 2:
            SetWeightField(buf.ReadString());
            break;
        }
    }

    public String toString(String indent)
    {
        String str = new String();
        str = str + stringToString("flowField", flowField, indent) + "\n";
        str = str + boolToString("weight", weight, indent) + "\n";
        str = str + stringToString("weightField", weightField, indent) + "\n";
        return str;
    }


    // Attributes
    private String  flowField;
    private boolean weight;
    private String  weightField;
}

