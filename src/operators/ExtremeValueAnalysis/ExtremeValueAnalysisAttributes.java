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
// Class: ExtremeValueAnalysisAttributes
//
// Purpose:
//    Attributes for ExtremeValueAnalysis operator
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

public class ExtremeValueAnalysisAttributes extends AttributeSubject implements Plugin
{
    private static int ExtremeValueAnalysisAttributes_numAdditionalAtts = 2;

    // Enum values
    public final static int COMPUTEMAXES_MONTHLY = 0;
    public final static int COMPUTEMAXES_YEARLY = 1;

    public final static int DISPLAYVALUES_A = 0;
    public final static int DISPLAYVALUES_B = 1;
    public final static int DISPLAYVALUES_C = 2;
    public final static int DISPLAYVALUES_D = 3;


    public ExtremeValueAnalysisAttributes()
    {
        super(ExtremeValueAnalysisAttributes_numAdditionalAtts);

        computeMaxes = COMPUTEMAXES_MONTHLY;
        DisplayValue = DISPLAYVALUES_A;
    }

    public ExtremeValueAnalysisAttributes(int nMoreFields)
    {
        super(ExtremeValueAnalysisAttributes_numAdditionalAtts + nMoreFields);

        computeMaxes = COMPUTEMAXES_MONTHLY;
        DisplayValue = DISPLAYVALUES_A;
    }

    public ExtremeValueAnalysisAttributes(ExtremeValueAnalysisAttributes obj)
    {
        super(ExtremeValueAnalysisAttributes_numAdditionalAtts);

        computeMaxes = obj.computeMaxes;
        DisplayValue = obj.DisplayValue;

        SelectAll();
    }

    public int Offset()
    {
        return super.Offset() + super.GetNumAdditionalAttributes();
    }

    public int GetNumAdditionalAttributes()
    {
        return ExtremeValueAnalysisAttributes_numAdditionalAtts;
    }

    public boolean equals(ExtremeValueAnalysisAttributes obj)
    {
        // Create the return value
        return ((computeMaxes == obj.computeMaxes) &&
                (DisplayValue == obj.DisplayValue));
    }

    public String GetName() { return "ExtremeValueAnalysis"; }
    public String GetVersion() { return "1.0"; }

    // Property setting methods
    public void SetComputeMaxes(int computeMaxes_)
    {
        computeMaxes = computeMaxes_;
        Select(0);
    }

    public void SetDisplayValue(int DisplayValue_)
    {
        DisplayValue = DisplayValue_;
        Select(1);
    }

    // Property getting methods
    public int GetComputeMaxes() { return computeMaxes; }
    public int GetDisplayValue() { return DisplayValue; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteInt(computeMaxes);
        if(WriteSelect(1, buf))
            buf.WriteInt(DisplayValue);
    }

    public void ReadAtts(int index, CommunicationBuffer buf)
    {
        switch(index)
        {
        case 0:
            SetComputeMaxes(buf.ReadInt());
            break;
        case 1:
            SetDisplayValue(buf.ReadInt());
            break;
        }
    }

    public String toString(String indent)
    {
        String str = new String();
        str = str + indent + "computeMaxes = ";
        if(computeMaxes == COMPUTEMAXES_MONTHLY)
            str = str + "COMPUTEMAXES_MONTHLY";
        if(computeMaxes == COMPUTEMAXES_YEARLY)
            str = str + "COMPUTEMAXES_YEARLY";
        str = str + "\n";
        str = str + indent + "DisplayValue = ";
        if(DisplayValue == DISPLAYVALUES_A)
            str = str + "DISPLAYVALUES_A";
        if(DisplayValue == DISPLAYVALUES_B)
            str = str + "DISPLAYVALUES_B";
        if(DisplayValue == DISPLAYVALUES_C)
            str = str + "DISPLAYVALUES_C";
        if(DisplayValue == DISPLAYVALUES_D)
            str = str + "DISPLAYVALUES_D";
        str = str + "\n";
        return str;
    }


    // Attributes
    private int computeMaxes;
    private int DisplayValue;
}

