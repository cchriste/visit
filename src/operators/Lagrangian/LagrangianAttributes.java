// ***************************************************************************
//
// Copyright (c) 2000 - 2011, Lawrence Livermore National Security, LLC
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
// Class: LagrangianAttributes
//
// Purpose:
//    Attributes for Lagrangian operator
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

public class LagrangianAttributes extends AttributeSubject implements Plugin
{
    private static int LagrangianAttributes_numAdditionalAtts = 5;

    // Enum values
    public final static int SAMPLETYPE_STEP = 0;
    public final static int SAMPLETYPE_TIME = 1;
    public final static int SAMPLETYPE_ARCLENGTH = 2;
    public final static int SAMPLETYPE_SPEED = 3;
    public final static int SAMPLETYPE_VORTICITY = 4;
    public final static int SAMPLETYPE_VARIABLE = 5;


    public LagrangianAttributes()
    {
        super(LagrangianAttributes_numAdditionalAtts);

        seedPoint = new double[3];
        seedPoint[0] = 0;
        seedPoint[1] = 0;
        seedPoint[2] = 0;
        numSteps = 1000;
        XAxisSample = SAMPLETYPE_STEP;
        YAxisSample = SAMPLETYPE_STEP;
        variable = new String("");
    }

    public LagrangianAttributes(int nMoreFields)
    {
        super(LagrangianAttributes_numAdditionalAtts + nMoreFields);

        seedPoint = new double[3];
        seedPoint[0] = 0;
        seedPoint[1] = 0;
        seedPoint[2] = 0;
        numSteps = 1000;
        XAxisSample = SAMPLETYPE_STEP;
        YAxisSample = SAMPLETYPE_STEP;
        variable = new String("");
    }

    public LagrangianAttributes(LagrangianAttributes obj)
    {
        super(LagrangianAttributes_numAdditionalAtts);

        int i;

        seedPoint = new double[3];
        seedPoint[0] = obj.seedPoint[0];
        seedPoint[1] = obj.seedPoint[1];
        seedPoint[2] = obj.seedPoint[2];

        numSteps = obj.numSteps;
        XAxisSample = obj.XAxisSample;
        YAxisSample = obj.YAxisSample;
        variable = new String(obj.variable);

        SelectAll();
    }

    public int Offset()
    {
        return super.Offset() + super.GetNumAdditionalAttributes();
    }

    public int GetNumAdditionalAttributes()
    {
        return LagrangianAttributes_numAdditionalAtts;
    }

    public boolean equals(LagrangianAttributes obj)
    {
        int i;

        // Compare the seedPoint arrays.
        boolean seedPoint_equal = true;
        for(i = 0; i < 3 && seedPoint_equal; ++i)
            seedPoint_equal = (seedPoint[i] == obj.seedPoint[i]);

        // Create the return value
        return (seedPoint_equal &&
                (numSteps == obj.numSteps) &&
                (XAxisSample == obj.XAxisSample) &&
                (YAxisSample == obj.YAxisSample) &&
                (variable.equals(obj.variable)));
    }

    public String GetName() { return "Lagrangian"; }
    public String GetVersion() { return "1.0"; }

    // Property setting methods
    public void SetSeedPoint(double[] seedPoint_)
    {
        seedPoint[0] = seedPoint_[0];
        seedPoint[1] = seedPoint_[1];
        seedPoint[2] = seedPoint_[2];
        Select(0);
    }

    public void SetSeedPoint(double e0, double e1, double e2)
    {
        seedPoint[0] = e0;
        seedPoint[1] = e1;
        seedPoint[2] = e2;
        Select(0);
    }

    public void SetNumSteps(int numSteps_)
    {
        numSteps = numSteps_;
        Select(1);
    }

    public void SetXAxisSample(int XAxisSample_)
    {
        XAxisSample = XAxisSample_;
        Select(2);
    }

    public void SetYAxisSample(int YAxisSample_)
    {
        YAxisSample = YAxisSample_;
        Select(3);
    }

    public void SetVariable(String variable_)
    {
        variable = variable_;
        Select(4);
    }

    // Property getting methods
    public double[] GetSeedPoint() { return seedPoint; }
    public int      GetNumSteps() { return numSteps; }
    public int      GetXAxisSample() { return XAxisSample; }
    public int      GetYAxisSample() { return YAxisSample; }
    public String   GetVariable() { return variable; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteDoubleArray(seedPoint);
        if(WriteSelect(1, buf))
            buf.WriteInt(numSteps);
        if(WriteSelect(2, buf))
            buf.WriteInt(XAxisSample);
        if(WriteSelect(3, buf))
            buf.WriteInt(YAxisSample);
        if(WriteSelect(4, buf))
            buf.WriteString(variable);
    }

    public void ReadAtts(int index, CommunicationBuffer buf)
    {
        switch(index)
        {
        case 0:
            SetSeedPoint(buf.ReadDoubleArray());
            break;
        case 1:
            SetNumSteps(buf.ReadInt());
            break;
        case 2:
            SetXAxisSample(buf.ReadInt());
            break;
        case 3:
            SetYAxisSample(buf.ReadInt());
            break;
        case 4:
            SetVariable(buf.ReadString());
            break;
        }
    }

    public String toString(String indent)
    {
        String str = new String();
        str = str + doubleArrayToString("seedPoint", seedPoint, indent) + "\n";
        str = str + intToString("numSteps", numSteps, indent) + "\n";
        str = str + indent + "XAxisSample = ";
        if(XAxisSample == SAMPLETYPE_STEP)
            str = str + "SAMPLETYPE_STEP";
        if(XAxisSample == SAMPLETYPE_TIME)
            str = str + "SAMPLETYPE_TIME";
        if(XAxisSample == SAMPLETYPE_ARCLENGTH)
            str = str + "SAMPLETYPE_ARCLENGTH";
        if(XAxisSample == SAMPLETYPE_SPEED)
            str = str + "SAMPLETYPE_SPEED";
        if(XAxisSample == SAMPLETYPE_VORTICITY)
            str = str + "SAMPLETYPE_VORTICITY";
        if(XAxisSample == SAMPLETYPE_VARIABLE)
            str = str + "SAMPLETYPE_VARIABLE";
        str = str + "\n";
        str = str + indent + "YAxisSample = ";
        if(YAxisSample == SAMPLETYPE_STEP)
            str = str + "SAMPLETYPE_STEP";
        if(YAxisSample == SAMPLETYPE_TIME)
            str = str + "SAMPLETYPE_TIME";
        if(YAxisSample == SAMPLETYPE_ARCLENGTH)
            str = str + "SAMPLETYPE_ARCLENGTH";
        if(YAxisSample == SAMPLETYPE_SPEED)
            str = str + "SAMPLETYPE_SPEED";
        if(YAxisSample == SAMPLETYPE_VORTICITY)
            str = str + "SAMPLETYPE_VORTICITY";
        if(YAxisSample == SAMPLETYPE_VARIABLE)
            str = str + "SAMPLETYPE_VARIABLE";
        str = str + "\n";
        str = str + stringToString("variable", variable, indent) + "\n";
        return str;
    }


    // Attributes
    private double[] seedPoint;
    private int      numSteps;
    private int      XAxisSample;
    private int      YAxisSample;
    private String   variable;
}

