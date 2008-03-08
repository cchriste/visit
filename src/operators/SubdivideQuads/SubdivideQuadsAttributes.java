// ***************************************************************************
//
// Copyright (c) 2000 - 2008, Lawrence Livermore National Security, LLC
// Produced at the Lawrence Livermore National Laboratory
// LLNL-CODE-400142
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
// Class: SubdivideQuadsAttributes
//
// Purpose:
//    Attributes for SubdivideQuads operator
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

public class SubdivideQuadsAttributes extends AttributeSubject implements Plugin
{
    public SubdivideQuadsAttributes()
    {
        super(5);

        threshold = 0.500002;
        maxSubdivs = 4;
        fanOutPoints = true;
        doTriangles = false;
        variable = new String("default");
    }

    public SubdivideQuadsAttributes(SubdivideQuadsAttributes obj)
    {
        super(5);

        threshold = obj.threshold;
        maxSubdivs = obj.maxSubdivs;
        fanOutPoints = obj.fanOutPoints;
        doTriangles = obj.doTriangles;
        variable = new String(obj.variable);

        SelectAll();
    }

    public boolean equals(SubdivideQuadsAttributes obj)
    {
        // Create the return value
        return ((threshold == obj.threshold) &&
                (maxSubdivs == obj.maxSubdivs) &&
                (fanOutPoints == obj.fanOutPoints) &&
                (doTriangles == obj.doTriangles) &&
                (variable.equals(obj.variable)));
    }

    public String GetName() { return "SubdivideQuads"; }
    public String GetVersion() { return "1.0"; }

    // Property setting methods
    public void SetThreshold(double threshold_)
    {
        threshold = threshold_;
        Select(0);
    }

    public void SetMaxSubdivs(int maxSubdivs_)
    {
        maxSubdivs = maxSubdivs_;
        Select(1);
    }

    public void SetFanOutPoints(boolean fanOutPoints_)
    {
        fanOutPoints = fanOutPoints_;
        Select(2);
    }

    public void SetDoTriangles(boolean doTriangles_)
    {
        doTriangles = doTriangles_;
        Select(3);
    }

    public void SetVariable(String variable_)
    {
        variable = variable_;
        Select(4);
    }

    // Property getting methods
    public double  GetThreshold() { return threshold; }
    public int     GetMaxSubdivs() { return maxSubdivs; }
    public boolean GetFanOutPoints() { return fanOutPoints; }
    public boolean GetDoTriangles() { return doTriangles; }
    public String  GetVariable() { return variable; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteDouble(threshold);
        if(WriteSelect(1, buf))
            buf.WriteInt(maxSubdivs);
        if(WriteSelect(2, buf))
            buf.WriteBool(fanOutPoints);
        if(WriteSelect(3, buf))
            buf.WriteBool(doTriangles);
        if(WriteSelect(4, buf))
            buf.WriteString(variable);
    }

    public void ReadAtts(int n, CommunicationBuffer buf)
    {
        for(int i = 0; i < n; ++i)
        {
            int index = (int)buf.ReadByte();
            switch(index)
            {
            case 0:
                SetThreshold(buf.ReadDouble());
                break;
            case 1:
                SetMaxSubdivs(buf.ReadInt());
                break;
            case 2:
                SetFanOutPoints(buf.ReadBool());
                break;
            case 3:
                SetDoTriangles(buf.ReadBool());
                break;
            case 4:
                SetVariable(buf.ReadString());
                break;
            }
        }
    }

    public String toString(String indent)
    {
        String str = new String();
        str = str + doubleToString("threshold", threshold, indent) + "\n";
        str = str + intToString("maxSubdivs", maxSubdivs, indent) + "\n";
        str = str + boolToString("fanOutPoints", fanOutPoints, indent) + "\n";
        str = str + boolToString("doTriangles", doTriangles, indent) + "\n";
        str = str + stringToString("variable", variable, indent) + "\n";
        return str;
    }


    // Attributes
    private double  threshold;
    private int     maxSubdivs;
    private boolean fanOutPoints;
    private boolean doTriangles;
    private String  variable;
}

