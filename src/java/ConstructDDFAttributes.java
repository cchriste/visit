// ***************************************************************************
//
// Copyright (c) 2000 - 2009, Lawrence Livermore National Security, LLC
// Produced at the Lawrence Livermore National Laboratory
// LLNL-CODE-400124
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

package llnl.visit;

import java.util.Vector;
import java.lang.Double;
import java.lang.Integer;

// ****************************************************************************
// Class: ConstructDDFAttributes
//
// Purpose:
//    Attributes for constructing DDFs
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

public class ConstructDDFAttributes extends AttributeSubject
{
    // Enum values
    public final static int BINNINGSCHEME_UNIFORM = 0;
    public final static int BINNINGSCHEME_UNKNOWN = 1;

    public final static int STATISTICALOPERATOR_AVERAGE = 0;
    public final static int STATISTICALOPERATOR_MINIMUM = 1;
    public final static int STATISTICALOPERATOR_MAXIMUM = 2;
    public final static int STATISTICALOPERATOR_STANDARDDEVIATION = 3;
    public final static int STATISTICALOPERATOR_VARIANCE = 4;
    public final static int STATISTICALOPERATOR_PERCENTILE = 5;
    public final static int STATISTICALOPERATOR_MEDIAN = 6;
    public final static int STATISTICALOPERATOR_INTERQUARTILERANGE = 7;
    public final static int STATISTICALOPERATOR_SUM = 8;
    public final static int STATISTICALOPERATOR_COUNT = 9;


    public ConstructDDFAttributes()
    {
        super(13);

        ddfName = new String("");
        varnames = new Vector();
        ranges = new Vector();
        codomainName = new String("");
        statisticalOperator = STATISTICALOPERATOR_AVERAGE;
        percentile = 90;
        undefinedValue = 0;
        binningScheme = BINNINGSCHEME_UNIFORM;
        numSamples = new Vector();
        overTime = false;
        timeStart = 0;
        timeEnd = 1;
        timeStride = 1;
    }

    public ConstructDDFAttributes(ConstructDDFAttributes obj)
    {
        super(13);

        int i;

        ddfName = new String(obj.ddfName);
        varnames = new Vector(obj.varnames.size());
        for(i = 0; i < obj.varnames.size(); ++i)
            varnames.addElement(new String((String)obj.varnames.elementAt(i)));

        ranges = new Vector(obj.ranges.size());
        for(i = 0; i < obj.ranges.size(); ++i)
        {
            Double dv = (Double)obj.ranges.elementAt(i);
            ranges.addElement(new Double(dv.doubleValue()));
        }

        codomainName = new String(obj.codomainName);
        statisticalOperator = obj.statisticalOperator;
        percentile = obj.percentile;
        undefinedValue = obj.undefinedValue;
        binningScheme = obj.binningScheme;
        numSamples = new Vector();
        for(i = 0; i < obj.numSamples.size(); ++i)
        {
            Integer iv = (Integer)obj.numSamples.elementAt(i);
            numSamples.addElement(new Integer(iv.intValue()));
        }
        overTime = obj.overTime;
        timeStart = obj.timeStart;
        timeEnd = obj.timeEnd;
        timeStride = obj.timeStride;

        SelectAll();
    }

    public boolean equals(ConstructDDFAttributes obj)
    {
        int i;

        // Compare the elements in the varnames vector.
        boolean varnames_equal = (obj.varnames.size() == varnames.size());
        for(i = 0; (i < varnames.size()) && varnames_equal; ++i)
        {
            // Make references to String from Object.
            String varnames1 = (String)varnames.elementAt(i);
            String varnames2 = (String)obj.varnames.elementAt(i);
            varnames_equal = varnames1.equals(varnames2);
        }
        // Compare the elements in the ranges vector.
        boolean ranges_equal = (obj.ranges.size() == ranges.size());
        for(i = 0; (i < ranges.size()) && ranges_equal; ++i)
        {
            // Make references to Double from Object.
            Double ranges1 = (Double)ranges.elementAt(i);
            Double ranges2 = (Double)obj.ranges.elementAt(i);
            ranges_equal = ranges1.equals(ranges2);
        }
        // Compare the elements in the numSamples vector.
        boolean numSamples_equal = (obj.numSamples.size() == numSamples.size());
        for(i = 0; (i < numSamples.size()) && numSamples_equal; ++i)
        {
            // Make references to Integer from Object.
            Integer numSamples1 = (Integer)numSamples.elementAt(i);
            Integer numSamples2 = (Integer)obj.numSamples.elementAt(i);
            numSamples_equal = numSamples1.equals(numSamples2);
        }
        // Create the return value
        return ((ddfName.equals(obj.ddfName)) &&
                varnames_equal &&
                ranges_equal &&
                (codomainName.equals(obj.codomainName)) &&
                (statisticalOperator == obj.statisticalOperator) &&
                (percentile == obj.percentile) &&
                (undefinedValue == obj.undefinedValue) &&
                (binningScheme == obj.binningScheme) &&
                numSamples_equal &&
                (overTime == obj.overTime) &&
                (timeStart == obj.timeStart) &&
                (timeEnd == obj.timeEnd) &&
                (timeStride == obj.timeStride));
    }

    // Property setting methods
    public void SetDdfName(String ddfName_)
    {
        ddfName = ddfName_;
        Select(0);
    }

    public void SetVarnames(Vector varnames_)
    {
        varnames = varnames_;
        Select(1);
    }

    public void SetRanges(Vector ranges_)
    {
        ranges = ranges_;
        Select(2);
    }

    public void SetCodomainName(String codomainName_)
    {
        codomainName = codomainName_;
        Select(3);
    }

    public void SetStatisticalOperator(int statisticalOperator_)
    {
        statisticalOperator = statisticalOperator_;
        Select(4);
    }

    public void SetPercentile(double percentile_)
    {
        percentile = percentile_;
        Select(5);
    }

    public void SetUndefinedValue(double undefinedValue_)
    {
        undefinedValue = undefinedValue_;
        Select(6);
    }

    public void SetBinningScheme(int binningScheme_)
    {
        binningScheme = binningScheme_;
        Select(7);
    }

    public void SetNumSamples(Vector numSamples_)
    {
        numSamples = numSamples_;
        Select(8);
    }

    public void SetOverTime(boolean overTime_)
    {
        overTime = overTime_;
        Select(9);
    }

    public void SetTimeStart(int timeStart_)
    {
        timeStart = timeStart_;
        Select(10);
    }

    public void SetTimeEnd(int timeEnd_)
    {
        timeEnd = timeEnd_;
        Select(11);
    }

    public void SetTimeStride(int timeStride_)
    {
        timeStride = timeStride_;
        Select(12);
    }

    // Property getting methods
    public String  GetDdfName() { return ddfName; }
    public Vector  GetVarnames() { return varnames; }
    public Vector  GetRanges() { return ranges; }
    public String  GetCodomainName() { return codomainName; }
    public int     GetStatisticalOperator() { return statisticalOperator; }
    public double  GetPercentile() { return percentile; }
    public double  GetUndefinedValue() { return undefinedValue; }
    public int     GetBinningScheme() { return binningScheme; }
    public Vector  GetNumSamples() { return numSamples; }
    public boolean GetOverTime() { return overTime; }
    public int     GetTimeStart() { return timeStart; }
    public int     GetTimeEnd() { return timeEnd; }
    public int     GetTimeStride() { return timeStride; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteString(ddfName);
        if(WriteSelect(1, buf))
            buf.WriteStringVector(varnames);
        if(WriteSelect(2, buf))
            buf.WriteDoubleVector(ranges);
        if(WriteSelect(3, buf))
            buf.WriteString(codomainName);
        if(WriteSelect(4, buf))
            buf.WriteInt(statisticalOperator);
        if(WriteSelect(5, buf))
            buf.WriteDouble(percentile);
        if(WriteSelect(6, buf))
            buf.WriteDouble(undefinedValue);
        if(WriteSelect(7, buf))
            buf.WriteInt(binningScheme);
        if(WriteSelect(8, buf))
            buf.WriteIntVector(numSamples);
        if(WriteSelect(9, buf))
            buf.WriteBool(overTime);
        if(WriteSelect(10, buf))
            buf.WriteInt(timeStart);
        if(WriteSelect(11, buf))
            buf.WriteInt(timeEnd);
        if(WriteSelect(12, buf))
            buf.WriteInt(timeStride);
    }

    public void ReadAtts(int n, CommunicationBuffer buf)
    {
        for(int i = 0; i < n; ++i)
        {
            int index = (int)buf.ReadByte();
            switch(index)
            {
            case 0:
                SetDdfName(buf.ReadString());
                break;
            case 1:
                SetVarnames(buf.ReadStringVector());
                break;
            case 2:
                SetRanges(buf.ReadDoubleVector());
                break;
            case 3:
                SetCodomainName(buf.ReadString());
                break;
            case 4:
                SetStatisticalOperator(buf.ReadInt());
                break;
            case 5:
                SetPercentile(buf.ReadDouble());
                break;
            case 6:
                SetUndefinedValue(buf.ReadDouble());
                break;
            case 7:
                SetBinningScheme(buf.ReadInt());
                break;
            case 8:
                SetNumSamples(buf.ReadIntVector());
                break;
            case 9:
                SetOverTime(buf.ReadBool());
                break;
            case 10:
                SetTimeStart(buf.ReadInt());
                break;
            case 11:
                SetTimeEnd(buf.ReadInt());
                break;
            case 12:
                SetTimeStride(buf.ReadInt());
                break;
            }
        }
    }

    public String toString(String indent)
    {
        String str = new String();
        str = str + stringToString("ddfName", ddfName, indent) + "\n";
        str = str + stringVectorToString("varnames", varnames, indent) + "\n";
        str = str + doubleVectorToString("ranges", ranges, indent) + "\n";
        str = str + stringToString("codomainName", codomainName, indent) + "\n";
        str = str + indent + "statisticalOperator = ";
        if(statisticalOperator == STATISTICALOPERATOR_AVERAGE)
            str = str + "STATISTICALOPERATOR_AVERAGE";
        if(statisticalOperator == STATISTICALOPERATOR_MINIMUM)
            str = str + "STATISTICALOPERATOR_MINIMUM";
        if(statisticalOperator == STATISTICALOPERATOR_MAXIMUM)
            str = str + "STATISTICALOPERATOR_MAXIMUM";
        if(statisticalOperator == STATISTICALOPERATOR_STANDARDDEVIATION)
            str = str + "STATISTICALOPERATOR_STANDARDDEVIATION";
        if(statisticalOperator == STATISTICALOPERATOR_VARIANCE)
            str = str + "STATISTICALOPERATOR_VARIANCE";
        if(statisticalOperator == STATISTICALOPERATOR_PERCENTILE)
            str = str + "STATISTICALOPERATOR_PERCENTILE";
        if(statisticalOperator == STATISTICALOPERATOR_MEDIAN)
            str = str + "STATISTICALOPERATOR_MEDIAN";
        if(statisticalOperator == STATISTICALOPERATOR_INTERQUARTILERANGE)
            str = str + "STATISTICALOPERATOR_INTERQUARTILERANGE";
        if(statisticalOperator == STATISTICALOPERATOR_SUM)
            str = str + "STATISTICALOPERATOR_SUM";
        if(statisticalOperator == STATISTICALOPERATOR_COUNT)
            str = str + "STATISTICALOPERATOR_COUNT";
        str = str + "\n";
        str = str + doubleToString("percentile", percentile, indent) + "\n";
        str = str + doubleToString("undefinedValue", undefinedValue, indent) + "\n";
        str = str + indent + "binningScheme = ";
        if(binningScheme == BINNINGSCHEME_UNIFORM)
            str = str + "BINNINGSCHEME_UNIFORM";
        if(binningScheme == BINNINGSCHEME_UNKNOWN)
            str = str + "BINNINGSCHEME_UNKNOWN";
        str = str + "\n";
        str = str + intVectorToString("numSamples", numSamples, indent) + "\n";
        str = str + boolToString("overTime", overTime, indent) + "\n";
        str = str + intToString("timeStart", timeStart, indent) + "\n";
        str = str + intToString("timeEnd", timeEnd, indent) + "\n";
        str = str + intToString("timeStride", timeStride, indent) + "\n";
        return str;
    }


    // Attributes
    private String  ddfName;
    private Vector  varnames; // vector of String objects
    private Vector  ranges; // vector of Double objects
    private String  codomainName;
    private int     statisticalOperator;
    private double  percentile;
    private double  undefinedValue;
    private int     binningScheme;
    private Vector  numSamples; // vector of Integer objects
    private boolean overTime;
    private int     timeStart;
    private int     timeEnd;
    private int     timeStride;
}

