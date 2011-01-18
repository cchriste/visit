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

package llnl.visit;

import java.util.Vector;
import java.lang.Double;
import java.lang.Integer;

// ****************************************************************************
// Class: ConstructDataBinningAttributes
//
// Purpose:
//    Attributes for constructing a data binning
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

public class ConstructDataBinningAttributes extends AttributeSubject
{
    private static int ConstructDataBinningAttributes_numAdditionalAtts = 13;

    // Enum values
    public final static int BINNINGSCHEME_UNIFORM = 0;
    public final static int BINNINGSCHEME_UNKNOWN = 1;

    public final static int REDUCTIONOPERATOR_AVERAGE = 0;
    public final static int REDUCTIONOPERATOR_MINIMUM = 1;
    public final static int REDUCTIONOPERATOR_MAXIMUM = 2;
    public final static int REDUCTIONOPERATOR_STANDARDDEVIATION = 3;
    public final static int REDUCTIONOPERATOR_VARIANCE = 4;
    public final static int REDUCTIONOPERATOR_SUM = 5;
    public final static int REDUCTIONOPERATOR_COUNT = 6;
    public final static int REDUCTIONOPERATOR_RMS = 7;
    public final static int REDUCTIONOPERATOR_PDF = 8;

    public final static int OUTOFBOUNDSBEHAVIOR_CLAMP = 0;
    public final static int OUTOFBOUNDSBEHAVIOR_DISCARD = 1;


    public ConstructDataBinningAttributes()
    {
        super(ConstructDataBinningAttributes_numAdditionalAtts);

        name = new String("");
        varnames = new Vector();
        binBoundaries = new Vector();
        reductionOperator = REDUCTIONOPERATOR_AVERAGE;
        varForReductionOperator = new String("");
        undefinedValue = 0;
        binningScheme = BINNINGSCHEME_UNIFORM;
        numBins = new Vector();
        overTime = false;
        timeStart = 0;
        timeEnd = 1;
        timeStride = 1;
        outOfBoundsBehavior = OUTOFBOUNDSBEHAVIOR_CLAMP;
    }

    public ConstructDataBinningAttributes(int nMoreFields)
    {
        super(ConstructDataBinningAttributes_numAdditionalAtts + nMoreFields);

        name = new String("");
        varnames = new Vector();
        binBoundaries = new Vector();
        reductionOperator = REDUCTIONOPERATOR_AVERAGE;
        varForReductionOperator = new String("");
        undefinedValue = 0;
        binningScheme = BINNINGSCHEME_UNIFORM;
        numBins = new Vector();
        overTime = false;
        timeStart = 0;
        timeEnd = 1;
        timeStride = 1;
        outOfBoundsBehavior = OUTOFBOUNDSBEHAVIOR_CLAMP;
    }

    public ConstructDataBinningAttributes(ConstructDataBinningAttributes obj)
    {
        super(ConstructDataBinningAttributes_numAdditionalAtts);

        int i;

        name = new String(obj.name);
        varnames = new Vector(obj.varnames.size());
        for(i = 0; i < obj.varnames.size(); ++i)
            varnames.addElement(new String((String)obj.varnames.elementAt(i)));

        binBoundaries = new Vector(obj.binBoundaries.size());
        for(i = 0; i < obj.binBoundaries.size(); ++i)
        {
            Double dv = (Double)obj.binBoundaries.elementAt(i);
            binBoundaries.addElement(new Double(dv.doubleValue()));
        }

        reductionOperator = obj.reductionOperator;
        varForReductionOperator = new String(obj.varForReductionOperator);
        undefinedValue = obj.undefinedValue;
        binningScheme = obj.binningScheme;
        numBins = new Vector();
        for(i = 0; i < obj.numBins.size(); ++i)
        {
            Integer iv = (Integer)obj.numBins.elementAt(i);
            numBins.addElement(new Integer(iv.intValue()));
        }
        overTime = obj.overTime;
        timeStart = obj.timeStart;
        timeEnd = obj.timeEnd;
        timeStride = obj.timeStride;
        outOfBoundsBehavior = obj.outOfBoundsBehavior;

        SelectAll();
    }

    public int Offset()
    {
        return super.Offset() + super.GetNumAdditionalAttributes();
    }

    public int GetNumAdditionalAttributes()
    {
        return ConstructDataBinningAttributes_numAdditionalAtts;
    }

    public boolean equals(ConstructDataBinningAttributes obj)
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
        // Compare the elements in the binBoundaries vector.
        boolean binBoundaries_equal = (obj.binBoundaries.size() == binBoundaries.size());
        for(i = 0; (i < binBoundaries.size()) && binBoundaries_equal; ++i)
        {
            // Make references to Double from Object.
            Double binBoundaries1 = (Double)binBoundaries.elementAt(i);
            Double binBoundaries2 = (Double)obj.binBoundaries.elementAt(i);
            binBoundaries_equal = binBoundaries1.equals(binBoundaries2);
        }
        // Compare the elements in the numBins vector.
        boolean numBins_equal = (obj.numBins.size() == numBins.size());
        for(i = 0; (i < numBins.size()) && numBins_equal; ++i)
        {
            // Make references to Integer from Object.
            Integer numBins1 = (Integer)numBins.elementAt(i);
            Integer numBins2 = (Integer)obj.numBins.elementAt(i);
            numBins_equal = numBins1.equals(numBins2);
        }
        // Create the return value
        return ((name.equals(obj.name)) &&
                varnames_equal &&
                binBoundaries_equal &&
                (reductionOperator == obj.reductionOperator) &&
                (varForReductionOperator.equals(obj.varForReductionOperator)) &&
                (undefinedValue == obj.undefinedValue) &&
                (binningScheme == obj.binningScheme) &&
                numBins_equal &&
                (overTime == obj.overTime) &&
                (timeStart == obj.timeStart) &&
                (timeEnd == obj.timeEnd) &&
                (timeStride == obj.timeStride) &&
                (outOfBoundsBehavior == obj.outOfBoundsBehavior));
    }

    // Property setting methods
    public void SetName(String name_)
    {
        name = name_;
        Select(0);
    }

    public void SetVarnames(Vector varnames_)
    {
        varnames = varnames_;
        Select(1);
    }

    public void SetBinBoundaries(Vector binBoundaries_)
    {
        binBoundaries = binBoundaries_;
        Select(2);
    }

    public void SetReductionOperator(int reductionOperator_)
    {
        reductionOperator = reductionOperator_;
        Select(3);
    }

    public void SetVarForReductionOperator(String varForReductionOperator_)
    {
        varForReductionOperator = varForReductionOperator_;
        Select(4);
    }

    public void SetUndefinedValue(double undefinedValue_)
    {
        undefinedValue = undefinedValue_;
        Select(5);
    }

    public void SetBinningScheme(int binningScheme_)
    {
        binningScheme = binningScheme_;
        Select(6);
    }

    public void SetNumBins(Vector numBins_)
    {
        numBins = numBins_;
        Select(7);
    }

    public void SetOverTime(boolean overTime_)
    {
        overTime = overTime_;
        Select(8);
    }

    public void SetTimeStart(int timeStart_)
    {
        timeStart = timeStart_;
        Select(9);
    }

    public void SetTimeEnd(int timeEnd_)
    {
        timeEnd = timeEnd_;
        Select(10);
    }

    public void SetTimeStride(int timeStride_)
    {
        timeStride = timeStride_;
        Select(11);
    }

    public void SetOutOfBoundsBehavior(int outOfBoundsBehavior_)
    {
        outOfBoundsBehavior = outOfBoundsBehavior_;
        Select(12);
    }

    // Property getting methods
    public String  GetName() { return name; }
    public Vector  GetVarnames() { return varnames; }
    public Vector  GetBinBoundaries() { return binBoundaries; }
    public int     GetReductionOperator() { return reductionOperator; }
    public String  GetVarForReductionOperator() { return varForReductionOperator; }
    public double  GetUndefinedValue() { return undefinedValue; }
    public int     GetBinningScheme() { return binningScheme; }
    public Vector  GetNumBins() { return numBins; }
    public boolean GetOverTime() { return overTime; }
    public int     GetTimeStart() { return timeStart; }
    public int     GetTimeEnd() { return timeEnd; }
    public int     GetTimeStride() { return timeStride; }
    public int     GetOutOfBoundsBehavior() { return outOfBoundsBehavior; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteString(name);
        if(WriteSelect(1, buf))
            buf.WriteStringVector(varnames);
        if(WriteSelect(2, buf))
            buf.WriteDoubleVector(binBoundaries);
        if(WriteSelect(3, buf))
            buf.WriteInt(reductionOperator);
        if(WriteSelect(4, buf))
            buf.WriteString(varForReductionOperator);
        if(WriteSelect(5, buf))
            buf.WriteDouble(undefinedValue);
        if(WriteSelect(6, buf))
            buf.WriteInt(binningScheme);
        if(WriteSelect(7, buf))
            buf.WriteIntVector(numBins);
        if(WriteSelect(8, buf))
            buf.WriteBool(overTime);
        if(WriteSelect(9, buf))
            buf.WriteInt(timeStart);
        if(WriteSelect(10, buf))
            buf.WriteInt(timeEnd);
        if(WriteSelect(11, buf))
            buf.WriteInt(timeStride);
        if(WriteSelect(12, buf))
            buf.WriteInt(outOfBoundsBehavior);
    }

    public void ReadAtts(int index, CommunicationBuffer buf)
    {
        switch(index)
        {
        case 0:
            SetName(buf.ReadString());
            break;
        case 1:
            SetVarnames(buf.ReadStringVector());
            break;
        case 2:
            SetBinBoundaries(buf.ReadDoubleVector());
            break;
        case 3:
            SetReductionOperator(buf.ReadInt());
            break;
        case 4:
            SetVarForReductionOperator(buf.ReadString());
            break;
        case 5:
            SetUndefinedValue(buf.ReadDouble());
            break;
        case 6:
            SetBinningScheme(buf.ReadInt());
            break;
        case 7:
            SetNumBins(buf.ReadIntVector());
            break;
        case 8:
            SetOverTime(buf.ReadBool());
            break;
        case 9:
            SetTimeStart(buf.ReadInt());
            break;
        case 10:
            SetTimeEnd(buf.ReadInt());
            break;
        case 11:
            SetTimeStride(buf.ReadInt());
            break;
        case 12:
            SetOutOfBoundsBehavior(buf.ReadInt());
            break;
        }
    }

    public String toString(String indent)
    {
        String str = new String();
        str = str + stringToString("name", name, indent) + "\n";
        str = str + stringVectorToString("varnames", varnames, indent) + "\n";
        str = str + doubleVectorToString("binBoundaries", binBoundaries, indent) + "\n";
        str = str + indent + "reductionOperator = ";
        if(reductionOperator == REDUCTIONOPERATOR_AVERAGE)
            str = str + "REDUCTIONOPERATOR_AVERAGE";
        if(reductionOperator == REDUCTIONOPERATOR_MINIMUM)
            str = str + "REDUCTIONOPERATOR_MINIMUM";
        if(reductionOperator == REDUCTIONOPERATOR_MAXIMUM)
            str = str + "REDUCTIONOPERATOR_MAXIMUM";
        if(reductionOperator == REDUCTIONOPERATOR_STANDARDDEVIATION)
            str = str + "REDUCTIONOPERATOR_STANDARDDEVIATION";
        if(reductionOperator == REDUCTIONOPERATOR_VARIANCE)
            str = str + "REDUCTIONOPERATOR_VARIANCE";
        if(reductionOperator == REDUCTIONOPERATOR_SUM)
            str = str + "REDUCTIONOPERATOR_SUM";
        if(reductionOperator == REDUCTIONOPERATOR_COUNT)
            str = str + "REDUCTIONOPERATOR_COUNT";
        if(reductionOperator == REDUCTIONOPERATOR_RMS)
            str = str + "REDUCTIONOPERATOR_RMS";
        if(reductionOperator == REDUCTIONOPERATOR_PDF)
            str = str + "REDUCTIONOPERATOR_PDF";
        str = str + "\n";
        str = str + stringToString("varForReductionOperator", varForReductionOperator, indent) + "\n";
        str = str + doubleToString("undefinedValue", undefinedValue, indent) + "\n";
        str = str + indent + "binningScheme = ";
        if(binningScheme == BINNINGSCHEME_UNIFORM)
            str = str + "BINNINGSCHEME_UNIFORM";
        if(binningScheme == BINNINGSCHEME_UNKNOWN)
            str = str + "BINNINGSCHEME_UNKNOWN";
        str = str + "\n";
        str = str + intVectorToString("numBins", numBins, indent) + "\n";
        str = str + boolToString("overTime", overTime, indent) + "\n";
        str = str + intToString("timeStart", timeStart, indent) + "\n";
        str = str + intToString("timeEnd", timeEnd, indent) + "\n";
        str = str + intToString("timeStride", timeStride, indent) + "\n";
        str = str + indent + "outOfBoundsBehavior = ";
        if(outOfBoundsBehavior == OUTOFBOUNDSBEHAVIOR_CLAMP)
            str = str + "OUTOFBOUNDSBEHAVIOR_CLAMP";
        if(outOfBoundsBehavior == OUTOFBOUNDSBEHAVIOR_DISCARD)
            str = str + "OUTOFBOUNDSBEHAVIOR_DISCARD";
        str = str + "\n";
        return str;
    }


    // Attributes
    private String  name;
    private Vector  varnames; // vector of String objects
    private Vector  binBoundaries; // vector of Double objects
    private int     reductionOperator;
    private String  varForReductionOperator;
    private double  undefinedValue;
    private int     binningScheme;
    private Vector  numBins; // vector of Integer objects
    private boolean overTime;
    private int     timeStart;
    private int     timeEnd;
    private int     timeStride;
    private int     outOfBoundsBehavior;
}

