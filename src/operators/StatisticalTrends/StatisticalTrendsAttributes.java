// ***************************************************************************
//
// Copyright (c) 2000 - 2016, Lawrence Livermore National Security, LLC
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
// Class: StatisticalTrendsAttributes
//
// Purpose:
//    This class contains attributes for the StatisticalTrends operator.
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

public class StatisticalTrendsAttributes extends AttributeSubject implements Plugin
{
    private static int StatisticalTrendsAttributes_numAdditionalAtts = 8;

    // Enum values
    public final static int TRENDTYPEENUM_ABSOLUTE = 0;
    public final static int TRENDTYPEENUM_RELATIVE = 1;

    public final static int STATISTICTYPEENUM_SUM = 0;
    public final static int STATISTICTYPEENUM_MEAN = 1;
    public final static int STATISTICTYPEENUM_VARIANCE = 2;
    public final static int STATISTICTYPEENUM_STANDARDDEVIATION = 3;
    public final static int STATISTICTYPEENUM_SLOPE = 4;
    public final static int STATISTICTYPEENUM_RESIDUALS = 5;

    public final static int TRENDAXISENUM_STEP = 0;
    public final static int TRENDAXISENUM_TIME = 1;
    public final static int TRENDAXISENUM_CYCLE = 2;

    public final static int VARIABLESOURCEENUM_DEFAULT = 0;
    public final static int VARIABLESOURCEENUM_OPERATOREXPRESSION = 1;


    public StatisticalTrendsAttributes()
    {
        super(StatisticalTrendsAttributes_numAdditionalAtts);

        startIndex = 0;
        stopIndex = 1;
        stride = 1;
        startTrendType = TRENDTYPEENUM_ABSOLUTE;
        stopTrendType = TRENDTYPEENUM_ABSOLUTE;
        statisticType = STATISTICTYPEENUM_MEAN;
        trendAxis = TRENDAXISENUM_STEP;
        variableSource = VARIABLESOURCEENUM_DEFAULT;
    }

    public StatisticalTrendsAttributes(int nMoreFields)
    {
        super(StatisticalTrendsAttributes_numAdditionalAtts + nMoreFields);

        startIndex = 0;
        stopIndex = 1;
        stride = 1;
        startTrendType = TRENDTYPEENUM_ABSOLUTE;
        stopTrendType = TRENDTYPEENUM_ABSOLUTE;
        statisticType = STATISTICTYPEENUM_MEAN;
        trendAxis = TRENDAXISENUM_STEP;
        variableSource = VARIABLESOURCEENUM_DEFAULT;
    }

    public StatisticalTrendsAttributes(StatisticalTrendsAttributes obj)
    {
        super(StatisticalTrendsAttributes_numAdditionalAtts);

        startIndex = obj.startIndex;
        stopIndex = obj.stopIndex;
        stride = obj.stride;
        startTrendType = obj.startTrendType;
        stopTrendType = obj.stopTrendType;
        statisticType = obj.statisticType;
        trendAxis = obj.trendAxis;
        variableSource = obj.variableSource;

        SelectAll();
    }

    public int Offset()
    {
        return super.Offset() + super.GetNumAdditionalAttributes();
    }

    public int GetNumAdditionalAttributes()
    {
        return StatisticalTrendsAttributes_numAdditionalAtts;
    }

    public boolean equals(StatisticalTrendsAttributes obj)
    {
        // Create the return value
        return ((startIndex == obj.startIndex) &&
                (stopIndex == obj.stopIndex) &&
                (stride == obj.stride) &&
                (startTrendType == obj.startTrendType) &&
                (stopTrendType == obj.stopTrendType) &&
                (statisticType == obj.statisticType) &&
                (trendAxis == obj.trendAxis) &&
                (variableSource == obj.variableSource));
    }

    public String GetName() { return "StatisticalTrends"; }
    public String GetVersion() { return "1.0"; }

    // Property setting methods
    public void SetStartIndex(int startIndex_)
    {
        startIndex = startIndex_;
        Select(0);
    }

    public void SetStopIndex(int stopIndex_)
    {
        stopIndex = stopIndex_;
        Select(1);
    }

    public void SetStride(int stride_)
    {
        stride = stride_;
        Select(2);
    }

    public void SetStartTrendType(int startTrendType_)
    {
        startTrendType = startTrendType_;
        Select(3);
    }

    public void SetStopTrendType(int stopTrendType_)
    {
        stopTrendType = stopTrendType_;
        Select(4);
    }

    public void SetStatisticType(int statisticType_)
    {
        statisticType = statisticType_;
        Select(5);
    }

    public void SetTrendAxis(int trendAxis_)
    {
        trendAxis = trendAxis_;
        Select(6);
    }

    public void SetVariableSource(int variableSource_)
    {
        variableSource = variableSource_;
        Select(7);
    }

    // Property getting methods
    public int GetStartIndex() { return startIndex; }
    public int GetStopIndex() { return stopIndex; }
    public int GetStride() { return stride; }
    public int GetStartTrendType() { return startTrendType; }
    public int GetStopTrendType() { return stopTrendType; }
    public int GetStatisticType() { return statisticType; }
    public int GetTrendAxis() { return trendAxis; }
    public int GetVariableSource() { return variableSource; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteInt(startIndex);
        if(WriteSelect(1, buf))
            buf.WriteInt(stopIndex);
        if(WriteSelect(2, buf))
            buf.WriteInt(stride);
        if(WriteSelect(3, buf))
            buf.WriteInt(startTrendType);
        if(WriteSelect(4, buf))
            buf.WriteInt(stopTrendType);
        if(WriteSelect(5, buf))
            buf.WriteInt(statisticType);
        if(WriteSelect(6, buf))
            buf.WriteInt(trendAxis);
        if(WriteSelect(7, buf))
            buf.WriteInt(variableSource);
    }

    public void ReadAtts(int index, CommunicationBuffer buf)
    {
        switch(index)
        {
        case 0:
            SetStartIndex(buf.ReadInt());
            break;
        case 1:
            SetStopIndex(buf.ReadInt());
            break;
        case 2:
            SetStride(buf.ReadInt());
            break;
        case 3:
            SetStartTrendType(buf.ReadInt());
            break;
        case 4:
            SetStopTrendType(buf.ReadInt());
            break;
        case 5:
            SetStatisticType(buf.ReadInt());
            break;
        case 6:
            SetTrendAxis(buf.ReadInt());
            break;
        case 7:
            SetVariableSource(buf.ReadInt());
            break;
        }
    }

    public String toString(String indent)
    {
        String str = new String();
        str = str + intToString("startIndex", startIndex, indent) + "\n";
        str = str + intToString("stopIndex", stopIndex, indent) + "\n";
        str = str + intToString("stride", stride, indent) + "\n";
        str = str + indent + "startTrendType = ";
        if(startTrendType == TRENDTYPEENUM_ABSOLUTE)
            str = str + "TRENDTYPEENUM_ABSOLUTE";
        if(startTrendType == TRENDTYPEENUM_RELATIVE)
            str = str + "TRENDTYPEENUM_RELATIVE";
        str = str + "\n";
        str = str + indent + "stopTrendType = ";
        if(stopTrendType == TRENDTYPEENUM_ABSOLUTE)
            str = str + "TRENDTYPEENUM_ABSOLUTE";
        if(stopTrendType == TRENDTYPEENUM_RELATIVE)
            str = str + "TRENDTYPEENUM_RELATIVE";
        str = str + "\n";
        str = str + indent + "statisticType = ";
        if(statisticType == STATISTICTYPEENUM_SUM)
            str = str + "STATISTICTYPEENUM_SUM";
        if(statisticType == STATISTICTYPEENUM_MEAN)
            str = str + "STATISTICTYPEENUM_MEAN";
        if(statisticType == STATISTICTYPEENUM_VARIANCE)
            str = str + "STATISTICTYPEENUM_VARIANCE";
        if(statisticType == STATISTICTYPEENUM_STANDARDDEVIATION)
            str = str + "STATISTICTYPEENUM_STANDARDDEVIATION";
        if(statisticType == STATISTICTYPEENUM_SLOPE)
            str = str + "STATISTICTYPEENUM_SLOPE";
        if(statisticType == STATISTICTYPEENUM_RESIDUALS)
            str = str + "STATISTICTYPEENUM_RESIDUALS";
        str = str + "\n";
        str = str + indent + "trendAxis = ";
        if(trendAxis == TRENDAXISENUM_STEP)
            str = str + "TRENDAXISENUM_STEP";
        if(trendAxis == TRENDAXISENUM_TIME)
            str = str + "TRENDAXISENUM_TIME";
        if(trendAxis == TRENDAXISENUM_CYCLE)
            str = str + "TRENDAXISENUM_CYCLE";
        str = str + "\n";
        str = str + indent + "variableSource = ";
        if(variableSource == VARIABLESOURCEENUM_DEFAULT)
            str = str + "VARIABLESOURCEENUM_DEFAULT";
        if(variableSource == VARIABLESOURCEENUM_OPERATOREXPRESSION)
            str = str + "VARIABLESOURCEENUM_OPERATOREXPRESSION";
        str = str + "\n";
        return str;
    }


    // Attributes
    private int startIndex;
    private int stopIndex;
    private int stride;
    private int startTrendType;
    private int stopTrendType;
    private int statisticType;
    private int trendAxis;
    private int variableSource;
}

