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
import java.lang.Integer;
import java.util.Vector;

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
    private static int ExtremeValueAnalysisAttributes_numAdditionalAtts = 23;

    // Enum values
    public final static int AGGREGATIONTYPE_ANNUAL = 0;
    public final static int AGGREGATIONTYPE_SEASONAL = 1;
    public final static int AGGREGATIONTYPE_MONTHLY = 2;

    public final static int MONTHTYPE_JANUARY = 0;
    public final static int MONTHTYPE_FEBRUARY = 1;
    public final static int MONTHTYPE_MARCH = 2;
    public final static int MONTHTYPE_APRIL = 3;
    public final static int MONTHTYPE_MAY = 4;
    public final static int MONTHTYPE_JUNE = 5;
    public final static int MONTHTYPE_JULY = 6;
    public final static int MONTHTYPE_AUGUST = 7;
    public final static int MONTHTYPE_SEPTEMBER = 8;
    public final static int MONTHTYPE_OCTOBER = 9;
    public final static int MONTHTYPE_NOVEMBER = 10;
    public final static int MONTHTYPE_DECEMBER = 11;

    public final static int SEASONTYPE_WINTER = 0;
    public final static int SEASONTYPE_SPRING = 1;
    public final static int SEASONTYPE_SUMMER = 2;
    public final static int SEASONTYPE_FALL = 3;

    public final static int OPTIMIZATIONTYPE_NELDER_MEAD = 0;
    public final static int OPTIMIZATIONTYPE_BFGS = 1;

    public final static int EXTREMETYPE_MINIMA = 0;
    public final static int EXTREMETYPE_MAXIMA = 1;


    public ExtremeValueAnalysisAttributes()
    {
        super(ExtremeValueAnalysisAttributes_numAdditionalAtts);

        dataYearBegin = 1;
        dataAnalysisYearRangeEnabled = false;
        dataAnalysisYear1 = 1;
        dataAnalysisYear2 = 1;
        ensemble = false;
        numEnsembles = 1;
        dataScaling = 1;
        extremeMethod = EXTREMETYPE_MAXIMA;
        optimizationMethod = OPTIMIZATIONTYPE_NELDER_MEAD;
        aggregation = AGGREGATIONTYPE_ANNUAL;
        covariateModelScale = false;
        covariateModelLocation = false;
        covariateModelShape = false;
        computeReturnValues = false;
        returnValues = new Vector();
        returnValues.addElement(new Integer(1));
        computeRVDifferences = false;
        rvDifference1 = 1;
        rvDifference2 = 1;
        displayMonth = MONTHTYPE_JANUARY;
        displaySeason = SEASONTYPE_WINTER;
        computeParamValues = false;
        dumpData = true;
        dumpDebug = false;
    }

    public ExtremeValueAnalysisAttributes(int nMoreFields)
    {
        super(ExtremeValueAnalysisAttributes_numAdditionalAtts + nMoreFields);

        dataYearBegin = 1;
        dataAnalysisYearRangeEnabled = false;
        dataAnalysisYear1 = 1;
        dataAnalysisYear2 = 1;
        ensemble = false;
        numEnsembles = 1;
        dataScaling = 1;
        extremeMethod = EXTREMETYPE_MAXIMA;
        optimizationMethod = OPTIMIZATIONTYPE_NELDER_MEAD;
        aggregation = AGGREGATIONTYPE_ANNUAL;
        covariateModelScale = false;
        covariateModelLocation = false;
        covariateModelShape = false;
        computeReturnValues = false;
        returnValues = new Vector();
        returnValues.addElement(new Integer(1));
        computeRVDifferences = false;
        rvDifference1 = 1;
        rvDifference2 = 1;
        displayMonth = MONTHTYPE_JANUARY;
        displaySeason = SEASONTYPE_WINTER;
        computeParamValues = false;
        dumpData = true;
        dumpDebug = false;
    }

    public ExtremeValueAnalysisAttributes(ExtremeValueAnalysisAttributes obj)
    {
        super(ExtremeValueAnalysisAttributes_numAdditionalAtts);

        int i;

        dataYearBegin = obj.dataYearBegin;
        dataAnalysisYearRangeEnabled = obj.dataAnalysisYearRangeEnabled;
        dataAnalysisYear1 = obj.dataAnalysisYear1;
        dataAnalysisYear2 = obj.dataAnalysisYear2;
        ensemble = obj.ensemble;
        numEnsembles = obj.numEnsembles;
        dataScaling = obj.dataScaling;
        extremeMethod = obj.extremeMethod;
        optimizationMethod = obj.optimizationMethod;
        aggregation = obj.aggregation;
        covariateModelScale = obj.covariateModelScale;
        covariateModelLocation = obj.covariateModelLocation;
        covariateModelShape = obj.covariateModelShape;
        computeReturnValues = obj.computeReturnValues;
        returnValues = new Vector();
        for(i = 0; i < obj.returnValues.size(); ++i)
        {
            Integer iv = (Integer)obj.returnValues.elementAt(i);
            returnValues.addElement(new Integer(iv.intValue()));
        }
        computeRVDifferences = obj.computeRVDifferences;
        rvDifference1 = obj.rvDifference1;
        rvDifference2 = obj.rvDifference2;
        displayMonth = obj.displayMonth;
        displaySeason = obj.displaySeason;
        computeParamValues = obj.computeParamValues;
        dumpData = obj.dumpData;
        dumpDebug = obj.dumpDebug;

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
        int i;

        // Compare the elements in the returnValues vector.
        boolean returnValues_equal = (obj.returnValues.size() == returnValues.size());
        for(i = 0; (i < returnValues.size()) && returnValues_equal; ++i)
        {
            // Make references to Integer from Object.
            Integer returnValues1 = (Integer)returnValues.elementAt(i);
            Integer returnValues2 = (Integer)obj.returnValues.elementAt(i);
            returnValues_equal = returnValues1.equals(returnValues2);
        }
        // Create the return value
        return ((dataYearBegin == obj.dataYearBegin) &&
                (dataAnalysisYearRangeEnabled == obj.dataAnalysisYearRangeEnabled) &&
                (dataAnalysisYear1 == obj.dataAnalysisYear1) &&
                (dataAnalysisYear2 == obj.dataAnalysisYear2) &&
                (ensemble == obj.ensemble) &&
                (numEnsembles == obj.numEnsembles) &&
                (dataScaling == obj.dataScaling) &&
                (extremeMethod == obj.extremeMethod) &&
                (optimizationMethod == obj.optimizationMethod) &&
                (aggregation == obj.aggregation) &&
                (covariateModelScale == obj.covariateModelScale) &&
                (covariateModelLocation == obj.covariateModelLocation) &&
                (covariateModelShape == obj.covariateModelShape) &&
                (computeReturnValues == obj.computeReturnValues) &&
                returnValues_equal &&
                (computeRVDifferences == obj.computeRVDifferences) &&
                (rvDifference1 == obj.rvDifference1) &&
                (rvDifference2 == obj.rvDifference2) &&
                (displayMonth == obj.displayMonth) &&
                (displaySeason == obj.displaySeason) &&
                (computeParamValues == obj.computeParamValues) &&
                (dumpData == obj.dumpData) &&
                (dumpDebug == obj.dumpDebug));
    }

    public String GetName() { return "ExtremeValueAnalysis"; }
    public String GetVersion() { return "1.0"; }

    // Property setting methods
    public void SetDataYearBegin(int dataYearBegin_)
    {
        dataYearBegin = dataYearBegin_;
        Select(0);
    }

    public void SetDataAnalysisYearRangeEnabled(boolean dataAnalysisYearRangeEnabled_)
    {
        dataAnalysisYearRangeEnabled = dataAnalysisYearRangeEnabled_;
        Select(1);
    }

    public void SetDataAnalysisYear1(int dataAnalysisYear1_)
    {
        dataAnalysisYear1 = dataAnalysisYear1_;
        Select(2);
    }

    public void SetDataAnalysisYear2(int dataAnalysisYear2_)
    {
        dataAnalysisYear2 = dataAnalysisYear2_;
        Select(3);
    }

    public void SetEnsemble(boolean ensemble_)
    {
        ensemble = ensemble_;
        Select(4);
    }

    public void SetNumEnsembles(int numEnsembles_)
    {
        numEnsembles = numEnsembles_;
        Select(5);
    }

    public void SetDataScaling(double dataScaling_)
    {
        dataScaling = dataScaling_;
        Select(6);
    }

    public void SetExtremeMethod(int extremeMethod_)
    {
        extremeMethod = extremeMethod_;
        Select(7);
    }

    public void SetOptimizationMethod(int optimizationMethod_)
    {
        optimizationMethod = optimizationMethod_;
        Select(8);
    }

    public void SetAggregation(int aggregation_)
    {
        aggregation = aggregation_;
        Select(9);
    }

    public void SetCovariateModelScale(boolean covariateModelScale_)
    {
        covariateModelScale = covariateModelScale_;
        Select(10);
    }

    public void SetCovariateModelLocation(boolean covariateModelLocation_)
    {
        covariateModelLocation = covariateModelLocation_;
        Select(11);
    }

    public void SetCovariateModelShape(boolean covariateModelShape_)
    {
        covariateModelShape = covariateModelShape_;
        Select(12);
    }

    public void SetComputeReturnValues(boolean computeReturnValues_)
    {
        computeReturnValues = computeReturnValues_;
        Select(13);
    }

    public void SetReturnValues(Vector returnValues_)
    {
        returnValues = returnValues_;
        Select(14);
    }

    public void SetComputeRVDifferences(boolean computeRVDifferences_)
    {
        computeRVDifferences = computeRVDifferences_;
        Select(15);
    }

    public void SetRvDifference1(int rvDifference1_)
    {
        rvDifference1 = rvDifference1_;
        Select(16);
    }

    public void SetRvDifference2(int rvDifference2_)
    {
        rvDifference2 = rvDifference2_;
        Select(17);
    }

    public void SetDisplayMonth(int displayMonth_)
    {
        displayMonth = displayMonth_;
        Select(18);
    }

    public void SetDisplaySeason(int displaySeason_)
    {
        displaySeason = displaySeason_;
        Select(19);
    }

    public void SetComputeParamValues(boolean computeParamValues_)
    {
        computeParamValues = computeParamValues_;
        Select(20);
    }

    public void SetDumpData(boolean dumpData_)
    {
        dumpData = dumpData_;
        Select(21);
    }

    public void SetDumpDebug(boolean dumpDebug_)
    {
        dumpDebug = dumpDebug_;
        Select(22);
    }

    // Property getting methods
    public int     GetDataYearBegin() { return dataYearBegin; }
    public boolean GetDataAnalysisYearRangeEnabled() { return dataAnalysisYearRangeEnabled; }
    public int     GetDataAnalysisYear1() { return dataAnalysisYear1; }
    public int     GetDataAnalysisYear2() { return dataAnalysisYear2; }
    public boolean GetEnsemble() { return ensemble; }
    public int     GetNumEnsembles() { return numEnsembles; }
    public double  GetDataScaling() { return dataScaling; }
    public int     GetExtremeMethod() { return extremeMethod; }
    public int     GetOptimizationMethod() { return optimizationMethod; }
    public int     GetAggregation() { return aggregation; }
    public boolean GetCovariateModelScale() { return covariateModelScale; }
    public boolean GetCovariateModelLocation() { return covariateModelLocation; }
    public boolean GetCovariateModelShape() { return covariateModelShape; }
    public boolean GetComputeReturnValues() { return computeReturnValues; }
    public Vector  GetReturnValues() { return returnValues; }
    public boolean GetComputeRVDifferences() { return computeRVDifferences; }
    public int     GetRvDifference1() { return rvDifference1; }
    public int     GetRvDifference2() { return rvDifference2; }
    public int     GetDisplayMonth() { return displayMonth; }
    public int     GetDisplaySeason() { return displaySeason; }
    public boolean GetComputeParamValues() { return computeParamValues; }
    public boolean GetDumpData() { return dumpData; }
    public boolean GetDumpDebug() { return dumpDebug; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteInt(dataYearBegin);
        if(WriteSelect(1, buf))
            buf.WriteBool(dataAnalysisYearRangeEnabled);
        if(WriteSelect(2, buf))
            buf.WriteInt(dataAnalysisYear1);
        if(WriteSelect(3, buf))
            buf.WriteInt(dataAnalysisYear2);
        if(WriteSelect(4, buf))
            buf.WriteBool(ensemble);
        if(WriteSelect(5, buf))
            buf.WriteInt(numEnsembles);
        if(WriteSelect(6, buf))
            buf.WriteDouble(dataScaling);
        if(WriteSelect(7, buf))
            buf.WriteInt(extremeMethod);
        if(WriteSelect(8, buf))
            buf.WriteInt(optimizationMethod);
        if(WriteSelect(9, buf))
            buf.WriteInt(aggregation);
        if(WriteSelect(10, buf))
            buf.WriteBool(covariateModelScale);
        if(WriteSelect(11, buf))
            buf.WriteBool(covariateModelLocation);
        if(WriteSelect(12, buf))
            buf.WriteBool(covariateModelShape);
        if(WriteSelect(13, buf))
            buf.WriteBool(computeReturnValues);
        if(WriteSelect(14, buf))
            buf.WriteIntVector(returnValues);
        if(WriteSelect(15, buf))
            buf.WriteBool(computeRVDifferences);
        if(WriteSelect(16, buf))
            buf.WriteInt(rvDifference1);
        if(WriteSelect(17, buf))
            buf.WriteInt(rvDifference2);
        if(WriteSelect(18, buf))
            buf.WriteInt(displayMonth);
        if(WriteSelect(19, buf))
            buf.WriteInt(displaySeason);
        if(WriteSelect(20, buf))
            buf.WriteBool(computeParamValues);
        if(WriteSelect(21, buf))
            buf.WriteBool(dumpData);
        if(WriteSelect(22, buf))
            buf.WriteBool(dumpDebug);
    }

    public void ReadAtts(int index, CommunicationBuffer buf)
    {
        switch(index)
        {
        case 0:
            SetDataYearBegin(buf.ReadInt());
            break;
        case 1:
            SetDataAnalysisYearRangeEnabled(buf.ReadBool());
            break;
        case 2:
            SetDataAnalysisYear1(buf.ReadInt());
            break;
        case 3:
            SetDataAnalysisYear2(buf.ReadInt());
            break;
        case 4:
            SetEnsemble(buf.ReadBool());
            break;
        case 5:
            SetNumEnsembles(buf.ReadInt());
            break;
        case 6:
            SetDataScaling(buf.ReadDouble());
            break;
        case 7:
            SetExtremeMethod(buf.ReadInt());
            break;
        case 8:
            SetOptimizationMethod(buf.ReadInt());
            break;
        case 9:
            SetAggregation(buf.ReadInt());
            break;
        case 10:
            SetCovariateModelScale(buf.ReadBool());
            break;
        case 11:
            SetCovariateModelLocation(buf.ReadBool());
            break;
        case 12:
            SetCovariateModelShape(buf.ReadBool());
            break;
        case 13:
            SetComputeReturnValues(buf.ReadBool());
            break;
        case 14:
            SetReturnValues(buf.ReadIntVector());
            break;
        case 15:
            SetComputeRVDifferences(buf.ReadBool());
            break;
        case 16:
            SetRvDifference1(buf.ReadInt());
            break;
        case 17:
            SetRvDifference2(buf.ReadInt());
            break;
        case 18:
            SetDisplayMonth(buf.ReadInt());
            break;
        case 19:
            SetDisplaySeason(buf.ReadInt());
            break;
        case 20:
            SetComputeParamValues(buf.ReadBool());
            break;
        case 21:
            SetDumpData(buf.ReadBool());
            break;
        case 22:
            SetDumpDebug(buf.ReadBool());
            break;
        }
    }

    public String toString(String indent)
    {
        String str = new String();
        str = str + intToString("dataYearBegin", dataYearBegin, indent) + "\n";
        str = str + boolToString("dataAnalysisYearRangeEnabled", dataAnalysisYearRangeEnabled, indent) + "\n";
        str = str + intToString("dataAnalysisYear1", dataAnalysisYear1, indent) + "\n";
        str = str + intToString("dataAnalysisYear2", dataAnalysisYear2, indent) + "\n";
        str = str + boolToString("ensemble", ensemble, indent) + "\n";
        str = str + intToString("numEnsembles", numEnsembles, indent) + "\n";
        str = str + doubleToString("dataScaling", dataScaling, indent) + "\n";
        str = str + indent + "extremeMethod = ";
        if(extremeMethod == EXTREMETYPE_MINIMA)
            str = str + "EXTREMETYPE_MINIMA";
        if(extremeMethod == EXTREMETYPE_MAXIMA)
            str = str + "EXTREMETYPE_MAXIMA";
        str = str + "\n";
        str = str + indent + "optimizationMethod = ";
        if(optimizationMethod == OPTIMIZATIONTYPE_NELDER_MEAD)
            str = str + "OPTIMIZATIONTYPE_NELDER_MEAD";
        if(optimizationMethod == OPTIMIZATIONTYPE_BFGS)
            str = str + "OPTIMIZATIONTYPE_BFGS";
        str = str + "\n";
        str = str + indent + "aggregation = ";
        if(aggregation == AGGREGATIONTYPE_ANNUAL)
            str = str + "AGGREGATIONTYPE_ANNUAL";
        if(aggregation == AGGREGATIONTYPE_SEASONAL)
            str = str + "AGGREGATIONTYPE_SEASONAL";
        if(aggregation == AGGREGATIONTYPE_MONTHLY)
            str = str + "AGGREGATIONTYPE_MONTHLY";
        str = str + "\n";
        str = str + boolToString("covariateModelScale", covariateModelScale, indent) + "\n";
        str = str + boolToString("covariateModelLocation", covariateModelLocation, indent) + "\n";
        str = str + boolToString("covariateModelShape", covariateModelShape, indent) + "\n";
        str = str + boolToString("computeReturnValues", computeReturnValues, indent) + "\n";
        str = str + intVectorToString("returnValues", returnValues, indent) + "\n";
        str = str + boolToString("computeRVDifferences", computeRVDifferences, indent) + "\n";
        str = str + intToString("rvDifference1", rvDifference1, indent) + "\n";
        str = str + intToString("rvDifference2", rvDifference2, indent) + "\n";
        str = str + indent + "displayMonth = ";
        if(displayMonth == MONTHTYPE_JANUARY)
            str = str + "MONTHTYPE_JANUARY";
        if(displayMonth == MONTHTYPE_FEBRUARY)
            str = str + "MONTHTYPE_FEBRUARY";
        if(displayMonth == MONTHTYPE_MARCH)
            str = str + "MONTHTYPE_MARCH";
        if(displayMonth == MONTHTYPE_APRIL)
            str = str + "MONTHTYPE_APRIL";
        if(displayMonth == MONTHTYPE_MAY)
            str = str + "MONTHTYPE_MAY";
        if(displayMonth == MONTHTYPE_JUNE)
            str = str + "MONTHTYPE_JUNE";
        if(displayMonth == MONTHTYPE_JULY)
            str = str + "MONTHTYPE_JULY";
        if(displayMonth == MONTHTYPE_AUGUST)
            str = str + "MONTHTYPE_AUGUST";
        if(displayMonth == MONTHTYPE_SEPTEMBER)
            str = str + "MONTHTYPE_SEPTEMBER";
        if(displayMonth == MONTHTYPE_OCTOBER)
            str = str + "MONTHTYPE_OCTOBER";
        if(displayMonth == MONTHTYPE_NOVEMBER)
            str = str + "MONTHTYPE_NOVEMBER";
        if(displayMonth == MONTHTYPE_DECEMBER)
            str = str + "MONTHTYPE_DECEMBER";
        str = str + "\n";
        str = str + indent + "displaySeason = ";
        if(displaySeason == SEASONTYPE_WINTER)
            str = str + "SEASONTYPE_WINTER";
        if(displaySeason == SEASONTYPE_SPRING)
            str = str + "SEASONTYPE_SPRING";
        if(displaySeason == SEASONTYPE_SUMMER)
            str = str + "SEASONTYPE_SUMMER";
        if(displaySeason == SEASONTYPE_FALL)
            str = str + "SEASONTYPE_FALL";
        str = str + "\n";
        str = str + boolToString("computeParamValues", computeParamValues, indent) + "\n";
        str = str + boolToString("dumpData", dumpData, indent) + "\n";
        str = str + boolToString("dumpDebug", dumpDebug, indent) + "\n";
        return str;
    }


    // Attributes
    private int     dataYearBegin;
    private boolean dataAnalysisYearRangeEnabled;
    private int     dataAnalysisYear1;
    private int     dataAnalysisYear2;
    private boolean ensemble;
    private int     numEnsembles;
    private double  dataScaling;
    private int     extremeMethod;
    private int     optimizationMethod;
    private int     aggregation;
    private boolean covariateModelScale;
    private boolean covariateModelLocation;
    private boolean covariateModelShape;
    private boolean computeReturnValues;
    private Vector  returnValues; // vector of Integer objects
    private boolean computeRVDifferences;
    private int     rvDifference1;
    private int     rvDifference2;
    private int     displayMonth;
    private int     displaySeason;
    private boolean computeParamValues;
    private boolean dumpData;
    private boolean dumpDebug;
}

