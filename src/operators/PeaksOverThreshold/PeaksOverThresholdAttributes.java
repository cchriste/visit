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
import java.lang.Integer;
import java.util.Vector;

// ****************************************************************************
// Class: PeaksOverThresholdAttributes
//
// Purpose:
//    Attributes for PeaksOverThreshold operator
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

public class PeaksOverThresholdAttributes extends AttributeSubject implements Plugin
{
    private static int PeaksOverThresholdAttributes_numAdditionalAtts = 20;

    // Enum values
    public final static int AGGREGATIONTYPE_ANNUAL = 0;
    public final static int AGGREGATIONTYPE_SEASONAL = 1;
    public final static int AGGREGATIONTYPE_MONTHLY = 2;

    public final static int SEASONTYPE_WINTER = 0;
    public final static int SEASONTYPE_SPRING = 1;
    public final static int SEASONTYPE_SUMMER = 2;
    public final static int SEASONTYPE_FALL = 3;

    public final static int MONTHTYPE_JAN = 0;
    public final static int MONTHTYPE_FEB = 1;
    public final static int MONTHTYPE_MAR = 2;
    public final static int MONTHTYPE_APR = 3;
    public final static int MONTHTYPE_MAY = 4;
    public final static int MONTHTYPE_JUN = 5;
    public final static int MONTHTYPE_JUL = 6;
    public final static int MONTHTYPE_AUG = 7;
    public final static int MONTHTYPE_SEP = 8;
    public final static int MONTHTYPE_OCT = 9;
    public final static int MONTHTYPE_NOV = 10;
    public final static int MONTHTYPE_DEC = 11;


    public PeaksOverThresholdAttributes()
    {
        super(PeaksOverThresholdAttributes_numAdditionalAtts);

        dataYearBegin = 1;
        dataAnalysisYearRangeEnabled = false;
        dataAnalysisYearRange = new int[2];
        dataAnalysisYearRange[0] = 0;
        dataAnalysisYearRange[1] = 0;
        aggregation = AGGREGATIONTYPE_ANNUAL;
        annualPercentile = 0.9;
        seasonalPercentile = new double[4];
        seasonalPercentile[0] = 0.9;
        seasonalPercentile[1] = 0.9;
        seasonalPercentile[2] = 0.9;
        seasonalPercentile[3] = 0.9;
        monthlyPercentile = new double[12];
        monthlyPercentile[0] = 0.9;
        monthlyPercentile[1] = 0.9;
        monthlyPercentile[2] = 0.9;
        monthlyPercentile[3] = 0.9;
        monthlyPercentile[4] = 0.9;
        monthlyPercentile[5] = 0.9;
        monthlyPercentile[6] = 0.9;
        monthlyPercentile[7] = 0.9;
        monthlyPercentile[8] = 0.9;
        monthlyPercentile[9] = 0.9;
        monthlyPercentile[10] = 0.9;
        monthlyPercentile[11] = 0.9;
        displaySeason = SEASONTYPE_WINTER;
        displayMonth = MONTHTYPE_JAN;
        dataYearBegin = 0;
        cutoff = 0f;
        computeParamValues = false;
        computeCovariates = false;
        covariateReturnYears = new Vector();
        covariateReturnYears.addElement(new Integer(1));
        covariateModelLocation = false;
        covariateModelShape = false;
        covariateModelScale = false;
        computeRVDifferences = false;
        rvDifferences = new int[2];
        rvDifferences[0] = 0;
        rvDifferences[1] = 0;
        dataScaling = 86500;
        dumpData = false;
    }

    public PeaksOverThresholdAttributes(int nMoreFields)
    {
        super(PeaksOverThresholdAttributes_numAdditionalAtts + nMoreFields);

        dataYearBegin = 1;
        dataAnalysisYearRangeEnabled = false;
        dataAnalysisYearRange = new int[2];
        dataAnalysisYearRange[0] = 0;
        dataAnalysisYearRange[1] = 0;
        aggregation = AGGREGATIONTYPE_ANNUAL;
        annualPercentile = 0.9;
        seasonalPercentile = new double[4];
        seasonalPercentile[0] = 0.9;
        seasonalPercentile[1] = 0.9;
        seasonalPercentile[2] = 0.9;
        seasonalPercentile[3] = 0.9;
        monthlyPercentile = new double[12];
        monthlyPercentile[0] = 0.9;
        monthlyPercentile[1] = 0.9;
        monthlyPercentile[2] = 0.9;
        monthlyPercentile[3] = 0.9;
        monthlyPercentile[4] = 0.9;
        monthlyPercentile[5] = 0.9;
        monthlyPercentile[6] = 0.9;
        monthlyPercentile[7] = 0.9;
        monthlyPercentile[8] = 0.9;
        monthlyPercentile[9] = 0.9;
        monthlyPercentile[10] = 0.9;
        monthlyPercentile[11] = 0.9;
        displaySeason = SEASONTYPE_WINTER;
        displayMonth = MONTHTYPE_JAN;
        dataYearBegin = 0;
        cutoff = 0f;
        computeParamValues = false;
        computeCovariates = false;
        covariateReturnYears = new Vector();
        covariateReturnYears.addElement(new Integer(1));
        covariateModelLocation = false;
        covariateModelShape = false;
        covariateModelScale = false;
        computeRVDifferences = false;
        rvDifferences = new int[2];
        rvDifferences[0] = 0;
        rvDifferences[1] = 0;
        dataScaling = 86500;
        dumpData = false;
    }

    public PeaksOverThresholdAttributes(PeaksOverThresholdAttributes obj)
    {
        super(PeaksOverThresholdAttributes_numAdditionalAtts);

        int i;

        dataYearBegin = obj.dataYearBegin;
        dataAnalysisYearRangeEnabled = obj.dataAnalysisYearRangeEnabled;
        dataAnalysisYearRange = new int[2];
        dataAnalysisYearRange[0] = obj.dataAnalysisYearRange[0];
        dataAnalysisYearRange[1] = obj.dataAnalysisYearRange[1];

        aggregation = obj.aggregation;
        annualPercentile = obj.annualPercentile;
        seasonalPercentile = new double[4];
        for(i = 0; i < obj.seasonalPercentile.length; ++i)
            seasonalPercentile[i] = obj.seasonalPercentile[i];

        monthlyPercentile = new double[12];
        for(i = 0; i < obj.monthlyPercentile.length; ++i)
            monthlyPercentile[i] = obj.monthlyPercentile[i];

        displaySeason = obj.displaySeason;
        displayMonth = obj.displayMonth;
        dataYearBegin = obj.dataYearBegin;
        cutoff = obj.cutoff;
        computeParamValues = obj.computeParamValues;
        computeCovariates = obj.computeCovariates;
        covariateReturnYears = new Vector();
        for(i = 0; i < obj.covariateReturnYears.size(); ++i)
        {
            Integer iv = (Integer)obj.covariateReturnYears.elementAt(i);
            covariateReturnYears.addElement(new Integer(iv.intValue()));
        }
        covariateModelLocation = obj.covariateModelLocation;
        covariateModelShape = obj.covariateModelShape;
        covariateModelScale = obj.covariateModelScale;
        computeRVDifferences = obj.computeRVDifferences;
        rvDifferences = new int[2];
        rvDifferences[0] = obj.rvDifferences[0];
        rvDifferences[1] = obj.rvDifferences[1];

        dataScaling = obj.dataScaling;
        dumpData = obj.dumpData;

        SelectAll();
    }

    public int Offset()
    {
        return super.Offset() + super.GetNumAdditionalAttributes();
    }

    public int GetNumAdditionalAttributes()
    {
        return PeaksOverThresholdAttributes_numAdditionalAtts;
    }

    public boolean equals(PeaksOverThresholdAttributes obj)
    {
        int i;

        // Compare the dataAnalysisYearRange arrays.
        boolean dataAnalysisYearRange_equal = true;
        for(i = 0; i < 2 && dataAnalysisYearRange_equal; ++i)
            dataAnalysisYearRange_equal = (dataAnalysisYearRange[i] == obj.dataAnalysisYearRange[i]);

        // Compare the seasonalPercentile arrays.
        boolean seasonalPercentile_equal = true;
        for(i = 0; i < 4 && seasonalPercentile_equal; ++i)
            seasonalPercentile_equal = (seasonalPercentile[i] == obj.seasonalPercentile[i]);

        // Compare the monthlyPercentile arrays.
        boolean monthlyPercentile_equal = true;
        for(i = 0; i < 12 && monthlyPercentile_equal; ++i)
            monthlyPercentile_equal = (monthlyPercentile[i] == obj.monthlyPercentile[i]);

        // Compare the elements in the covariateReturnYears vector.
        boolean covariateReturnYears_equal = (obj.covariateReturnYears.size() == covariateReturnYears.size());
        for(i = 0; (i < covariateReturnYears.size()) && covariateReturnYears_equal; ++i)
        {
            // Make references to Integer from Object.
            Integer covariateReturnYears1 = (Integer)covariateReturnYears.elementAt(i);
            Integer covariateReturnYears2 = (Integer)obj.covariateReturnYears.elementAt(i);
            covariateReturnYears_equal = covariateReturnYears1.equals(covariateReturnYears2);
        }
        // Compare the rvDifferences arrays.
        boolean rvDifferences_equal = true;
        for(i = 0; i < 2 && rvDifferences_equal; ++i)
            rvDifferences_equal = (rvDifferences[i] == obj.rvDifferences[i]);

        // Create the return value
        return ((dataYearBegin == obj.dataYearBegin) &&
                (dataAnalysisYearRangeEnabled == obj.dataAnalysisYearRangeEnabled) &&
                dataAnalysisYearRange_equal &&
                (aggregation == obj.aggregation) &&
                (annualPercentile == obj.annualPercentile) &&
                seasonalPercentile_equal &&
                monthlyPercentile_equal &&
                (displaySeason == obj.displaySeason) &&
                (displayMonth == obj.displayMonth) &&
                (dataYearBegin == obj.dataYearBegin) &&
                (cutoff == obj.cutoff) &&
                (computeParamValues == obj.computeParamValues) &&
                (computeCovariates == obj.computeCovariates) &&
                covariateReturnYears_equal &&
                (covariateModelLocation == obj.covariateModelLocation) &&
                (covariateModelShape == obj.covariateModelShape) &&
                (covariateModelScale == obj.covariateModelScale) &&
                (computeRVDifferences == obj.computeRVDifferences) &&
                rvDifferences_equal &&
                (dataScaling == obj.dataScaling) &&
                (dumpData == obj.dumpData));
    }

    public String GetName() { return "PeaksOverThreshold"; }
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

    public void SetDataAnalysisYearRange(int[] dataAnalysisYearRange_)
    {
        dataAnalysisYearRange[0] = dataAnalysisYearRange_[0];
        dataAnalysisYearRange[1] = dataAnalysisYearRange_[1];
        Select(2);
    }

    public void SetDataAnalysisYearRange(int e0, int e1)
    {
        dataAnalysisYearRange[0] = e0;
        dataAnalysisYearRange[1] = e1;
        Select(2);
    }

    public void SetAggregation(int aggregation_)
    {
        aggregation = aggregation_;
        Select(3);
    }

    public void SetAnnualPercentile(double annualPercentile_)
    {
        annualPercentile = annualPercentile_;
        Select(4);
    }

    public void SetSeasonalPercentile(double[] seasonalPercentile_)
    {
        seasonalPercentile[0] = seasonalPercentile_[0];
        seasonalPercentile[1] = seasonalPercentile_[1];
        seasonalPercentile[2] = seasonalPercentile_[2];
        seasonalPercentile[3] = seasonalPercentile_[3];
        Select(5);
    }

    public void SetSeasonalPercentile(double e0, double e1, double e2, double e3)
    {
        seasonalPercentile[0] = e0;
        seasonalPercentile[1] = e1;
        seasonalPercentile[2] = e2;
        seasonalPercentile[3] = e3;
        Select(5);
    }

    public void SetMonthlyPercentile(double[] monthlyPercentile_)
    {
        for(int i = 0; i < 12; ++i)
             monthlyPercentile[i] = monthlyPercentile_[i];
        Select(6);
    }

    public void SetDisplaySeason(int displaySeason_)
    {
        displaySeason = displaySeason_;
        Select(7);
    }

    public void SetDisplayMonth(int displayMonth_)
    {
        displayMonth = displayMonth_;
        Select(8);
    }

    public void SetDataYearBegin(int dataYearBegin_)
    {
        dataYearBegin = dataYearBegin_;
        Select(6);
    }

    public void SetCutoff(float cutoff_)
    {
        cutoff = cutoff_;
        Select(9);
    }

    public void SetComputeParamValues(boolean computeParamValues_)
    {
        computeParamValues = computeParamValues_;
        Select(10);
    }

    public void SetComputeCovariates(boolean computeCovariates_)
    {
        computeCovariates = computeCovariates_;
        Select(11);
    }

    public void SetCovariateReturnYears(Vector covariateReturnYears_)
    {
        covariateReturnYears = covariateReturnYears_;
        Select(12);
    }

    public void SetCovariateModelLocation(boolean covariateModelLocation_)
    {
        covariateModelLocation = covariateModelLocation_;
        Select(13);
    }

    public void SetCovariateModelShape(boolean covariateModelShape_)
    {
        covariateModelShape = covariateModelShape_;
        Select(14);
    }

    public void SetCovariateModelScale(boolean covariateModelScale_)
    {
        covariateModelScale = covariateModelScale_;
        Select(15);
    }

    public void SetComputeRVDifferences(boolean computeRVDifferences_)
    {
        computeRVDifferences = computeRVDifferences_;
        Select(16);
    }

    public void SetRvDifferences(int[] rvDifferences_)
    {
        rvDifferences[0] = rvDifferences_[0];
        rvDifferences[1] = rvDifferences_[1];
        Select(17);
    }

    public void SetRvDifferences(int e0, int e1)
    {
        rvDifferences[0] = e0;
        rvDifferences[1] = e1;
        Select(17);
    }

    public void SetDataScaling(double dataScaling_)
    {
        dataScaling = dataScaling_;
        Select(18);
    }

    public void SetDumpData(boolean dumpData_)
    {
        dumpData = dumpData_;
        Select(19);
    }

    // Property getting methods
    public int      GetDataYearBegin() { return dataYearBegin; }
    public boolean  GetDataAnalysisYearRangeEnabled() { return dataAnalysisYearRangeEnabled; }
    public int[]    GetDataAnalysisYearRange() { return dataAnalysisYearRange; }
    public int      GetAggregation() { return aggregation; }
    public double   GetAnnualPercentile() { return annualPercentile; }
    public double[] GetSeasonalPercentile() { return seasonalPercentile; }
    public double[] GetMonthlyPercentile() { return monthlyPercentile; }
    public int      GetDisplaySeason() { return displaySeason; }
    public int      GetDisplayMonth() { return displayMonth; }
    public int      GetDataYearBegin() { return dataYearBegin; }
    public float    GetCutoff() { return cutoff; }
    public boolean  GetComputeParamValues() { return computeParamValues; }
    public boolean  GetComputeCovariates() { return computeCovariates; }
    public Vector   GetCovariateReturnYears() { return covariateReturnYears; }
    public boolean  GetCovariateModelLocation() { return covariateModelLocation; }
    public boolean  GetCovariateModelShape() { return covariateModelShape; }
    public boolean  GetCovariateModelScale() { return covariateModelScale; }
    public boolean  GetComputeRVDifferences() { return computeRVDifferences; }
    public int[]    GetRvDifferences() { return rvDifferences; }
    public double   GetDataScaling() { return dataScaling; }
    public boolean  GetDumpData() { return dumpData; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteInt(dataYearBegin);
        if(WriteSelect(1, buf))
            buf.WriteBool(dataAnalysisYearRangeEnabled);
        if(WriteSelect(2, buf))
            buf.WriteIntArray(dataAnalysisYearRange);
        if(WriteSelect(3, buf))
            buf.WriteInt(aggregation);
        if(WriteSelect(4, buf))
            buf.WriteDouble(annualPercentile);
        if(WriteSelect(5, buf))
            buf.WriteDoubleArray(seasonalPercentile);
        if(WriteSelect(6, buf))
            buf.WriteDoubleArray(monthlyPercentile);
        if(WriteSelect(8, buf))
            buf.WriteInt(displaySeason);
        if(WriteSelect(9, buf))
            buf.WriteInt(displayMonth);
        if(WriteSelect(10, buf))
            buf.WriteFloat(cutoff);
        if(WriteSelect(10, buf))
            buf.WriteBool(computeParamValues);
        if(WriteSelect(11, buf))
            buf.WriteBool(computeCovariates);
        if(WriteSelect(12, buf))
            buf.WriteIntVector(covariateReturnYears);
        if(WriteSelect(13, buf))
            buf.WriteBool(covariateModelLocation);
        if(WriteSelect(14, buf))
            buf.WriteBool(covariateModelShape);
        if(WriteSelect(15, buf))
            buf.WriteBool(covariateModelScale);
        if(WriteSelect(16, buf))
            buf.WriteBool(computeRVDifferences);
        if(WriteSelect(17, buf))
            buf.WriteIntArray(rvDifferences);
        if(WriteSelect(18, buf))
            buf.WriteDouble(dataScaling);
        if(WriteSelect(19, buf))
            buf.WriteBool(dumpData);
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
            SetDataAnalysisYearRange(buf.ReadIntArray());
            break;
        case 3:
            SetAggregation(buf.ReadInt());
            break;
        case 4:
            SetAnnualPercentile(buf.ReadDouble());
            break;
        case 5:
            SetSeasonalPercentile(buf.ReadDoubleArray());
            break;
        case 6:
            SetMonthlyPercentile(buf.ReadDoubleArray());
            break;
        case 8:
            SetDisplaySeason(buf.ReadInt());
            break;
        case 9:
            SetDisplayMonth(buf.ReadInt());
            break;
        case 10:
            SetCutoff(buf.ReadFloat());
            break;
        case 10:
            SetComputeParamValues(buf.ReadBool());
            break;
        case 11:
            SetComputeCovariates(buf.ReadBool());
            break;
        case 12:
            SetCovariateReturnYears(buf.ReadIntVector());
            break;
        case 13:
            SetCovariateModelLocation(buf.ReadBool());
            break;
        case 14:
            SetCovariateModelShape(buf.ReadBool());
            break;
        case 15:
            SetCovariateModelScale(buf.ReadBool());
            break;
        case 16:
            SetComputeRVDifferences(buf.ReadBool());
            break;
        case 17:
            SetRvDifferences(buf.ReadIntArray());
            break;
        case 18:
            SetDataScaling(buf.ReadDouble());
            break;
        case 19:
            SetDumpData(buf.ReadBool());
            break;
        }
    }

    public String toString(String indent)
    {
        String str = new String();
        str = str + intToString("dataYearBegin", dataYearBegin, indent) + "\n";
        str = str + boolToString("dataAnalysisYearRangeEnabled", dataAnalysisYearRangeEnabled, indent) + "\n";
        str = str + intArrayToString("dataAnalysisYearRange", dataAnalysisYearRange, indent) + "\n";
        str = str + indent + "aggregation = ";
        if(aggregation == AGGREGATIONTYPE_ANNUAL)
            str = str + "AGGREGATIONTYPE_ANNUAL";
        if(aggregation == AGGREGATIONTYPE_SEASONAL)
            str = str + "AGGREGATIONTYPE_SEASONAL";
        if(aggregation == AGGREGATIONTYPE_MONTHLY)
            str = str + "AGGREGATIONTYPE_MONTHLY";
        str = str + "\n";
        str = str + doubleToString("annualPercentile", annualPercentile, indent) + "\n";
        str = str + doubleArrayToString("seasonalPercentile", seasonalPercentile, indent) + "\n";
        str = str + doubleArrayToString("monthlyPercentile", monthlyPercentile, indent) + "\n";
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
        str = str + indent + "displayMonth = ";
        if(displayMonth == MONTHTYPE_JAN)
            str = str + "MONTHTYPE_JAN";
        if(displayMonth == MONTHTYPE_FEB)
            str = str + "MONTHTYPE_FEB";
        if(displayMonth == MONTHTYPE_MAR)
            str = str + "MONTHTYPE_MAR";
        if(displayMonth == MONTHTYPE_APR)
            str = str + "MONTHTYPE_APR";
        if(displayMonth == MONTHTYPE_MAY)
            str = str + "MONTHTYPE_MAY";
        if(displayMonth == MONTHTYPE_JUN)
            str = str + "MONTHTYPE_JUN";
        if(displayMonth == MONTHTYPE_JUL)
            str = str + "MONTHTYPE_JUL";
        if(displayMonth == MONTHTYPE_AUG)
            str = str + "MONTHTYPE_AUG";
        if(displayMonth == MONTHTYPE_SEP)
            str = str + "MONTHTYPE_SEP";
        if(displayMonth == MONTHTYPE_OCT)
            str = str + "MONTHTYPE_OCT";
        if(displayMonth == MONTHTYPE_NOV)
            str = str + "MONTHTYPE_NOV";
        if(displayMonth == MONTHTYPE_DEC)
            str = str + "MONTHTYPE_DEC";
        str = str + "\n";
        str = str + intToString("dataYearBegin", dataYearBegin, indent) + "\n";
        str = str + floatToString("cutoff", cutoff, indent) + "\n";
        str = str + boolToString("computeParamValues", computeParamValues, indent) + "\n";
        str = str + boolToString("computeCovariates", computeCovariates, indent) + "\n";
        str = str + intVectorToString("covariateReturnYears", covariateReturnYears, indent) + "\n";
        str = str + boolToString("covariateModelLocation", covariateModelLocation, indent) + "\n";
        str = str + boolToString("covariateModelShape", covariateModelShape, indent) + "\n";
        str = str + boolToString("covariateModelScale", covariateModelScale, indent) + "\n";
        str = str + boolToString("computeRVDifferences", computeRVDifferences, indent) + "\n";
        str = str + intArrayToString("rvDifferences", rvDifferences, indent) + "\n";
        str = str + doubleToString("dataScaling", dataScaling, indent) + "\n";
        str = str + boolToString("dumpData", dumpData, indent) + "\n";
        return str;
    }


    // Attributes
    private int      dataYearBegin;
    private boolean  dataAnalysisYearRangeEnabled;
    private int[]    dataAnalysisYearRange;
    private int      aggregation;
    private double   annualPercentile;
    private double[] seasonalPercentile;
    private double[] monthlyPercentile;
    private int      displaySeason;
    private int      displayMonth;
    private int      dataYearBegin;
    private float    cutoff;
    private boolean  computeParamValues;
    private boolean  computeCovariates;
    private Vector   covariateReturnYears; // vector of Integer objects
    private boolean  covariateModelLocation;
    private boolean  covariateModelShape;
    private boolean  covariateModelScale;
    private boolean  computeRVDifferences;
    private int[]    rvDifferences;
    private double   dataScaling;
    private boolean  dumpData;
}

