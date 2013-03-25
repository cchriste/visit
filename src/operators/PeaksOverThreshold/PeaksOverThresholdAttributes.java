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
    private static int PeaksOverThresholdAttributes_numAdditionalAtts = 30;

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

    public final static int CUTOFFMODETYPE_UPPER_TAIL = 0;
    public final static int CUTOFFMODETYPE_LOWER_TAIL = 1;

    public final static int OPTIMIZATIONTYPE_NELDER_MEAD = 0;
    public final static int OPTIMIZATIONTYPE_BFGS = 1;


    public PeaksOverThresholdAttributes()
    {
        super(PeaksOverThresholdAttributes_numAdditionalAtts);

        dataYearBegin = 1;
        dataAnalysisYearRangeEnabled = false;
        dataAnalysisYear1 = 0;
        dataAnalysisYear2 = 0;
        ensemble = false;
        numEnsembles = 1;
        cutoff = 0f;
        cutoffMode = CUTOFFMODETYPE_UPPER_TAIL;
        noConsecutiveDay = false;
        optimizationMethod = OPTIMIZATIONTYPE_NELDER_MEAD;
        dataScaling = 1;
        aggregation = AGGREGATIONTYPE_ANNUAL;
        annualPercentile = 0.95;
        seasonalPercentile = new double[4];
        seasonalPercentile[0] = 0.95;
        seasonalPercentile[1] = 0.95;
        seasonalPercentile[2] = 0.95;
        seasonalPercentile[3] = 0.95;
        monthlyPercentile = new double[12];
        monthlyPercentile[0] = 0.95;
        monthlyPercentile[1] = 0.95;
        monthlyPercentile[2] = 0.95;
        monthlyPercentile[3] = 0.95;
        monthlyPercentile[4] = 0.95;
        monthlyPercentile[5] = 0.95;
        monthlyPercentile[6] = 0.95;
        monthlyPercentile[7] = 0.95;
        monthlyPercentile[8] = 0.95;
        monthlyPercentile[9] = 0.95;
        monthlyPercentile[10] = 0.95;
        monthlyPercentile[11] = 0.95;
        daysPerYear = 365;
        daysPerMonth = new int[12];
        daysPerMonth[0] = 31;
        daysPerMonth[1] = 28;
        daysPerMonth[2] = 31;
        daysPerMonth[3] = 30;
        daysPerMonth[4] = 31;
        daysPerMonth[5] = 30;
        daysPerMonth[6] = 31;
        daysPerMonth[7] = 31;
        daysPerMonth[8] = 30;
        daysPerMonth[9] = 31;
        daysPerMonth[10] = 30;
        daysPerMonth[11] = 31;
        covariateModelScale = false;
        covariateModelLocation = false;
        covariateModelShape = false;
        computeCovariates = false;
        covariateReturnYears = new Vector();
        covariateReturnYears.addElement(new Integer(1));
        computeRVDifferences = false;
        rvDifference1 = 0;
        rvDifference2 = 0;
        computeParamValues = false;
        displaySeason = SEASONTYPE_WINTER;
        displayMonth = MONTHTYPE_JAN;
        dumpData = true;
        dumpDebug = false;
    }

    public PeaksOverThresholdAttributes(int nMoreFields)
    {
        super(PeaksOverThresholdAttributes_numAdditionalAtts + nMoreFields);

        dataYearBegin = 1;
        dataAnalysisYearRangeEnabled = false;
        dataAnalysisYear1 = 0;
        dataAnalysisYear2 = 0;
        ensemble = false;
        numEnsembles = 1;
        cutoff = 0f;
        cutoffMode = CUTOFFMODETYPE_UPPER_TAIL;
        noConsecutiveDay = false;
        optimizationMethod = OPTIMIZATIONTYPE_NELDER_MEAD;
        dataScaling = 1;
        aggregation = AGGREGATIONTYPE_ANNUAL;
        annualPercentile = 0.95;
        seasonalPercentile = new double[4];
        seasonalPercentile[0] = 0.95;
        seasonalPercentile[1] = 0.95;
        seasonalPercentile[2] = 0.95;
        seasonalPercentile[3] = 0.95;
        monthlyPercentile = new double[12];
        monthlyPercentile[0] = 0.95;
        monthlyPercentile[1] = 0.95;
        monthlyPercentile[2] = 0.95;
        monthlyPercentile[3] = 0.95;
        monthlyPercentile[4] = 0.95;
        monthlyPercentile[5] = 0.95;
        monthlyPercentile[6] = 0.95;
        monthlyPercentile[7] = 0.95;
        monthlyPercentile[8] = 0.95;
        monthlyPercentile[9] = 0.95;
        monthlyPercentile[10] = 0.95;
        monthlyPercentile[11] = 0.95;
        daysPerYear = 365;
        daysPerMonth = new int[12];
        daysPerMonth[0] = 31;
        daysPerMonth[1] = 28;
        daysPerMonth[2] = 31;
        daysPerMonth[3] = 30;
        daysPerMonth[4] = 31;
        daysPerMonth[5] = 30;
        daysPerMonth[6] = 31;
        daysPerMonth[7] = 31;
        daysPerMonth[8] = 30;
        daysPerMonth[9] = 31;
        daysPerMonth[10] = 30;
        daysPerMonth[11] = 31;
        covariateModelScale = false;
        covariateModelLocation = false;
        covariateModelShape = false;
        computeCovariates = false;
        covariateReturnYears = new Vector();
        covariateReturnYears.addElement(new Integer(1));
        computeRVDifferences = false;
        rvDifference1 = 0;
        rvDifference2 = 0;
        computeParamValues = false;
        displaySeason = SEASONTYPE_WINTER;
        displayMonth = MONTHTYPE_JAN;
        dumpData = true;
        dumpDebug = false;
    }

    public PeaksOverThresholdAttributes(PeaksOverThresholdAttributes obj)
    {
        super(PeaksOverThresholdAttributes_numAdditionalAtts);

        int i;

        dataYearBegin = obj.dataYearBegin;
        dataAnalysisYearRangeEnabled = obj.dataAnalysisYearRangeEnabled;
        dataAnalysisYear1 = obj.dataAnalysisYear1;
        dataAnalysisYear2 = obj.dataAnalysisYear2;
        ensemble = obj.ensemble;
        numEnsembles = obj.numEnsembles;
        cutoff = obj.cutoff;
        cutoffMode = obj.cutoffMode;
        noConsecutiveDay = obj.noConsecutiveDay;
        optimizationMethod = obj.optimizationMethod;
        dataScaling = obj.dataScaling;
        aggregation = obj.aggregation;
        annualPercentile = obj.annualPercentile;
        seasonalPercentile = new double[4];
        for(i = 0; i < obj.seasonalPercentile.length; ++i)
            seasonalPercentile[i] = obj.seasonalPercentile[i];

        monthlyPercentile = new double[12];
        for(i = 0; i < obj.monthlyPercentile.length; ++i)
            monthlyPercentile[i] = obj.monthlyPercentile[i];

        daysPerYear = obj.daysPerYear;
        daysPerMonth = new int[12];
        for(i = 0; i < obj.daysPerMonth.length; ++i)
            daysPerMonth[i] = obj.daysPerMonth[i];

        covariateModelScale = obj.covariateModelScale;
        covariateModelLocation = obj.covariateModelLocation;
        covariateModelShape = obj.covariateModelShape;
        computeCovariates = obj.computeCovariates;
        covariateReturnYears = new Vector();
        for(i = 0; i < obj.covariateReturnYears.size(); ++i)
        {
            Integer iv = (Integer)obj.covariateReturnYears.elementAt(i);
            covariateReturnYears.addElement(new Integer(iv.intValue()));
        }
        computeRVDifferences = obj.computeRVDifferences;
        rvDifference1 = obj.rvDifference1;
        rvDifference2 = obj.rvDifference2;
        computeParamValues = obj.computeParamValues;
        displaySeason = obj.displaySeason;
        displayMonth = obj.displayMonth;
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
        return PeaksOverThresholdAttributes_numAdditionalAtts;
    }

    public boolean equals(PeaksOverThresholdAttributes obj)
    {
        int i;

        // Compare the seasonalPercentile arrays.
        boolean seasonalPercentile_equal = true;
        for(i = 0; i < 4 && seasonalPercentile_equal; ++i)
            seasonalPercentile_equal = (seasonalPercentile[i] == obj.seasonalPercentile[i]);

        // Compare the monthlyPercentile arrays.
        boolean monthlyPercentile_equal = true;
        for(i = 0; i < 12 && monthlyPercentile_equal; ++i)
            monthlyPercentile_equal = (monthlyPercentile[i] == obj.monthlyPercentile[i]);

        // Compare the daysPerMonth arrays.
        boolean daysPerMonth_equal = true;
        for(i = 0; i < 12 && daysPerMonth_equal; ++i)
            daysPerMonth_equal = (daysPerMonth[i] == obj.daysPerMonth[i]);

        // Compare the elements in the covariateReturnYears vector.
        boolean covariateReturnYears_equal = (obj.covariateReturnYears.size() == covariateReturnYears.size());
        for(i = 0; (i < covariateReturnYears.size()) && covariateReturnYears_equal; ++i)
        {
            // Make references to Integer from Object.
            Integer covariateReturnYears1 = (Integer)covariateReturnYears.elementAt(i);
            Integer covariateReturnYears2 = (Integer)obj.covariateReturnYears.elementAt(i);
            covariateReturnYears_equal = covariateReturnYears1.equals(covariateReturnYears2);
        }
        // Create the return value
        return ((dataYearBegin == obj.dataYearBegin) &&
                (dataAnalysisYearRangeEnabled == obj.dataAnalysisYearRangeEnabled) &&
                (dataAnalysisYear1 == obj.dataAnalysisYear1) &&
                (dataAnalysisYear2 == obj.dataAnalysisYear2) &&
                (ensemble == obj.ensemble) &&
                (numEnsembles == obj.numEnsembles) &&
                (cutoff == obj.cutoff) &&
                (cutoffMode == obj.cutoffMode) &&
                (noConsecutiveDay == obj.noConsecutiveDay) &&
                (optimizationMethod == obj.optimizationMethod) &&
                (dataScaling == obj.dataScaling) &&
                (aggregation == obj.aggregation) &&
                (annualPercentile == obj.annualPercentile) &&
                seasonalPercentile_equal &&
                monthlyPercentile_equal &&
                (daysPerYear == obj.daysPerYear) &&
                daysPerMonth_equal &&
                (covariateModelScale == obj.covariateModelScale) &&
                (covariateModelLocation == obj.covariateModelLocation) &&
                (covariateModelShape == obj.covariateModelShape) &&
                (computeCovariates == obj.computeCovariates) &&
                covariateReturnYears_equal &&
                (computeRVDifferences == obj.computeRVDifferences) &&
                (rvDifference1 == obj.rvDifference1) &&
                (rvDifference2 == obj.rvDifference2) &&
                (computeParamValues == obj.computeParamValues) &&
                (displaySeason == obj.displaySeason) &&
                (displayMonth == obj.displayMonth) &&
                (dumpData == obj.dumpData) &&
                (dumpDebug == obj.dumpDebug));
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

    public void SetCutoff(float cutoff_)
    {
        cutoff = cutoff_;
        Select(6);
    }

    public void SetCutoffMode(int cutoffMode_)
    {
        cutoffMode = cutoffMode_;
        Select(7);
    }

    public void SetNoConsecutiveDay(boolean noConsecutiveDay_)
    {
        noConsecutiveDay = noConsecutiveDay_;
        Select(8);
    }

    public void SetOptimizationMethod(int optimizationMethod_)
    {
        optimizationMethod = optimizationMethod_;
        Select(9);
    }

    public void SetDataScaling(double dataScaling_)
    {
        dataScaling = dataScaling_;
        Select(10);
    }

    public void SetAggregation(int aggregation_)
    {
        aggregation = aggregation_;
        Select(11);
    }

    public void SetAnnualPercentile(double annualPercentile_)
    {
        annualPercentile = annualPercentile_;
        Select(12);
    }

    public void SetSeasonalPercentile(double[] seasonalPercentile_)
    {
        seasonalPercentile[0] = seasonalPercentile_[0];
        seasonalPercentile[1] = seasonalPercentile_[1];
        seasonalPercentile[2] = seasonalPercentile_[2];
        seasonalPercentile[3] = seasonalPercentile_[3];
        Select(13);
    }

    public void SetSeasonalPercentile(double e0, double e1, double e2, double e3)
    {
        seasonalPercentile[0] = e0;
        seasonalPercentile[1] = e1;
        seasonalPercentile[2] = e2;
        seasonalPercentile[3] = e3;
        Select(13);
    }

    public void SetMonthlyPercentile(double[] monthlyPercentile_)
    {
        for(int i = 0; i < 12; ++i)
             monthlyPercentile[i] = monthlyPercentile_[i];
        Select(14);
    }

    public void SetDaysPerYear(int daysPerYear_)
    {
        daysPerYear = daysPerYear_;
        Select(15);
    }

    public void SetDaysPerMonth(int[] daysPerMonth_)
    {
        for(int i = 0; i < 12; ++i)
             daysPerMonth[i] = daysPerMonth_[i];
        Select(16);
    }

    public void SetCovariateModelScale(boolean covariateModelScale_)
    {
        covariateModelScale = covariateModelScale_;
        Select(17);
    }

    public void SetCovariateModelLocation(boolean covariateModelLocation_)
    {
        covariateModelLocation = covariateModelLocation_;
        Select(18);
    }

    public void SetCovariateModelShape(boolean covariateModelShape_)
    {
        covariateModelShape = covariateModelShape_;
        Select(19);
    }

    public void SetComputeCovariates(boolean computeCovariates_)
    {
        computeCovariates = computeCovariates_;
        Select(20);
    }

    public void SetCovariateReturnYears(Vector covariateReturnYears_)
    {
        covariateReturnYears = covariateReturnYears_;
        Select(21);
    }

    public void SetComputeRVDifferences(boolean computeRVDifferences_)
    {
        computeRVDifferences = computeRVDifferences_;
        Select(22);
    }

    public void SetRvDifference1(int rvDifference1_)
    {
        rvDifference1 = rvDifference1_;
        Select(23);
    }

    public void SetRvDifference2(int rvDifference2_)
    {
        rvDifference2 = rvDifference2_;
        Select(24);
    }

    public void SetComputeParamValues(boolean computeParamValues_)
    {
        computeParamValues = computeParamValues_;
        Select(25);
    }

    public void SetDisplaySeason(int displaySeason_)
    {
        displaySeason = displaySeason_;
        Select(26);
    }

    public void SetDisplayMonth(int displayMonth_)
    {
        displayMonth = displayMonth_;
        Select(27);
    }

    public void SetDumpData(boolean dumpData_)
    {
        dumpData = dumpData_;
        Select(28);
    }

    public void SetDumpDebug(boolean dumpDebug_)
    {
        dumpDebug = dumpDebug_;
        Select(29);
    }

    // Property getting methods
    public int      GetDataYearBegin() { return dataYearBegin; }
    public boolean  GetDataAnalysisYearRangeEnabled() { return dataAnalysisYearRangeEnabled; }
    public int      GetDataAnalysisYear1() { return dataAnalysisYear1; }
    public int      GetDataAnalysisYear2() { return dataAnalysisYear2; }
    public boolean  GetEnsemble() { return ensemble; }
    public int      GetNumEnsembles() { return numEnsembles; }
    public float    GetCutoff() { return cutoff; }
    public int      GetCutoffMode() { return cutoffMode; }
    public boolean  GetNoConsecutiveDay() { return noConsecutiveDay; }
    public int      GetOptimizationMethod() { return optimizationMethod; }
    public double   GetDataScaling() { return dataScaling; }
    public int      GetAggregation() { return aggregation; }
    public double   GetAnnualPercentile() { return annualPercentile; }
    public double[] GetSeasonalPercentile() { return seasonalPercentile; }
    public double[] GetMonthlyPercentile() { return monthlyPercentile; }
    public int      GetDaysPerYear() { return daysPerYear; }
    public int[]    GetDaysPerMonth() { return daysPerMonth; }
    public boolean  GetCovariateModelScale() { return covariateModelScale; }
    public boolean  GetCovariateModelLocation() { return covariateModelLocation; }
    public boolean  GetCovariateModelShape() { return covariateModelShape; }
    public boolean  GetComputeCovariates() { return computeCovariates; }
    public Vector   GetCovariateReturnYears() { return covariateReturnYears; }
    public boolean  GetComputeRVDifferences() { return computeRVDifferences; }
    public int      GetRvDifference1() { return rvDifference1; }
    public int      GetRvDifference2() { return rvDifference2; }
    public boolean  GetComputeParamValues() { return computeParamValues; }
    public int      GetDisplaySeason() { return displaySeason; }
    public int      GetDisplayMonth() { return displayMonth; }
    public boolean  GetDumpData() { return dumpData; }
    public boolean  GetDumpDebug() { return dumpDebug; }

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
            buf.WriteFloat(cutoff);
        if(WriteSelect(7, buf))
            buf.WriteInt(cutoffMode);
        if(WriteSelect(8, buf))
            buf.WriteBool(noConsecutiveDay);
        if(WriteSelect(9, buf))
            buf.WriteInt(optimizationMethod);
        if(WriteSelect(10, buf))
            buf.WriteDouble(dataScaling);
        if(WriteSelect(11, buf))
            buf.WriteInt(aggregation);
        if(WriteSelect(12, buf))
            buf.WriteDouble(annualPercentile);
        if(WriteSelect(13, buf))
            buf.WriteDoubleArray(seasonalPercentile);
        if(WriteSelect(14, buf))
            buf.WriteDoubleArray(monthlyPercentile);
        if(WriteSelect(15, buf))
            buf.WriteInt(daysPerYear);
        if(WriteSelect(16, buf))
            buf.WriteIntArray(daysPerMonth);
        if(WriteSelect(17, buf))
            buf.WriteBool(covariateModelScale);
        if(WriteSelect(18, buf))
            buf.WriteBool(covariateModelLocation);
        if(WriteSelect(19, buf))
            buf.WriteBool(covariateModelShape);
        if(WriteSelect(20, buf))
            buf.WriteBool(computeCovariates);
        if(WriteSelect(21, buf))
            buf.WriteIntVector(covariateReturnYears);
        if(WriteSelect(22, buf))
            buf.WriteBool(computeRVDifferences);
        if(WriteSelect(23, buf))
            buf.WriteInt(rvDifference1);
        if(WriteSelect(24, buf))
            buf.WriteInt(rvDifference2);
        if(WriteSelect(25, buf))
            buf.WriteBool(computeParamValues);
        if(WriteSelect(26, buf))
            buf.WriteInt(displaySeason);
        if(WriteSelect(27, buf))
            buf.WriteInt(displayMonth);
        if(WriteSelect(28, buf))
            buf.WriteBool(dumpData);
        if(WriteSelect(29, buf))
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
            SetCutoff(buf.ReadFloat());
            break;
        case 7:
            SetCutoffMode(buf.ReadInt());
            break;
        case 8:
            SetNoConsecutiveDay(buf.ReadBool());
            break;
        case 9:
            SetOptimizationMethod(buf.ReadInt());
            break;
        case 10:
            SetDataScaling(buf.ReadDouble());
            break;
        case 11:
            SetAggregation(buf.ReadInt());
            break;
        case 12:
            SetAnnualPercentile(buf.ReadDouble());
            break;
        case 13:
            SetSeasonalPercentile(buf.ReadDoubleArray());
            break;
        case 14:
            SetMonthlyPercentile(buf.ReadDoubleArray());
            break;
        case 15:
            SetDaysPerYear(buf.ReadInt());
            break;
        case 16:
            SetDaysPerMonth(buf.ReadIntArray());
            break;
        case 17:
            SetCovariateModelScale(buf.ReadBool());
            break;
        case 18:
            SetCovariateModelLocation(buf.ReadBool());
            break;
        case 19:
            SetCovariateModelShape(buf.ReadBool());
            break;
        case 20:
            SetComputeCovariates(buf.ReadBool());
            break;
        case 21:
            SetCovariateReturnYears(buf.ReadIntVector());
            break;
        case 22:
            SetComputeRVDifferences(buf.ReadBool());
            break;
        case 23:
            SetRvDifference1(buf.ReadInt());
            break;
        case 24:
            SetRvDifference2(buf.ReadInt());
            break;
        case 25:
            SetComputeParamValues(buf.ReadBool());
            break;
        case 26:
            SetDisplaySeason(buf.ReadInt());
            break;
        case 27:
            SetDisplayMonth(buf.ReadInt());
            break;
        case 28:
            SetDumpData(buf.ReadBool());
            break;
        case 29:
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
        str = str + floatToString("cutoff", cutoff, indent) + "\n";
        str = str + indent + "cutoffMode = ";
        if(cutoffMode == CUTOFFMODETYPE_UPPER_TAIL)
            str = str + "CUTOFFMODETYPE_UPPER_TAIL";
        if(cutoffMode == CUTOFFMODETYPE_LOWER_TAIL)
            str = str + "CUTOFFMODETYPE_LOWER_TAIL";
        str = str + "\n";
        str = str + boolToString("noConsecutiveDay", noConsecutiveDay, indent) + "\n";
        str = str + indent + "optimizationMethod = ";
        if(optimizationMethod == OPTIMIZATIONTYPE_NELDER_MEAD)
            str = str + "OPTIMIZATIONTYPE_NELDER_MEAD";
        if(optimizationMethod == OPTIMIZATIONTYPE_BFGS)
            str = str + "OPTIMIZATIONTYPE_BFGS";
        str = str + "\n";
        str = str + doubleToString("dataScaling", dataScaling, indent) + "\n";
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
        str = str + intToString("daysPerYear", daysPerYear, indent) + "\n";
        str = str + intArrayToString("daysPerMonth", daysPerMonth, indent) + "\n";
        str = str + boolToString("covariateModelScale", covariateModelScale, indent) + "\n";
        str = str + boolToString("covariateModelLocation", covariateModelLocation, indent) + "\n";
        str = str + boolToString("covariateModelShape", covariateModelShape, indent) + "\n";
        str = str + boolToString("computeCovariates", computeCovariates, indent) + "\n";
        str = str + intVectorToString("covariateReturnYears", covariateReturnYears, indent) + "\n";
        str = str + boolToString("computeRVDifferences", computeRVDifferences, indent) + "\n";
        str = str + intToString("rvDifference1", rvDifference1, indent) + "\n";
        str = str + intToString("rvDifference2", rvDifference2, indent) + "\n";
        str = str + boolToString("computeParamValues", computeParamValues, indent) + "\n";
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
        str = str + boolToString("dumpData", dumpData, indent) + "\n";
        str = str + boolToString("dumpDebug", dumpDebug, indent) + "\n";
        return str;
    }


    // Attributes
    private int      dataYearBegin;
    private boolean  dataAnalysisYearRangeEnabled;
    private int      dataAnalysisYear1;
    private int      dataAnalysisYear2;
    private boolean  ensemble;
    private int      numEnsembles;
    private float    cutoff;
    private int      cutoffMode;
    private boolean  noConsecutiveDay;
    private int      optimizationMethod;
    private double   dataScaling;
    private int      aggregation;
    private double   annualPercentile;
    private double[] seasonalPercentile;
    private double[] monthlyPercentile;
    private int      daysPerYear;
    private int[]    daysPerMonth;
    private boolean  covariateModelScale;
    private boolean  covariateModelLocation;
    private boolean  covariateModelShape;
    private boolean  computeCovariates;
    private Vector   covariateReturnYears; // vector of Integer objects
    private boolean  computeRVDifferences;
    private int      rvDifference1;
    private int      rvDifference2;
    private boolean  computeParamValues;
    private int      displaySeason;
    private int      displayMonth;
    private boolean  dumpData;
    private boolean  dumpDebug;
}

