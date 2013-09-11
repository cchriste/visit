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
// Class: FTLEAttributes
//
// Purpose:
//    Attributes for FTLE
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

public class FTLEAttributes extends AttributeSubject implements Plugin
{
    private static int FTLEAttributes_numAdditionalAtts = 39;

    // Enum values
    public final static int SOURCETYPE_NATIVERESOLUTIONOFMESH = 0;
    public final static int SOURCETYPE_REGULARGRID = 1;

    public final static int EXTENTS_FULL = 0;
    public final static int EXTENTS_SUBSET = 1;

    public final static int INTEGRATIONDIRECTION_FORWARD = 0;
    public final static int INTEGRATIONDIRECTION_BACKWARD = 1;
    public final static int INTEGRATIONDIRECTION_BOTH = 2;

    public final static int FIELDTYPE_DEFAULT = 0;
    public final static int FIELDTYPE_FLASHFIELD = 1;
    public final static int FIELDTYPE_M3DC12DFIELD = 2;
    public final static int FIELDTYPE_M3DC13DFIELD = 3;
    public final static int FIELDTYPE_NEK5000FIELD = 4;
    public final static int FIELDTYPE_NIMRODFIELD = 5;

    public final static int INTEGRATIONTYPE_EULER = 0;
    public final static int INTEGRATIONTYPE_LEAPFROG = 1;
    public final static int INTEGRATIONTYPE_DORMANDPRINCE = 2;
    public final static int INTEGRATIONTYPE_ADAMSBASHFORTH = 3;
    public final static int INTEGRATIONTYPE_RK4 = 4;
    public final static int INTEGRATIONTYPE_M3DC12DINTEGRATOR = 5;

    public final static int SIZETYPE_ABSOLUTE = 0;
    public final static int SIZETYPE_FRACTIONOFBBOX = 1;

    public final static int PARALLELIZATIONALGORITHMTYPE_LOADONDEMAND = 0;
    public final static int PARALLELIZATIONALGORITHMTYPE_PARALLELSTATICDOMAINS = 1;
    public final static int PARALLELIZATIONALGORITHMTYPE_MASTERSLAVE = 2;
    public final static int PARALLELIZATIONALGORITHMTYPE_VISITSELECTS = 3;

    public final static int TERMINATIONTYPE_TIME = 0;
    public final static int TERMINATIONTYPE_DISTANCE = 1;
    public final static int TERMINATIONTYPE_SIZE = 2;

    public final static int PATHLINESCMFE_CONN_CMFE = 0;
    public final static int PATHLINESCMFE_POS_CMFE = 1;


    public FTLEAttributes()
    {
        super(FTLEAttributes_numAdditionalAtts);

        sourceType = SOURCETYPE_NATIVERESOLUTIONOFMESH;
        Resolution = new int[3];
        Resolution[0] = 10;
        Resolution[1] = 10;
        Resolution[2] = 10;
        UseDataSetStart = EXTENTS_FULL;
        StartPosition = new double[3];
        StartPosition[0] = 0;
        StartPosition[1] = 0;
        StartPosition[2] = 0;
        UseDataSetEnd = EXTENTS_FULL;
        EndPosition = new double[3];
        EndPosition[0] = 1;
        EndPosition[1] = 1;
        EndPosition[2] = 1;
        integrationDirection = INTEGRATIONDIRECTION_FORWARD;
        maxSteps = 1000;
        terminationType = TERMINATIONTYPE_TIME;
        terminateBySize = false;
        termSize = 10;
        terminateByDistance = false;
        termDistance = 10;
        terminateByTime = false;
        termTime = 10;
        maxStepLength = 0.1;
        limitMaximumTimestep = false;
        maxTimeStep = 0.1;
        relTol = 0.0001;
        absTolSizeType = SIZETYPE_FRACTIONOFBBOX;
        absTolAbsolute = 1e-06;
        absTolBBox = 1e-06;
        fieldType = FIELDTYPE_DEFAULT;
        fieldConstant = 1;
        velocitySource = new double[3];
        velocitySource[0] = 0;
        velocitySource[1] = 0;
        velocitySource[2] = 0;
        integrationType = INTEGRATIONTYPE_DORMANDPRINCE;
        parallelizationAlgorithmType = PARALLELIZATIONALGORITHMTYPE_VISITSELECTS;
        maxProcessCount = 10;
        maxDomainCacheSize = 3;
        workGroupSize = 32;
        pathlines = false;
        pathlinesOverrideStartingTimeFlag = false;
        pathlinesOverrideStartingTime = 0;
        pathlinesCMFE = PATHLINESCMFE_POS_CMFE;
        forceNodeCenteredData = false;
        issueTerminationWarnings = true;
        issueStiffnessWarnings = true;
        issueCriticalPointsWarnings = true;
        criticalPointThreshold = 0.001;
    }

    public FTLEAttributes(int nMoreFields)
    {
        super(FTLEAttributes_numAdditionalAtts + nMoreFields);

        sourceType = SOURCETYPE_NATIVERESOLUTIONOFMESH;
        Resolution = new int[3];
        Resolution[0] = 10;
        Resolution[1] = 10;
        Resolution[2] = 10;
        UseDataSetStart = EXTENTS_FULL;
        StartPosition = new double[3];
        StartPosition[0] = 0;
        StartPosition[1] = 0;
        StartPosition[2] = 0;
        UseDataSetEnd = EXTENTS_FULL;
        EndPosition = new double[3];
        EndPosition[0] = 1;
        EndPosition[1] = 1;
        EndPosition[2] = 1;
        integrationDirection = INTEGRATIONDIRECTION_FORWARD;
        maxSteps = 1000;
        terminationType = TERMINATIONTYPE_TIME;
        terminateBySize = false;
        termSize = 10;
        terminateByDistance = false;
        termDistance = 10;
        terminateByTime = false;
        termTime = 10;
        maxStepLength = 0.1;
        limitMaximumTimestep = false;
        maxTimeStep = 0.1;
        relTol = 0.0001;
        absTolSizeType = SIZETYPE_FRACTIONOFBBOX;
        absTolAbsolute = 1e-06;
        absTolBBox = 1e-06;
        fieldType = FIELDTYPE_DEFAULT;
        fieldConstant = 1;
        velocitySource = new double[3];
        velocitySource[0] = 0;
        velocitySource[1] = 0;
        velocitySource[2] = 0;
        integrationType = INTEGRATIONTYPE_DORMANDPRINCE;
        parallelizationAlgorithmType = PARALLELIZATIONALGORITHMTYPE_VISITSELECTS;
        maxProcessCount = 10;
        maxDomainCacheSize = 3;
        workGroupSize = 32;
        pathlines = false;
        pathlinesOverrideStartingTimeFlag = false;
        pathlinesOverrideStartingTime = 0;
        pathlinesCMFE = PATHLINESCMFE_POS_CMFE;
        forceNodeCenteredData = false;
        issueTerminationWarnings = true;
        issueStiffnessWarnings = true;
        issueCriticalPointsWarnings = true;
        criticalPointThreshold = 0.001;
    }

    public FTLEAttributes(FTLEAttributes obj)
    {
        super(FTLEAttributes_numAdditionalAtts);

        int i;

        sourceType = obj.sourceType;
        Resolution = new int[3];
        Resolution[0] = obj.Resolution[0];
        Resolution[1] = obj.Resolution[1];
        Resolution[2] = obj.Resolution[2];

        UseDataSetStart = obj.UseDataSetStart;
        StartPosition = new double[3];
        StartPosition[0] = obj.StartPosition[0];
        StartPosition[1] = obj.StartPosition[1];
        StartPosition[2] = obj.StartPosition[2];

        UseDataSetEnd = obj.UseDataSetEnd;
        EndPosition = new double[3];
        EndPosition[0] = obj.EndPosition[0];
        EndPosition[1] = obj.EndPosition[1];
        EndPosition[2] = obj.EndPosition[2];

        integrationDirection = obj.integrationDirection;
        maxSteps = obj.maxSteps;
        terminationType = obj.terminationType;
        terminateBySize = obj.terminateBySize;
        termSize = obj.termSize;
        terminateByDistance = obj.terminateByDistance;
        termDistance = obj.termDistance;
        terminateByTime = obj.terminateByTime;
        termTime = obj.termTime;
        maxStepLength = obj.maxStepLength;
        limitMaximumTimestep = obj.limitMaximumTimestep;
        maxTimeStep = obj.maxTimeStep;
        relTol = obj.relTol;
        absTolSizeType = obj.absTolSizeType;
        absTolAbsolute = obj.absTolAbsolute;
        absTolBBox = obj.absTolBBox;
        fieldType = obj.fieldType;
        fieldConstant = obj.fieldConstant;
        velocitySource = new double[3];
        velocitySource[0] = obj.velocitySource[0];
        velocitySource[1] = obj.velocitySource[1];
        velocitySource[2] = obj.velocitySource[2];

        integrationType = obj.integrationType;
        parallelizationAlgorithmType = obj.parallelizationAlgorithmType;
        maxProcessCount = obj.maxProcessCount;
        maxDomainCacheSize = obj.maxDomainCacheSize;
        workGroupSize = obj.workGroupSize;
        pathlines = obj.pathlines;
        pathlinesOverrideStartingTimeFlag = obj.pathlinesOverrideStartingTimeFlag;
        pathlinesOverrideStartingTime = obj.pathlinesOverrideStartingTime;
        pathlinesCMFE = obj.pathlinesCMFE;
        forceNodeCenteredData = obj.forceNodeCenteredData;
        issueTerminationWarnings = obj.issueTerminationWarnings;
        issueStiffnessWarnings = obj.issueStiffnessWarnings;
        issueCriticalPointsWarnings = obj.issueCriticalPointsWarnings;
        criticalPointThreshold = obj.criticalPointThreshold;

        SelectAll();
    }

    public int Offset()
    {
        return super.Offset() + super.GetNumAdditionalAttributes();
    }

    public int GetNumAdditionalAttributes()
    {
        return FTLEAttributes_numAdditionalAtts;
    }

    public boolean equals(FTLEAttributes obj)
    {
        int i;

        // Compare the Resolution arrays.
        boolean Resolution_equal = true;
        for(i = 0; i < 3 && Resolution_equal; ++i)
            Resolution_equal = (Resolution[i] == obj.Resolution[i]);

        // Compare the StartPosition arrays.
        boolean StartPosition_equal = true;
        for(i = 0; i < 3 && StartPosition_equal; ++i)
            StartPosition_equal = (StartPosition[i] == obj.StartPosition[i]);

        // Compare the EndPosition arrays.
        boolean EndPosition_equal = true;
        for(i = 0; i < 3 && EndPosition_equal; ++i)
            EndPosition_equal = (EndPosition[i] == obj.EndPosition[i]);

        // Compare the velocitySource arrays.
        boolean velocitySource_equal = true;
        for(i = 0; i < 3 && velocitySource_equal; ++i)
            velocitySource_equal = (velocitySource[i] == obj.velocitySource[i]);

        // Create the return value
        return ((sourceType == obj.sourceType) &&
                Resolution_equal &&
                (UseDataSetStart == obj.UseDataSetStart) &&
                StartPosition_equal &&
                (UseDataSetEnd == obj.UseDataSetEnd) &&
                EndPosition_equal &&
                (integrationDirection == obj.integrationDirection) &&
                (maxSteps == obj.maxSteps) &&
                (terminationType == obj.terminationType) &&
                (terminateBySize == obj.terminateBySize) &&
                (termSize == obj.termSize) &&
                (terminateByDistance == obj.terminateByDistance) &&
                (termDistance == obj.termDistance) &&
                (terminateByTime == obj.terminateByTime) &&
                (termTime == obj.termTime) &&
                (maxStepLength == obj.maxStepLength) &&
                (limitMaximumTimestep == obj.limitMaximumTimestep) &&
                (maxTimeStep == obj.maxTimeStep) &&
                (relTol == obj.relTol) &&
                (absTolSizeType == obj.absTolSizeType) &&
                (absTolAbsolute == obj.absTolAbsolute) &&
                (absTolBBox == obj.absTolBBox) &&
                (fieldType == obj.fieldType) &&
                (fieldConstant == obj.fieldConstant) &&
                velocitySource_equal &&
                (integrationType == obj.integrationType) &&
                (parallelizationAlgorithmType == obj.parallelizationAlgorithmType) &&
                (maxProcessCount == obj.maxProcessCount) &&
                (maxDomainCacheSize == obj.maxDomainCacheSize) &&
                (workGroupSize == obj.workGroupSize) &&
                (pathlines == obj.pathlines) &&
                (pathlinesOverrideStartingTimeFlag == obj.pathlinesOverrideStartingTimeFlag) &&
                (pathlinesOverrideStartingTime == obj.pathlinesOverrideStartingTime) &&
                (pathlinesCMFE == obj.pathlinesCMFE) &&
                (forceNodeCenteredData == obj.forceNodeCenteredData) &&
                (issueTerminationWarnings == obj.issueTerminationWarnings) &&
                (issueStiffnessWarnings == obj.issueStiffnessWarnings) &&
                (issueCriticalPointsWarnings == obj.issueCriticalPointsWarnings) &&
                (criticalPointThreshold == obj.criticalPointThreshold));
    }

    public String GetName() { return "FTLE"; }
    public String GetVersion() { return "1.0"; }

    // Property setting methods
    public void SetSourceType(int sourceType_)
    {
        sourceType = sourceType_;
        Select(0);
    }

    public void SetResolution(int[] Resolution_)
    {
        Resolution[0] = Resolution_[0];
        Resolution[1] = Resolution_[1];
        Resolution[2] = Resolution_[2];
        Select(1);
    }

    public void SetResolution(int e0, int e1, int e2)
    {
        Resolution[0] = e0;
        Resolution[1] = e1;
        Resolution[2] = e2;
        Select(1);
    }

    public void SetUseDataSetStart(int UseDataSetStart_)
    {
        UseDataSetStart = UseDataSetStart_;
        Select(2);
    }

    public void SetStartPosition(double[] StartPosition_)
    {
        StartPosition[0] = StartPosition_[0];
        StartPosition[1] = StartPosition_[1];
        StartPosition[2] = StartPosition_[2];
        Select(3);
    }

    public void SetStartPosition(double e0, double e1, double e2)
    {
        StartPosition[0] = e0;
        StartPosition[1] = e1;
        StartPosition[2] = e2;
        Select(3);
    }

    public void SetUseDataSetEnd(int UseDataSetEnd_)
    {
        UseDataSetEnd = UseDataSetEnd_;
        Select(4);
    }

    public void SetEndPosition(double[] EndPosition_)
    {
        EndPosition[0] = EndPosition_[0];
        EndPosition[1] = EndPosition_[1];
        EndPosition[2] = EndPosition_[2];
        Select(5);
    }

    public void SetEndPosition(double e0, double e1, double e2)
    {
        EndPosition[0] = e0;
        EndPosition[1] = e1;
        EndPosition[2] = e2;
        Select(5);
    }

    public void SetIntegrationDirection(int integrationDirection_)
    {
        integrationDirection = integrationDirection_;
        Select(6);
    }

    public void SetMaxSteps(int maxSteps_)
    {
        maxSteps = maxSteps_;
        Select(7);
    }

    public void SetTerminationType(int terminationType_)
    {
        terminationType = terminationType_;
        Select(8);
    }

    public void SetTerminateBySize(boolean terminateBySize_)
    {
        terminateBySize = terminateBySize_;
        Select(9);
    }

    public void SetTermSize(double termSize_)
    {
        termSize = termSize_;
        Select(10);
    }

    public void SetTerminateByDistance(boolean terminateByDistance_)
    {
        terminateByDistance = terminateByDistance_;
        Select(11);
    }

    public void SetTermDistance(double termDistance_)
    {
        termDistance = termDistance_;
        Select(12);
    }

    public void SetTerminateByTime(boolean terminateByTime_)
    {
        terminateByTime = terminateByTime_;
        Select(13);
    }

    public void SetTermTime(double termTime_)
    {
        termTime = termTime_;
        Select(14);
    }

    public void SetMaxStepLength(double maxStepLength_)
    {
        maxStepLength = maxStepLength_;
        Select(15);
    }

    public void SetLimitMaximumTimestep(boolean limitMaximumTimestep_)
    {
        limitMaximumTimestep = limitMaximumTimestep_;
        Select(16);
    }

    public void SetMaxTimeStep(double maxTimeStep_)
    {
        maxTimeStep = maxTimeStep_;
        Select(17);
    }

    public void SetRelTol(double relTol_)
    {
        relTol = relTol_;
        Select(18);
    }

    public void SetAbsTolSizeType(int absTolSizeType_)
    {
        absTolSizeType = absTolSizeType_;
        Select(19);
    }

    public void SetAbsTolAbsolute(double absTolAbsolute_)
    {
        absTolAbsolute = absTolAbsolute_;
        Select(20);
    }

    public void SetAbsTolBBox(double absTolBBox_)
    {
        absTolBBox = absTolBBox_;
        Select(21);
    }

    public void SetFieldType(int fieldType_)
    {
        fieldType = fieldType_;
        Select(22);
    }

    public void SetFieldConstant(double fieldConstant_)
    {
        fieldConstant = fieldConstant_;
        Select(23);
    }

    public void SetVelocitySource(double[] velocitySource_)
    {
        velocitySource[0] = velocitySource_[0];
        velocitySource[1] = velocitySource_[1];
        velocitySource[2] = velocitySource_[2];
        Select(24);
    }

    public void SetVelocitySource(double e0, double e1, double e2)
    {
        velocitySource[0] = e0;
        velocitySource[1] = e1;
        velocitySource[2] = e2;
        Select(24);
    }

    public void SetIntegrationType(int integrationType_)
    {
        integrationType = integrationType_;
        Select(25);
    }

    public void SetParallelizationAlgorithmType(int parallelizationAlgorithmType_)
    {
        parallelizationAlgorithmType = parallelizationAlgorithmType_;
        Select(26);
    }

    public void SetMaxProcessCount(int maxProcessCount_)
    {
        maxProcessCount = maxProcessCount_;
        Select(27);
    }

    public void SetMaxDomainCacheSize(int maxDomainCacheSize_)
    {
        maxDomainCacheSize = maxDomainCacheSize_;
        Select(28);
    }

    public void SetWorkGroupSize(int workGroupSize_)
    {
        workGroupSize = workGroupSize_;
        Select(29);
    }

    public void SetPathlines(boolean pathlines_)
    {
        pathlines = pathlines_;
        Select(30);
    }

    public void SetPathlinesOverrideStartingTimeFlag(boolean pathlinesOverrideStartingTimeFlag_)
    {
        pathlinesOverrideStartingTimeFlag = pathlinesOverrideStartingTimeFlag_;
        Select(31);
    }

    public void SetPathlinesOverrideStartingTime(double pathlinesOverrideStartingTime_)
    {
        pathlinesOverrideStartingTime = pathlinesOverrideStartingTime_;
        Select(32);
    }

    public void SetPathlinesCMFE(int pathlinesCMFE_)
    {
        pathlinesCMFE = pathlinesCMFE_;
        Select(33);
    }

    public void SetForceNodeCenteredData(boolean forceNodeCenteredData_)
    {
        forceNodeCenteredData = forceNodeCenteredData_;
        Select(34);
    }

    public void SetIssueTerminationWarnings(boolean issueTerminationWarnings_)
    {
        issueTerminationWarnings = issueTerminationWarnings_;
        Select(35);
    }

    public void SetIssueStiffnessWarnings(boolean issueStiffnessWarnings_)
    {
        issueStiffnessWarnings = issueStiffnessWarnings_;
        Select(36);
    }

    public void SetIssueCriticalPointsWarnings(boolean issueCriticalPointsWarnings_)
    {
        issueCriticalPointsWarnings = issueCriticalPointsWarnings_;
        Select(37);
    }

    public void SetCriticalPointThreshold(double criticalPointThreshold_)
    {
        criticalPointThreshold = criticalPointThreshold_;
        Select(38);
    }

    // Property getting methods
    public int      GetSourceType() { return sourceType; }
    public int[]    GetResolution() { return Resolution; }
    public int      GetUseDataSetStart() { return UseDataSetStart; }
    public double[] GetStartPosition() { return StartPosition; }
    public int      GetUseDataSetEnd() { return UseDataSetEnd; }
    public double[] GetEndPosition() { return EndPosition; }
    public int      GetIntegrationDirection() { return integrationDirection; }
    public int      GetMaxSteps() { return maxSteps; }
    public int      GetTerminationType() { return terminationType; }
    public boolean  GetTerminateBySize() { return terminateBySize; }
    public double   GetTermSize() { return termSize; }
    public boolean  GetTerminateByDistance() { return terminateByDistance; }
    public double   GetTermDistance() { return termDistance; }
    public boolean  GetTerminateByTime() { return terminateByTime; }
    public double   GetTermTime() { return termTime; }
    public double   GetMaxStepLength() { return maxStepLength; }
    public boolean  GetLimitMaximumTimestep() { return limitMaximumTimestep; }
    public double   GetMaxTimeStep() { return maxTimeStep; }
    public double   GetRelTol() { return relTol; }
    public int      GetAbsTolSizeType() { return absTolSizeType; }
    public double   GetAbsTolAbsolute() { return absTolAbsolute; }
    public double   GetAbsTolBBox() { return absTolBBox; }
    public int      GetFieldType() { return fieldType; }
    public double   GetFieldConstant() { return fieldConstant; }
    public double[] GetVelocitySource() { return velocitySource; }
    public int      GetIntegrationType() { return integrationType; }
    public int      GetParallelizationAlgorithmType() { return parallelizationAlgorithmType; }
    public int      GetMaxProcessCount() { return maxProcessCount; }
    public int      GetMaxDomainCacheSize() { return maxDomainCacheSize; }
    public int      GetWorkGroupSize() { return workGroupSize; }
    public boolean  GetPathlines() { return pathlines; }
    public boolean  GetPathlinesOverrideStartingTimeFlag() { return pathlinesOverrideStartingTimeFlag; }
    public double   GetPathlinesOverrideStartingTime() { return pathlinesOverrideStartingTime; }
    public int      GetPathlinesCMFE() { return pathlinesCMFE; }
    public boolean  GetForceNodeCenteredData() { return forceNodeCenteredData; }
    public boolean  GetIssueTerminationWarnings() { return issueTerminationWarnings; }
    public boolean  GetIssueStiffnessWarnings() { return issueStiffnessWarnings; }
    public boolean  GetIssueCriticalPointsWarnings() { return issueCriticalPointsWarnings; }
    public double   GetCriticalPointThreshold() { return criticalPointThreshold; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteInt(sourceType);
        if(WriteSelect(1, buf))
            buf.WriteIntArray(Resolution);
        if(WriteSelect(2, buf))
            buf.WriteInt(UseDataSetStart);
        if(WriteSelect(3, buf))
            buf.WriteDoubleArray(StartPosition);
        if(WriteSelect(4, buf))
            buf.WriteInt(UseDataSetEnd);
        if(WriteSelect(5, buf))
            buf.WriteDoubleArray(EndPosition);
        if(WriteSelect(6, buf))
            buf.WriteInt(integrationDirection);
        if(WriteSelect(7, buf))
            buf.WriteInt(maxSteps);
        if(WriteSelect(8, buf))
            buf.WriteInt(terminationType);
        if(WriteSelect(9, buf))
            buf.WriteBool(terminateBySize);
        if(WriteSelect(10, buf))
            buf.WriteDouble(termSize);
        if(WriteSelect(11, buf))
            buf.WriteBool(terminateByDistance);
        if(WriteSelect(12, buf))
            buf.WriteDouble(termDistance);
        if(WriteSelect(13, buf))
            buf.WriteBool(terminateByTime);
        if(WriteSelect(14, buf))
            buf.WriteDouble(termTime);
        if(WriteSelect(15, buf))
            buf.WriteDouble(maxStepLength);
        if(WriteSelect(16, buf))
            buf.WriteBool(limitMaximumTimestep);
        if(WriteSelect(17, buf))
            buf.WriteDouble(maxTimeStep);
        if(WriteSelect(18, buf))
            buf.WriteDouble(relTol);
        if(WriteSelect(19, buf))
            buf.WriteInt(absTolSizeType);
        if(WriteSelect(20, buf))
            buf.WriteDouble(absTolAbsolute);
        if(WriteSelect(21, buf))
            buf.WriteDouble(absTolBBox);
        if(WriteSelect(22, buf))
            buf.WriteInt(fieldType);
        if(WriteSelect(23, buf))
            buf.WriteDouble(fieldConstant);
        if(WriteSelect(24, buf))
            buf.WriteDoubleArray(velocitySource);
        if(WriteSelect(25, buf))
            buf.WriteInt(integrationType);
        if(WriteSelect(26, buf))
            buf.WriteInt(parallelizationAlgorithmType);
        if(WriteSelect(27, buf))
            buf.WriteInt(maxProcessCount);
        if(WriteSelect(28, buf))
            buf.WriteInt(maxDomainCacheSize);
        if(WriteSelect(29, buf))
            buf.WriteInt(workGroupSize);
        if(WriteSelect(30, buf))
            buf.WriteBool(pathlines);
        if(WriteSelect(31, buf))
            buf.WriteBool(pathlinesOverrideStartingTimeFlag);
        if(WriteSelect(32, buf))
            buf.WriteDouble(pathlinesOverrideStartingTime);
        if(WriteSelect(33, buf))
            buf.WriteInt(pathlinesCMFE);
        if(WriteSelect(34, buf))
            buf.WriteBool(forceNodeCenteredData);
        if(WriteSelect(35, buf))
            buf.WriteBool(issueTerminationWarnings);
        if(WriteSelect(36, buf))
            buf.WriteBool(issueStiffnessWarnings);
        if(WriteSelect(37, buf))
            buf.WriteBool(issueCriticalPointsWarnings);
        if(WriteSelect(38, buf))
            buf.WriteDouble(criticalPointThreshold);
    }

    public void ReadAtts(int index, CommunicationBuffer buf)
    {
        switch(index)
        {
        case 0:
            SetSourceType(buf.ReadInt());
            break;
        case 1:
            SetResolution(buf.ReadIntArray());
            break;
        case 2:
            SetUseDataSetStart(buf.ReadInt());
            break;
        case 3:
            SetStartPosition(buf.ReadDoubleArray());
            break;
        case 4:
            SetUseDataSetEnd(buf.ReadInt());
            break;
        case 5:
            SetEndPosition(buf.ReadDoubleArray());
            break;
        case 6:
            SetIntegrationDirection(buf.ReadInt());
            break;
        case 7:
            SetMaxSteps(buf.ReadInt());
            break;
        case 8:
            SetTerminationType(buf.ReadInt());
            break;
        case 9:
            SetTerminateBySize(buf.ReadBool());
            break;
        case 10:
            SetTermSize(buf.ReadDouble());
            break;
        case 11:
            SetTerminateByDistance(buf.ReadBool());
            break;
        case 12:
            SetTermDistance(buf.ReadDouble());
            break;
        case 13:
            SetTerminateByTime(buf.ReadBool());
            break;
        case 14:
            SetTermTime(buf.ReadDouble());
            break;
        case 15:
            SetMaxStepLength(buf.ReadDouble());
            break;
        case 16:
            SetLimitMaximumTimestep(buf.ReadBool());
            break;
        case 17:
            SetMaxTimeStep(buf.ReadDouble());
            break;
        case 18:
            SetRelTol(buf.ReadDouble());
            break;
        case 19:
            SetAbsTolSizeType(buf.ReadInt());
            break;
        case 20:
            SetAbsTolAbsolute(buf.ReadDouble());
            break;
        case 21:
            SetAbsTolBBox(buf.ReadDouble());
            break;
        case 22:
            SetFieldType(buf.ReadInt());
            break;
        case 23:
            SetFieldConstant(buf.ReadDouble());
            break;
        case 24:
            SetVelocitySource(buf.ReadDoubleArray());
            break;
        case 25:
            SetIntegrationType(buf.ReadInt());
            break;
        case 26:
            SetParallelizationAlgorithmType(buf.ReadInt());
            break;
        case 27:
            SetMaxProcessCount(buf.ReadInt());
            break;
        case 28:
            SetMaxDomainCacheSize(buf.ReadInt());
            break;
        case 29:
            SetWorkGroupSize(buf.ReadInt());
            break;
        case 30:
            SetPathlines(buf.ReadBool());
            break;
        case 31:
            SetPathlinesOverrideStartingTimeFlag(buf.ReadBool());
            break;
        case 32:
            SetPathlinesOverrideStartingTime(buf.ReadDouble());
            break;
        case 33:
            SetPathlinesCMFE(buf.ReadInt());
            break;
        case 34:
            SetForceNodeCenteredData(buf.ReadBool());
            break;
        case 35:
            SetIssueTerminationWarnings(buf.ReadBool());
            break;
        case 36:
            SetIssueStiffnessWarnings(buf.ReadBool());
            break;
        case 37:
            SetIssueCriticalPointsWarnings(buf.ReadBool());
            break;
        case 38:
            SetCriticalPointThreshold(buf.ReadDouble());
            break;
        }
    }

    public String toString(String indent)
    {
        String str = new String();
        str = str + indent + "sourceType = ";
        if(sourceType == SOURCETYPE_NATIVERESOLUTIONOFMESH)
            str = str + "SOURCETYPE_NATIVERESOLUTIONOFMESH";
        if(sourceType == SOURCETYPE_REGULARGRID)
            str = str + "SOURCETYPE_REGULARGRID";
        str = str + "\n";
        str = str + intArrayToString("Resolution", Resolution, indent) + "\n";
        str = str + indent + "UseDataSetStart = ";
        if(UseDataSetStart == EXTENTS_FULL)
            str = str + "EXTENTS_FULL";
        if(UseDataSetStart == EXTENTS_SUBSET)
            str = str + "EXTENTS_SUBSET";
        str = str + "\n";
        str = str + doubleArrayToString("StartPosition", StartPosition, indent) + "\n";
        str = str + indent + "UseDataSetEnd = ";
        if(UseDataSetEnd == EXTENTS_FULL)
            str = str + "EXTENTS_FULL";
        if(UseDataSetEnd == EXTENTS_SUBSET)
            str = str + "EXTENTS_SUBSET";
        str = str + "\n";
        str = str + doubleArrayToString("EndPosition", EndPosition, indent) + "\n";
        str = str + indent + "integrationDirection = ";
        if(integrationDirection == INTEGRATIONDIRECTION_FORWARD)
            str = str + "INTEGRATIONDIRECTION_FORWARD";
        if(integrationDirection == INTEGRATIONDIRECTION_BACKWARD)
            str = str + "INTEGRATIONDIRECTION_BACKWARD";
        if(integrationDirection == INTEGRATIONDIRECTION_BOTH)
            str = str + "INTEGRATIONDIRECTION_BOTH";
        str = str + "\n";
        str = str + intToString("maxSteps", maxSteps, indent) + "\n";
        str = str + indent + "terminationType = ";
        if(terminationType == TERMINATIONTYPE_TIME)
            str = str + "TERMINATIONTYPE_TIME";
        if(terminationType == TERMINATIONTYPE_DISTANCE)
            str = str + "TERMINATIONTYPE_DISTANCE";
        if(terminationType == TERMINATIONTYPE_SIZE)
            str = str + "TERMINATIONTYPE_SIZE";
        str = str + "\n";
        str = str + boolToString("terminateBySize", terminateBySize, indent) + "\n";
        str = str + doubleToString("termSize", termSize, indent) + "\n";
        str = str + boolToString("terminateByDistance", terminateByDistance, indent) + "\n";
        str = str + doubleToString("termDistance", termDistance, indent) + "\n";
        str = str + boolToString("terminateByTime", terminateByTime, indent) + "\n";
        str = str + doubleToString("termTime", termTime, indent) + "\n";
        str = str + doubleToString("maxStepLength", maxStepLength, indent) + "\n";
        str = str + boolToString("limitMaximumTimestep", limitMaximumTimestep, indent) + "\n";
        str = str + doubleToString("maxTimeStep", maxTimeStep, indent) + "\n";
        str = str + doubleToString("relTol", relTol, indent) + "\n";
        str = str + indent + "absTolSizeType = ";
        if(absTolSizeType == SIZETYPE_ABSOLUTE)
            str = str + "SIZETYPE_ABSOLUTE";
        if(absTolSizeType == SIZETYPE_FRACTIONOFBBOX)
            str = str + "SIZETYPE_FRACTIONOFBBOX";
        str = str + "\n";
        str = str + doubleToString("absTolAbsolute", absTolAbsolute, indent) + "\n";
        str = str + doubleToString("absTolBBox", absTolBBox, indent) + "\n";
        str = str + indent + "fieldType = ";
        if(fieldType == FIELDTYPE_DEFAULT)
            str = str + "FIELDTYPE_DEFAULT";
        if(fieldType == FIELDTYPE_FLASHFIELD)
            str = str + "FIELDTYPE_FLASHFIELD";
        if(fieldType == FIELDTYPE_M3DC12DFIELD)
            str = str + "FIELDTYPE_M3DC12DFIELD";
        if(fieldType == FIELDTYPE_M3DC13DFIELD)
            str = str + "FIELDTYPE_M3DC13DFIELD";
        if(fieldType == FIELDTYPE_NEK5000FIELD)
            str = str + "FIELDTYPE_NEK5000FIELD";
        if(fieldType == FIELDTYPE_NIMRODFIELD)
            str = str + "FIELDTYPE_NIMRODFIELD";
        str = str + "\n";
        str = str + doubleToString("fieldConstant", fieldConstant, indent) + "\n";
        str = str + doubleArrayToString("velocitySource", velocitySource, indent) + "\n";
        str = str + indent + "integrationType = ";
        if(integrationType == INTEGRATIONTYPE_EULER)
            str = str + "INTEGRATIONTYPE_EULER";
        if(integrationType == INTEGRATIONTYPE_LEAPFROG)
            str = str + "INTEGRATIONTYPE_LEAPFROG";
        if(integrationType == INTEGRATIONTYPE_DORMANDPRINCE)
            str = str + "INTEGRATIONTYPE_DORMANDPRINCE";
        if(integrationType == INTEGRATIONTYPE_ADAMSBASHFORTH)
            str = str + "INTEGRATIONTYPE_ADAMSBASHFORTH";
        if(integrationType == INTEGRATIONTYPE_RK4)
            str = str + "INTEGRATIONTYPE_RK4";
        if(integrationType == INTEGRATIONTYPE_M3DC12DINTEGRATOR)
            str = str + "INTEGRATIONTYPE_M3DC12DINTEGRATOR";
        str = str + "\n";
        str = str + indent + "parallelizationAlgorithmType = ";
        if(parallelizationAlgorithmType == PARALLELIZATIONALGORITHMTYPE_LOADONDEMAND)
            str = str + "PARALLELIZATIONALGORITHMTYPE_LOADONDEMAND";
        if(parallelizationAlgorithmType == PARALLELIZATIONALGORITHMTYPE_PARALLELSTATICDOMAINS)
            str = str + "PARALLELIZATIONALGORITHMTYPE_PARALLELSTATICDOMAINS";
        if(parallelizationAlgorithmType == PARALLELIZATIONALGORITHMTYPE_MASTERSLAVE)
            str = str + "PARALLELIZATIONALGORITHMTYPE_MASTERSLAVE";
        if(parallelizationAlgorithmType == PARALLELIZATIONALGORITHMTYPE_VISITSELECTS)
            str = str + "PARALLELIZATIONALGORITHMTYPE_VISITSELECTS";
        str = str + "\n";
        str = str + intToString("maxProcessCount", maxProcessCount, indent) + "\n";
        str = str + intToString("maxDomainCacheSize", maxDomainCacheSize, indent) + "\n";
        str = str + intToString("workGroupSize", workGroupSize, indent) + "\n";
        str = str + boolToString("pathlines", pathlines, indent) + "\n";
        str = str + boolToString("pathlinesOverrideStartingTimeFlag", pathlinesOverrideStartingTimeFlag, indent) + "\n";
        str = str + doubleToString("pathlinesOverrideStartingTime", pathlinesOverrideStartingTime, indent) + "\n";
        str = str + indent + "pathlinesCMFE = ";
        if(pathlinesCMFE == PATHLINESCMFE_CONN_CMFE)
            str = str + "PATHLINESCMFE_CONN_CMFE";
        if(pathlinesCMFE == PATHLINESCMFE_POS_CMFE)
            str = str + "PATHLINESCMFE_POS_CMFE";
        str = str + "\n";
        str = str + boolToString("forceNodeCenteredData", forceNodeCenteredData, indent) + "\n";
        str = str + boolToString("issueTerminationWarnings", issueTerminationWarnings, indent) + "\n";
        str = str + boolToString("issueStiffnessWarnings", issueStiffnessWarnings, indent) + "\n";
        str = str + boolToString("issueCriticalPointsWarnings", issueCriticalPointsWarnings, indent) + "\n";
        str = str + doubleToString("criticalPointThreshold", criticalPointThreshold, indent) + "\n";
        return str;
    }


    // Attributes
    private int      sourceType;
    private int[]    Resolution;
    private int      UseDataSetStart;
    private double[] StartPosition;
    private int      UseDataSetEnd;
    private double[] EndPosition;
    private int      integrationDirection;
    private int      maxSteps;
    private int      terminationType;
    private boolean  terminateBySize;
    private double   termSize;
    private boolean  terminateByDistance;
    private double   termDistance;
    private boolean  terminateByTime;
    private double   termTime;
    private double   maxStepLength;
    private boolean  limitMaximumTimestep;
    private double   maxTimeStep;
    private double   relTol;
    private int      absTolSizeType;
    private double   absTolAbsolute;
    private double   absTolBBox;
    private int      fieldType;
    private double   fieldConstant;
    private double[] velocitySource;
    private int      integrationType;
    private int      parallelizationAlgorithmType;
    private int      maxProcessCount;
    private int      maxDomainCacheSize;
    private int      workGroupSize;
    private boolean  pathlines;
    private boolean  pathlinesOverrideStartingTimeFlag;
    private double   pathlinesOverrideStartingTime;
    private int      pathlinesCMFE;
    private boolean  forceNodeCenteredData;
    private boolean  issueTerminationWarnings;
    private boolean  issueStiffnessWarnings;
    private boolean  issueCriticalPointsWarnings;
    private double   criticalPointThreshold;
}

