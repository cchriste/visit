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

package llnl.visit;


// ****************************************************************************
// Class: InternalResampleAttributes
//
// Purpose:
//    This class contains attributes to specify a resampling.
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

public class InternalResampleAttributes extends AttributeSubject
{
    private static int InternalResampleAttributes_numAdditionalAtts = 18;

    public InternalResampleAttributes()
    {
        super(InternalResampleAttributes_numAdditionalAtts);

        useTargetVal = false;
        targetVal = 100000;
        width = 30;
        height = 30;
        depth = 30;
        prefersPowersOfTwo = false;
        defaultVal = -1e+38f;
        useBounds = false;
        minX = 0;
        minY = 0;
        minZ = 0;
        maxX = 1;
        maxY = 1;
        maxZ = 1;
        useArbitrator = false;
        arbitratorLessThan = false;
        arbitratorVarName = new String("default");
        distributedResample = false;
    }

    public InternalResampleAttributes(int nMoreFields)
    {
        super(InternalResampleAttributes_numAdditionalAtts + nMoreFields);

        useTargetVal = false;
        targetVal = 100000;
        width = 30;
        height = 30;
        depth = 30;
        prefersPowersOfTwo = false;
        defaultVal = -1e+38f;
        useBounds = false;
        minX = 0;
        minY = 0;
        minZ = 0;
        maxX = 1;
        maxY = 1;
        maxZ = 1;
        useArbitrator = false;
        arbitratorLessThan = false;
        arbitratorVarName = new String("default");
        distributedResample = false;
    }

    public InternalResampleAttributes(InternalResampleAttributes obj)
    {
        super(InternalResampleAttributes_numAdditionalAtts);

        useTargetVal = obj.useTargetVal;
        targetVal = obj.targetVal;
        width = obj.width;
        height = obj.height;
        depth = obj.depth;
        prefersPowersOfTwo = obj.prefersPowersOfTwo;
        defaultVal = obj.defaultVal;
        useBounds = obj.useBounds;
        minX = obj.minX;
        minY = obj.minY;
        minZ = obj.minZ;
        maxX = obj.maxX;
        maxY = obj.maxY;
        maxZ = obj.maxZ;
        useArbitrator = obj.useArbitrator;
        arbitratorLessThan = obj.arbitratorLessThan;
        arbitratorVarName = new String(obj.arbitratorVarName);
        distributedResample = obj.distributedResample;

        SelectAll();
    }

    public int Offset()
    {
        return super.Offset() + super.GetNumAdditionalAttributes();
    }

    public int GetNumAdditionalAttributes()
    {
        return InternalResampleAttributes_numAdditionalAtts;
    }

    public boolean equals(InternalResampleAttributes obj)
    {
        // Create the return value
        return ((useTargetVal == obj.useTargetVal) &&
                (targetVal == obj.targetVal) &&
                (width == obj.width) &&
                (height == obj.height) &&
                (depth == obj.depth) &&
                (prefersPowersOfTwo == obj.prefersPowersOfTwo) &&
                (defaultVal == obj.defaultVal) &&
                (useBounds == obj.useBounds) &&
                (minX == obj.minX) &&
                (minY == obj.minY) &&
                (minZ == obj.minZ) &&
                (maxX == obj.maxX) &&
                (maxY == obj.maxY) &&
                (maxZ == obj.maxZ) &&
                (useArbitrator == obj.useArbitrator) &&
                (arbitratorLessThan == obj.arbitratorLessThan) &&
                (arbitratorVarName.equals(obj.arbitratorVarName)) &&
                (distributedResample == obj.distributedResample));
    }

    // Property setting methods
    public void SetUseTargetVal(boolean useTargetVal_)
    {
        useTargetVal = useTargetVal_;
        Select(0);
    }

    public void SetTargetVal(int targetVal_)
    {
        targetVal = targetVal_;
        Select(1);
    }

    public void SetWidth(int width_)
    {
        width = width_;
        Select(2);
    }

    public void SetHeight(int height_)
    {
        height = height_;
        Select(3);
    }

    public void SetDepth(int depth_)
    {
        depth = depth_;
        Select(4);
    }

    public void SetPrefersPowersOfTwo(boolean prefersPowersOfTwo_)
    {
        prefersPowersOfTwo = prefersPowersOfTwo_;
        Select(5);
    }

    public void SetDefaultVal(float defaultVal_)
    {
        defaultVal = defaultVal_;
        Select(6);
    }

    public void SetUseBounds(boolean useBounds_)
    {
        useBounds = useBounds_;
        Select(7);
    }

    public void SetMinX(double minX_)
    {
        minX = minX_;
        Select(8);
    }

    public void SetMinY(double minY_)
    {
        minY = minY_;
        Select(9);
    }

    public void SetMinZ(double minZ_)
    {
        minZ = minZ_;
        Select(10);
    }

    public void SetMaxX(double maxX_)
    {
        maxX = maxX_;
        Select(11);
    }

    public void SetMaxY(double maxY_)
    {
        maxY = maxY_;
        Select(12);
    }

    public void SetMaxZ(double maxZ_)
    {
        maxZ = maxZ_;
        Select(13);
    }

    public void SetUseArbitrator(boolean useArbitrator_)
    {
        useArbitrator = useArbitrator_;
        Select(14);
    }

    public void SetArbitratorLessThan(boolean arbitratorLessThan_)
    {
        arbitratorLessThan = arbitratorLessThan_;
        Select(15);
    }

    public void SetArbitratorVarName(String arbitratorVarName_)
    {
        arbitratorVarName = arbitratorVarName_;
        Select(16);
    }

    public void SetDistributedResample(boolean distributedResample_)
    {
        distributedResample = distributedResample_;
        Select(17);
    }

    // Property getting methods
    public boolean GetUseTargetVal() { return useTargetVal; }
    public int     GetTargetVal() { return targetVal; }
    public int     GetWidth() { return width; }
    public int     GetHeight() { return height; }
    public int     GetDepth() { return depth; }
    public boolean GetPrefersPowersOfTwo() { return prefersPowersOfTwo; }
    public float   GetDefaultVal() { return defaultVal; }
    public boolean GetUseBounds() { return useBounds; }
    public double  GetMinX() { return minX; }
    public double  GetMinY() { return minY; }
    public double  GetMinZ() { return minZ; }
    public double  GetMaxX() { return maxX; }
    public double  GetMaxY() { return maxY; }
    public double  GetMaxZ() { return maxZ; }
    public boolean GetUseArbitrator() { return useArbitrator; }
    public boolean GetArbitratorLessThan() { return arbitratorLessThan; }
    public String  GetArbitratorVarName() { return arbitratorVarName; }
    public boolean GetDistributedResample() { return distributedResample; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteBool(useTargetVal);
        if(WriteSelect(1, buf))
            buf.WriteInt(targetVal);
        if(WriteSelect(2, buf))
            buf.WriteInt(width);
        if(WriteSelect(3, buf))
            buf.WriteInt(height);
        if(WriteSelect(4, buf))
            buf.WriteInt(depth);
        if(WriteSelect(5, buf))
            buf.WriteBool(prefersPowersOfTwo);
        if(WriteSelect(6, buf))
            buf.WriteFloat(defaultVal);
        if(WriteSelect(7, buf))
            buf.WriteBool(useBounds);
        if(WriteSelect(8, buf))
            buf.WriteDouble(minX);
        if(WriteSelect(9, buf))
            buf.WriteDouble(minY);
        if(WriteSelect(10, buf))
            buf.WriteDouble(minZ);
        if(WriteSelect(11, buf))
            buf.WriteDouble(maxX);
        if(WriteSelect(12, buf))
            buf.WriteDouble(maxY);
        if(WriteSelect(13, buf))
            buf.WriteDouble(maxZ);
        if(WriteSelect(14, buf))
            buf.WriteBool(useArbitrator);
        if(WriteSelect(15, buf))
            buf.WriteBool(arbitratorLessThan);
        if(WriteSelect(16, buf))
            buf.WriteString(arbitratorVarName);
        if(WriteSelect(17, buf))
            buf.WriteBool(distributedResample);
    }

    public void ReadAtts(int index, CommunicationBuffer buf)
    {
        switch(index)
        {
        case 0:
            SetUseTargetVal(buf.ReadBool());
            break;
        case 1:
            SetTargetVal(buf.ReadInt());
            break;
        case 2:
            SetWidth(buf.ReadInt());
            break;
        case 3:
            SetHeight(buf.ReadInt());
            break;
        case 4:
            SetDepth(buf.ReadInt());
            break;
        case 5:
            SetPrefersPowersOfTwo(buf.ReadBool());
            break;
        case 6:
            SetDefaultVal(buf.ReadFloat());
            break;
        case 7:
            SetUseBounds(buf.ReadBool());
            break;
        case 8:
            SetMinX(buf.ReadDouble());
            break;
        case 9:
            SetMinY(buf.ReadDouble());
            break;
        case 10:
            SetMinZ(buf.ReadDouble());
            break;
        case 11:
            SetMaxX(buf.ReadDouble());
            break;
        case 12:
            SetMaxY(buf.ReadDouble());
            break;
        case 13:
            SetMaxZ(buf.ReadDouble());
            break;
        case 14:
            SetUseArbitrator(buf.ReadBool());
            break;
        case 15:
            SetArbitratorLessThan(buf.ReadBool());
            break;
        case 16:
            SetArbitratorVarName(buf.ReadString());
            break;
        case 17:
            SetDistributedResample(buf.ReadBool());
            break;
        }
    }

    public String toString(String indent)
    {
        String str = new String();
        str = str + boolToString("useTargetVal", useTargetVal, indent) + "\n";
        str = str + intToString("targetVal", targetVal, indent) + "\n";
        str = str + intToString("width", width, indent) + "\n";
        str = str + intToString("height", height, indent) + "\n";
        str = str + intToString("depth", depth, indent) + "\n";
        str = str + boolToString("prefersPowersOfTwo", prefersPowersOfTwo, indent) + "\n";
        str = str + floatToString("defaultVal", defaultVal, indent) + "\n";
        str = str + boolToString("useBounds", useBounds, indent) + "\n";
        str = str + doubleToString("minX", minX, indent) + "\n";
        str = str + doubleToString("minY", minY, indent) + "\n";
        str = str + doubleToString("minZ", minZ, indent) + "\n";
        str = str + doubleToString("maxX", maxX, indent) + "\n";
        str = str + doubleToString("maxY", maxY, indent) + "\n";
        str = str + doubleToString("maxZ", maxZ, indent) + "\n";
        str = str + boolToString("useArbitrator", useArbitrator, indent) + "\n";
        str = str + boolToString("arbitratorLessThan", arbitratorLessThan, indent) + "\n";
        str = str + stringToString("arbitratorVarName", arbitratorVarName, indent) + "\n";
        str = str + boolToString("distributedResample", distributedResample, indent) + "\n";
        return str;
    }


    // Attributes
    private boolean useTargetVal;
    private int     targetVal;
    private int     width;
    private int     height;
    private int     depth;
    private boolean prefersPowersOfTwo;
    private float   defaultVal;
    private boolean useBounds;
    private double  minX;
    private double  minY;
    private double  minZ;
    private double  maxX;
    private double  maxY;
    private double  maxZ;
    private boolean useArbitrator;
    private boolean arbitratorLessThan;
    private String  arbitratorVarName;
    private boolean distributedResample;
}

