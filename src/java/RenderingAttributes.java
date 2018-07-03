// ***************************************************************************
//
// Copyright (c) 2000 - 2018, Lawrence Livermore National Security, LLC
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
// Class: RenderingAttributes
//
// Purpose:
//    This class contains special rendering attributes like antialiasing and stero settings.
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

public class RenderingAttributes extends AttributeSubject
{
    private static int RenderingAttributes_numAdditionalAtts = 35;

    // Enum values
    public final static int GEOMETRYREPRESENTATION_SURFACES = 0;
    public final static int GEOMETRYREPRESENTATION_WIREFRAME = 1;
    public final static int GEOMETRYREPRESENTATION_POINTS = 2;

    public final static int STEREOTYPES_REDBLUE = 0;
    public final static int STEREOTYPES_INTERLACED = 1;
    public final static int STEREOTYPES_CRYSTALEYES = 2;
    public final static int STEREOTYPES_REDGREEN = 3;

    public final static int TRISTATEMODE_NEVER = 0;
    public final static int TRISTATEMODE_ALWAYS = 1;
    public final static int TRISTATEMODE_AUTO = 2;

    // Constants
public final static int DEFAULT_SCALABLE_AUTO_THRESHOLD = 2000000;

public final static int DEFAULT_SCALABLE_ACTIVATION_MODE = TRISTATEMODE_AUTO;

public final static int DEFAULT_COMPACT_DOMAINS_ACTIVATION_MODE = TRISTATEMODE_AUTO;

public final static int DEFAULT_COMPACT_DOMAINS_AUTO_THRESHOLD = 256;


    public RenderingAttributes()
    {
        super(RenderingAttributes_numAdditionalAtts);

        antialiasing = false;
        orderComposite = true;
        depthCompositeThreads = 2;
        depthCompositeBlocking = 65536;
        alphaCompositeThreads = 2;
        alphaCompositeBlocking = 65536;
        depthPeeling = false;
        occlusionRatio = 0;
        numberOfPeels = 16;
        multiresolutionMode = false;
        multiresolutionCellSize = 0.002f;
        geometryRepresentation = GEOMETRYREPRESENTATION_SURFACES;
        stereoRendering = false;
        stereoType = STEREOTYPES_CRYSTALEYES;
        notifyForEachRender = false;
        scalableActivationMode = TRISTATEMODE_AUTO;
        scalableAutoThreshold = 2000000;
        specularFlag = false;
        specularCoeff = 0.6f;
        specularPower = 10f;
        specularColor = new ColorAttribute(255, 255, 255);
        doShadowing = false;
        shadowStrength = 0.5;
        doDepthCueing = false;
        depthCueingAutomatic = true;
        startCuePoint = new double[3];
        startCuePoint[0] = -10;
        startCuePoint[1] = 0;
        startCuePoint[2] = 0;
        endCuePoint = new double[3];
        endCuePoint[0] = 10;
        endCuePoint[1] = 0;
        endCuePoint[2] = 0;
        compressionActivationMode = TRISTATEMODE_NEVER;
        colorTexturingFlag = true;
        compactDomainsActivationMode = TRISTATEMODE_NEVER;
        compactDomainsAutoThreshold = 256;
        osprayRendering = false;
        ospraySPP = 1;
        osprayAO = 0;
        osprayShadows = false;
    }

    public RenderingAttributes(int nMoreFields)
    {
        super(RenderingAttributes_numAdditionalAtts + nMoreFields);

        antialiasing = false;
        orderComposite = true;
        depthCompositeThreads = 2;
        depthCompositeBlocking = 65536;
        alphaCompositeThreads = 2;
        alphaCompositeBlocking = 65536;
        depthPeeling = false;
        occlusionRatio = 0;
        numberOfPeels = 16;
        multiresolutionMode = false;
        multiresolutionCellSize = 0.002f;
        geometryRepresentation = GEOMETRYREPRESENTATION_SURFACES;
        stereoRendering = false;
        stereoType = STEREOTYPES_CRYSTALEYES;
        notifyForEachRender = false;
        scalableActivationMode = TRISTATEMODE_AUTO;
        scalableAutoThreshold = 2000000;
        specularFlag = false;
        specularCoeff = 0.6f;
        specularPower = 10f;
        specularColor = new ColorAttribute(255, 255, 255);
        doShadowing = false;
        shadowStrength = 0.5;
        doDepthCueing = false;
        depthCueingAutomatic = true;
        startCuePoint = new double[3];
        startCuePoint[0] = -10;
        startCuePoint[1] = 0;
        startCuePoint[2] = 0;
        endCuePoint = new double[3];
        endCuePoint[0] = 10;
        endCuePoint[1] = 0;
        endCuePoint[2] = 0;
        compressionActivationMode = TRISTATEMODE_NEVER;
        colorTexturingFlag = true;
        compactDomainsActivationMode = TRISTATEMODE_NEVER;
        compactDomainsAutoThreshold = 256;
        osprayRendering = false;
        ospraySPP = 1;
        osprayAO = 0;
        osprayShadows = false;
    }

    public RenderingAttributes(RenderingAttributes obj)
    {
        super(obj);

        int i;

        antialiasing = obj.antialiasing;
        orderComposite = obj.orderComposite;
        depthCompositeThreads = obj.depthCompositeThreads;
        depthCompositeBlocking = obj.depthCompositeBlocking;
        alphaCompositeThreads = obj.alphaCompositeThreads;
        alphaCompositeBlocking = obj.alphaCompositeBlocking;
        depthPeeling = obj.depthPeeling;
        occlusionRatio = obj.occlusionRatio;
        numberOfPeels = obj.numberOfPeels;
        multiresolutionMode = obj.multiresolutionMode;
        multiresolutionCellSize = obj.multiresolutionCellSize;
        geometryRepresentation = obj.geometryRepresentation;
        stereoRendering = obj.stereoRendering;
        stereoType = obj.stereoType;
        notifyForEachRender = obj.notifyForEachRender;
        scalableActivationMode = obj.scalableActivationMode;
        scalableAutoThreshold = obj.scalableAutoThreshold;
        specularFlag = obj.specularFlag;
        specularCoeff = obj.specularCoeff;
        specularPower = obj.specularPower;
        specularColor = new ColorAttribute(obj.specularColor);
        doShadowing = obj.doShadowing;
        shadowStrength = obj.shadowStrength;
        doDepthCueing = obj.doDepthCueing;
        depthCueingAutomatic = obj.depthCueingAutomatic;
        startCuePoint = new double[3];
        startCuePoint[0] = obj.startCuePoint[0];
        startCuePoint[1] = obj.startCuePoint[1];
        startCuePoint[2] = obj.startCuePoint[2];

        endCuePoint = new double[3];
        endCuePoint[0] = obj.endCuePoint[0];
        endCuePoint[1] = obj.endCuePoint[1];
        endCuePoint[2] = obj.endCuePoint[2];

        compressionActivationMode = obj.compressionActivationMode;
        colorTexturingFlag = obj.colorTexturingFlag;
        compactDomainsActivationMode = obj.compactDomainsActivationMode;
        compactDomainsAutoThreshold = obj.compactDomainsAutoThreshold;
        osprayRendering = obj.osprayRendering;
        ospraySPP = obj.ospraySPP;
        osprayAO = obj.osprayAO;
        osprayShadows = obj.osprayShadows;

        SelectAll();
    }

    public int Offset()
    {
        return super.Offset() + super.GetNumAdditionalAttributes();
    }

    public int GetNumAdditionalAttributes()
    {
        return RenderingAttributes_numAdditionalAtts;
    }

    public boolean equals(RenderingAttributes obj)
    {
        int i;

        // Compare the startCuePoint arrays.
        boolean startCuePoint_equal = true;
        for(i = 0; i < 3 && startCuePoint_equal; ++i)
            startCuePoint_equal = (startCuePoint[i] == obj.startCuePoint[i]);

        // Compare the endCuePoint arrays.
        boolean endCuePoint_equal = true;
        for(i = 0; i < 3 && endCuePoint_equal; ++i)
            endCuePoint_equal = (endCuePoint[i] == obj.endCuePoint[i]);

        // Create the return value
        return ((antialiasing == obj.antialiasing) &&
                (orderComposite == obj.orderComposite) &&
                (depthCompositeThreads == obj.depthCompositeThreads) &&
                (depthCompositeBlocking == obj.depthCompositeBlocking) &&
                (alphaCompositeThreads == obj.alphaCompositeThreads) &&
                (alphaCompositeBlocking == obj.alphaCompositeBlocking) &&
                (depthPeeling == obj.depthPeeling) &&
                (occlusionRatio == obj.occlusionRatio) &&
                (numberOfPeels == obj.numberOfPeels) &&
                (multiresolutionMode == obj.multiresolutionMode) &&
                (multiresolutionCellSize == obj.multiresolutionCellSize) &&
                (geometryRepresentation == obj.geometryRepresentation) &&
                (stereoRendering == obj.stereoRendering) &&
                (stereoType == obj.stereoType) &&
                (notifyForEachRender == obj.notifyForEachRender) &&
                (scalableActivationMode == obj.scalableActivationMode) &&
                (scalableAutoThreshold == obj.scalableAutoThreshold) &&
                (specularFlag == obj.specularFlag) &&
                (specularCoeff == obj.specularCoeff) &&
                (specularPower == obj.specularPower) &&
                (specularColor == obj.specularColor) &&
                (doShadowing == obj.doShadowing) &&
                (shadowStrength == obj.shadowStrength) &&
                (doDepthCueing == obj.doDepthCueing) &&
                (depthCueingAutomatic == obj.depthCueingAutomatic) &&
                startCuePoint_equal &&
                endCuePoint_equal &&
                (compressionActivationMode == obj.compressionActivationMode) &&
                (colorTexturingFlag == obj.colorTexturingFlag) &&
                (compactDomainsActivationMode == obj.compactDomainsActivationMode) &&
                (compactDomainsAutoThreshold == obj.compactDomainsAutoThreshold) &&
                (osprayRendering == obj.osprayRendering) &&
                (ospraySPP == obj.ospraySPP) &&
                (osprayAO == obj.osprayAO) &&
                (osprayShadows == obj.osprayShadows));
    }

    // Property setting methods
    public void SetAntialiasing(boolean antialiasing_)
    {
        antialiasing = antialiasing_;
        Select(0);
    }

    public void SetOrderComposite(boolean orderComposite_)
    {
        orderComposite = orderComposite_;
        Select(1);
    }

    public void SetDepthCompositeThreads(int depthCompositeThreads_)
    {
        depthCompositeThreads = depthCompositeThreads_;
        Select(2);
    }

    public void SetDepthCompositeBlocking(int depthCompositeBlocking_)
    {
        depthCompositeBlocking = depthCompositeBlocking_;
        Select(3);
    }

    public void SetAlphaCompositeThreads(int alphaCompositeThreads_)
    {
        alphaCompositeThreads = alphaCompositeThreads_;
        Select(4);
    }

    public void SetAlphaCompositeBlocking(int alphaCompositeBlocking_)
    {
        alphaCompositeBlocking = alphaCompositeBlocking_;
        Select(5);
    }

    public void SetDepthPeeling(boolean depthPeeling_)
    {
        depthPeeling = depthPeeling_;
        Select(6);
    }

    public void SetOcclusionRatio(double occlusionRatio_)
    {
        occlusionRatio = occlusionRatio_;
        Select(7);
    }

    public void SetNumberOfPeels(int numberOfPeels_)
    {
        numberOfPeels = numberOfPeels_;
        Select(8);
    }

    public void SetMultiresolutionMode(boolean multiresolutionMode_)
    {
        multiresolutionMode = multiresolutionMode_;
        Select(9);
    }

    public void SetMultiresolutionCellSize(float multiresolutionCellSize_)
    {
        multiresolutionCellSize = multiresolutionCellSize_;
        Select(10);
    }

    public void SetGeometryRepresentation(int geometryRepresentation_)
    {
        geometryRepresentation = geometryRepresentation_;
        Select(11);
    }

    public void SetStereoRendering(boolean stereoRendering_)
    {
        stereoRendering = stereoRendering_;
        Select(12);
    }

    public void SetStereoType(int stereoType_)
    {
        stereoType = stereoType_;
        Select(13);
    }

    public void SetNotifyForEachRender(boolean notifyForEachRender_)
    {
        notifyForEachRender = notifyForEachRender_;
        Select(14);
    }

    public void SetScalableActivationMode(int scalableActivationMode_)
    {
        scalableActivationMode = scalableActivationMode_;
        Select(15);
    }

    public void SetScalableAutoThreshold(int scalableAutoThreshold_)
    {
        scalableAutoThreshold = scalableAutoThreshold_;
        Select(16);
    }

    public void SetSpecularFlag(boolean specularFlag_)
    {
        specularFlag = specularFlag_;
        Select(17);
    }

    public void SetSpecularCoeff(float specularCoeff_)
    {
        specularCoeff = specularCoeff_;
        Select(18);
    }

    public void SetSpecularPower(float specularPower_)
    {
        specularPower = specularPower_;
        Select(19);
    }

    public void SetSpecularColor(ColorAttribute specularColor_)
    {
        specularColor = specularColor_;
        Select(20);
    }

    public void SetDoShadowing(boolean doShadowing_)
    {
        doShadowing = doShadowing_;
        Select(21);
    }

    public void SetShadowStrength(double shadowStrength_)
    {
        shadowStrength = shadowStrength_;
        Select(22);
    }

    public void SetDoDepthCueing(boolean doDepthCueing_)
    {
        doDepthCueing = doDepthCueing_;
        Select(23);
    }

    public void SetDepthCueingAutomatic(boolean depthCueingAutomatic_)
    {
        depthCueingAutomatic = depthCueingAutomatic_;
        Select(24);
    }

    public void SetStartCuePoint(double[] startCuePoint_)
    {
        startCuePoint[0] = startCuePoint_[0];
        startCuePoint[1] = startCuePoint_[1];
        startCuePoint[2] = startCuePoint_[2];
        Select(25);
    }

    public void SetStartCuePoint(double e0, double e1, double e2)
    {
        startCuePoint[0] = e0;
        startCuePoint[1] = e1;
        startCuePoint[2] = e2;
        Select(25);
    }

    public void SetEndCuePoint(double[] endCuePoint_)
    {
        endCuePoint[0] = endCuePoint_[0];
        endCuePoint[1] = endCuePoint_[1];
        endCuePoint[2] = endCuePoint_[2];
        Select(26);
    }

    public void SetEndCuePoint(double e0, double e1, double e2)
    {
        endCuePoint[0] = e0;
        endCuePoint[1] = e1;
        endCuePoint[2] = e2;
        Select(26);
    }

    public void SetCompressionActivationMode(int compressionActivationMode_)
    {
        compressionActivationMode = compressionActivationMode_;
        Select(27);
    }

    public void SetColorTexturingFlag(boolean colorTexturingFlag_)
    {
        colorTexturingFlag = colorTexturingFlag_;
        Select(28);
    }

    public void SetCompactDomainsActivationMode(int compactDomainsActivationMode_)
    {
        compactDomainsActivationMode = compactDomainsActivationMode_;
        Select(29);
    }

    public void SetCompactDomainsAutoThreshold(int compactDomainsAutoThreshold_)
    {
        compactDomainsAutoThreshold = compactDomainsAutoThreshold_;
        Select(30);
    }

    public void SetOsprayRendering(boolean osprayRendering_)
    {
        osprayRendering = osprayRendering_;
        Select(31);
    }

    public void SetOspraySPP(int ospraySPP_)
    {
        ospraySPP = ospraySPP_;
        Select(32);
    }

    public void SetOsprayAO(int osprayAO_)
    {
        osprayAO = osprayAO_;
        Select(33);
    }

    public void SetOsprayShadows(boolean osprayShadows_)
    {
        osprayShadows = osprayShadows_;
        Select(34);
    }

    // Property getting methods
    public boolean        GetAntialiasing() { return antialiasing; }
    public boolean        GetOrderComposite() { return orderComposite; }
    public int            GetDepthCompositeThreads() { return depthCompositeThreads; }
    public int            GetDepthCompositeBlocking() { return depthCompositeBlocking; }
    public int            GetAlphaCompositeThreads() { return alphaCompositeThreads; }
    public int            GetAlphaCompositeBlocking() { return alphaCompositeBlocking; }
    public boolean        GetDepthPeeling() { return depthPeeling; }
    public double         GetOcclusionRatio() { return occlusionRatio; }
    public int            GetNumberOfPeels() { return numberOfPeels; }
    public boolean        GetMultiresolutionMode() { return multiresolutionMode; }
    public float          GetMultiresolutionCellSize() { return multiresolutionCellSize; }
    public int            GetGeometryRepresentation() { return geometryRepresentation; }
    public boolean        GetStereoRendering() { return stereoRendering; }
    public int            GetStereoType() { return stereoType; }
    public boolean        GetNotifyForEachRender() { return notifyForEachRender; }
    public int            GetScalableActivationMode() { return scalableActivationMode; }
    public int            GetScalableAutoThreshold() { return scalableAutoThreshold; }
    public boolean        GetSpecularFlag() { return specularFlag; }
    public float          GetSpecularCoeff() { return specularCoeff; }
    public float          GetSpecularPower() { return specularPower; }
    public ColorAttribute GetSpecularColor() { return specularColor; }
    public boolean        GetDoShadowing() { return doShadowing; }
    public double         GetShadowStrength() { return shadowStrength; }
    public boolean        GetDoDepthCueing() { return doDepthCueing; }
    public boolean        GetDepthCueingAutomatic() { return depthCueingAutomatic; }
    public double[]       GetStartCuePoint() { return startCuePoint; }
    public double[]       GetEndCuePoint() { return endCuePoint; }
    public int            GetCompressionActivationMode() { return compressionActivationMode; }
    public boolean        GetColorTexturingFlag() { return colorTexturingFlag; }
    public int            GetCompactDomainsActivationMode() { return compactDomainsActivationMode; }
    public int            GetCompactDomainsAutoThreshold() { return compactDomainsAutoThreshold; }
    public boolean        GetOsprayRendering() { return osprayRendering; }
    public int            GetOspraySPP() { return ospraySPP; }
    public int            GetOsprayAO() { return osprayAO; }
    public boolean        GetOsprayShadows() { return osprayShadows; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteBool(antialiasing);
        if(WriteSelect(1, buf))
            buf.WriteBool(orderComposite);
        if(WriteSelect(2, buf))
            buf.WriteInt(depthCompositeThreads);
        if(WriteSelect(3, buf))
            buf.WriteInt(depthCompositeBlocking);
        if(WriteSelect(4, buf))
            buf.WriteInt(alphaCompositeThreads);
        if(WriteSelect(5, buf))
            buf.WriteInt(alphaCompositeBlocking);
        if(WriteSelect(6, buf))
            buf.WriteBool(depthPeeling);
        if(WriteSelect(7, buf))
            buf.WriteDouble(occlusionRatio);
        if(WriteSelect(8, buf))
            buf.WriteInt(numberOfPeels);
        if(WriteSelect(9, buf))
            buf.WriteBool(multiresolutionMode);
        if(WriteSelect(10, buf))
            buf.WriteFloat(multiresolutionCellSize);
        if(WriteSelect(11, buf))
            buf.WriteInt(geometryRepresentation);
        if(WriteSelect(12, buf))
            buf.WriteBool(stereoRendering);
        if(WriteSelect(13, buf))
            buf.WriteInt(stereoType);
        if(WriteSelect(14, buf))
            buf.WriteBool(notifyForEachRender);
        if(WriteSelect(15, buf))
            buf.WriteInt(scalableActivationMode);
        if(WriteSelect(16, buf))
            buf.WriteInt(scalableAutoThreshold);
        if(WriteSelect(17, buf))
            buf.WriteBool(specularFlag);
        if(WriteSelect(18, buf))
            buf.WriteFloat(specularCoeff);
        if(WriteSelect(19, buf))
            buf.WriteFloat(specularPower);
        if(WriteSelect(20, buf))
            specularColor.Write(buf);
        if(WriteSelect(21, buf))
            buf.WriteBool(doShadowing);
        if(WriteSelect(22, buf))
            buf.WriteDouble(shadowStrength);
        if(WriteSelect(23, buf))
            buf.WriteBool(doDepthCueing);
        if(WriteSelect(24, buf))
            buf.WriteBool(depthCueingAutomatic);
        if(WriteSelect(25, buf))
            buf.WriteDoubleArray(startCuePoint);
        if(WriteSelect(26, buf))
            buf.WriteDoubleArray(endCuePoint);
        if(WriteSelect(27, buf))
            buf.WriteInt(compressionActivationMode);
        if(WriteSelect(28, buf))
            buf.WriteBool(colorTexturingFlag);
        if(WriteSelect(29, buf))
            buf.WriteInt(compactDomainsActivationMode);
        if(WriteSelect(30, buf))
            buf.WriteInt(compactDomainsAutoThreshold);
        if(WriteSelect(31, buf))
            buf.WriteBool(osprayRendering);
        if(WriteSelect(32, buf))
            buf.WriteInt(ospraySPP);
        if(WriteSelect(33, buf))
            buf.WriteInt(osprayAO);
        if(WriteSelect(34, buf))
            buf.WriteBool(osprayShadows);
    }

    public void ReadAtts(int index, CommunicationBuffer buf)
    {
        switch(index)
        {
        case 0:
            SetAntialiasing(buf.ReadBool());
            break;
        case 1:
            SetOrderComposite(buf.ReadBool());
            break;
        case 2:
            SetDepthCompositeThreads(buf.ReadInt());
            break;
        case 3:
            SetDepthCompositeBlocking(buf.ReadInt());
            break;
        case 4:
            SetAlphaCompositeThreads(buf.ReadInt());
            break;
        case 5:
            SetAlphaCompositeBlocking(buf.ReadInt());
            break;
        case 6:
            SetDepthPeeling(buf.ReadBool());
            break;
        case 7:
            SetOcclusionRatio(buf.ReadDouble());
            break;
        case 8:
            SetNumberOfPeels(buf.ReadInt());
            break;
        case 9:
            SetMultiresolutionMode(buf.ReadBool());
            break;
        case 10:
            SetMultiresolutionCellSize(buf.ReadFloat());
            break;
        case 11:
            SetGeometryRepresentation(buf.ReadInt());
            break;
        case 12:
            SetStereoRendering(buf.ReadBool());
            break;
        case 13:
            SetStereoType(buf.ReadInt());
            break;
        case 14:
            SetNotifyForEachRender(buf.ReadBool());
            break;
        case 15:
            SetScalableActivationMode(buf.ReadInt());
            break;
        case 16:
            SetScalableAutoThreshold(buf.ReadInt());
            break;
        case 17:
            SetSpecularFlag(buf.ReadBool());
            break;
        case 18:
            SetSpecularCoeff(buf.ReadFloat());
            break;
        case 19:
            SetSpecularPower(buf.ReadFloat());
            break;
        case 20:
            specularColor.Read(buf);
            Select(20);
            break;
        case 21:
            SetDoShadowing(buf.ReadBool());
            break;
        case 22:
            SetShadowStrength(buf.ReadDouble());
            break;
        case 23:
            SetDoDepthCueing(buf.ReadBool());
            break;
        case 24:
            SetDepthCueingAutomatic(buf.ReadBool());
            break;
        case 25:
            SetStartCuePoint(buf.ReadDoubleArray());
            break;
        case 26:
            SetEndCuePoint(buf.ReadDoubleArray());
            break;
        case 27:
            SetCompressionActivationMode(buf.ReadInt());
            break;
        case 28:
            SetColorTexturingFlag(buf.ReadBool());
            break;
        case 29:
            SetCompactDomainsActivationMode(buf.ReadInt());
            break;
        case 30:
            SetCompactDomainsAutoThreshold(buf.ReadInt());
            break;
        case 31:
            SetOsprayRendering(buf.ReadBool());
            break;
        case 32:
            SetOspraySPP(buf.ReadInt());
            break;
        case 33:
            SetOsprayAO(buf.ReadInt());
            break;
        case 34:
            SetOsprayShadows(buf.ReadBool());
            break;
        }
    }

    public String toString(String indent)
    {
        String str = new String();
        str = str + boolToString("antialiasing", antialiasing, indent) + "\n";
        str = str + boolToString("orderComposite", orderComposite, indent) + "\n";
        str = str + intToString("depthCompositeThreads", depthCompositeThreads, indent) + "\n";
        str = str + intToString("depthCompositeBlocking", depthCompositeBlocking, indent) + "\n";
        str = str + intToString("alphaCompositeThreads", alphaCompositeThreads, indent) + "\n";
        str = str + intToString("alphaCompositeBlocking", alphaCompositeBlocking, indent) + "\n";
        str = str + boolToString("depthPeeling", depthPeeling, indent) + "\n";
        str = str + doubleToString("occlusionRatio", occlusionRatio, indent) + "\n";
        str = str + intToString("numberOfPeels", numberOfPeels, indent) + "\n";
        str = str + boolToString("multiresolutionMode", multiresolutionMode, indent) + "\n";
        str = str + floatToString("multiresolutionCellSize", multiresolutionCellSize, indent) + "\n";
        str = str + indent + "geometryRepresentation = ";
        if(geometryRepresentation == GEOMETRYREPRESENTATION_SURFACES)
            str = str + "GEOMETRYREPRESENTATION_SURFACES";
        if(geometryRepresentation == GEOMETRYREPRESENTATION_WIREFRAME)
            str = str + "GEOMETRYREPRESENTATION_WIREFRAME";
        if(geometryRepresentation == GEOMETRYREPRESENTATION_POINTS)
            str = str + "GEOMETRYREPRESENTATION_POINTS";
        str = str + "\n";
        str = str + boolToString("stereoRendering", stereoRendering, indent) + "\n";
        str = str + indent + "stereoType = ";
        if(stereoType == STEREOTYPES_REDBLUE)
            str = str + "STEREOTYPES_REDBLUE";
        if(stereoType == STEREOTYPES_INTERLACED)
            str = str + "STEREOTYPES_INTERLACED";
        if(stereoType == STEREOTYPES_CRYSTALEYES)
            str = str + "STEREOTYPES_CRYSTALEYES";
        if(stereoType == STEREOTYPES_REDGREEN)
            str = str + "STEREOTYPES_REDGREEN";
        str = str + "\n";
        str = str + boolToString("notifyForEachRender", notifyForEachRender, indent) + "\n";
        str = str + indent + "scalableActivationMode = ";
        if(scalableActivationMode == TRISTATEMODE_NEVER)
            str = str + "TRISTATEMODE_NEVER";
        if(scalableActivationMode == TRISTATEMODE_ALWAYS)
            str = str + "TRISTATEMODE_ALWAYS";
        if(scalableActivationMode == TRISTATEMODE_AUTO)
            str = str + "TRISTATEMODE_AUTO";
        str = str + "\n";
        str = str + intToString("scalableAutoThreshold", scalableAutoThreshold, indent) + "\n";
        str = str + boolToString("specularFlag", specularFlag, indent) + "\n";
        str = str + floatToString("specularCoeff", specularCoeff, indent) + "\n";
        str = str + floatToString("specularPower", specularPower, indent) + "\n";
        str = str + indent + "specularColor = {" + specularColor.Red() + ", " + specularColor.Green() + ", " + specularColor.Blue() + ", " + specularColor.Alpha() + "}\n";
        str = str + boolToString("doShadowing", doShadowing, indent) + "\n";
        str = str + doubleToString("shadowStrength", shadowStrength, indent) + "\n";
        str = str + boolToString("doDepthCueing", doDepthCueing, indent) + "\n";
        str = str + boolToString("depthCueingAutomatic", depthCueingAutomatic, indent) + "\n";
        str = str + doubleArrayToString("startCuePoint", startCuePoint, indent) + "\n";
        str = str + doubleArrayToString("endCuePoint", endCuePoint, indent) + "\n";
        str = str + indent + "compressionActivationMode = ";
        if(compressionActivationMode == TRISTATEMODE_NEVER)
            str = str + "TRISTATEMODE_NEVER";
        if(compressionActivationMode == TRISTATEMODE_ALWAYS)
            str = str + "TRISTATEMODE_ALWAYS";
        if(compressionActivationMode == TRISTATEMODE_AUTO)
            str = str + "TRISTATEMODE_AUTO";
        str = str + "\n";
        str = str + boolToString("colorTexturingFlag", colorTexturingFlag, indent) + "\n";
        str = str + indent + "compactDomainsActivationMode = ";
        if(compactDomainsActivationMode == TRISTATEMODE_NEVER)
            str = str + "TRISTATEMODE_NEVER";
        if(compactDomainsActivationMode == TRISTATEMODE_ALWAYS)
            str = str + "TRISTATEMODE_ALWAYS";
        if(compactDomainsActivationMode == TRISTATEMODE_AUTO)
            str = str + "TRISTATEMODE_AUTO";
        str = str + "\n";
        str = str + intToString("compactDomainsAutoThreshold", compactDomainsAutoThreshold, indent) + "\n";
        str = str + boolToString("osprayRendering", osprayRendering, indent) + "\n";
        str = str + intToString("ospraySPP", ospraySPP, indent) + "\n";
        str = str + intToString("osprayAO", osprayAO, indent) + "\n";
        str = str + boolToString("osprayShadows", osprayShadows, indent) + "\n";
        return str;
    }


    // Attributes
    private boolean        antialiasing;
    private boolean        orderComposite;
    private int            depthCompositeThreads;
    private int            depthCompositeBlocking;
    private int            alphaCompositeThreads;
    private int            alphaCompositeBlocking;
    private boolean        depthPeeling;
    private double         occlusionRatio;
    private int            numberOfPeels;
    private boolean        multiresolutionMode;
    private float          multiresolutionCellSize;
    private int            geometryRepresentation;
    private boolean        stereoRendering;
    private int            stereoType;
    private boolean        notifyForEachRender;
    private int            scalableActivationMode;
    private int            scalableAutoThreshold;
    private boolean        specularFlag;
    private float          specularCoeff;
    private float          specularPower;
    private ColorAttribute specularColor;
    private boolean        doShadowing;
    private double         shadowStrength;
    private boolean        doDepthCueing;
    private boolean        depthCueingAutomatic;
    private double[]       startCuePoint;
    private double[]       endCuePoint;
    private int            compressionActivationMode;
    private boolean        colorTexturingFlag;
    private int            compactDomainsActivationMode;
    private int            compactDomainsAutoThreshold;
    private boolean        osprayRendering;
    private int            ospraySPP;
    private int            osprayAO;
    private boolean        osprayShadows;
}

