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


// ****************************************************************************
// Class: MaterialAttributes
//
// Purpose:
//    Attributes to control material interface reconstruction
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

public class MaterialAttributes extends AttributeSubject
{
    // Enum values
    public final static int ALGORITHM_TETRAHEDRAL = 0;
    public final static int ALGORITHM_ZOOCLIPPING = 1;
    public final static int ALGORITHM_ISOVOLUME = 2;


    public MaterialAttributes()
    {
        super(11);

        smoothing = false;
        forceMIR = false;
        cleanZonesOnly = false;
        needValidConnectivity = false;
        algorithm = ALGORITHM_ZOOCLIPPING;
        iterationEnabled = false;
        numIterations = 5;
        iterationDamping = 0.4f;
        simplifyHeavilyMixedZones = false;
        maxMaterialsPerZone = 3;
        isoVolumeFraction = 0.5f;
    }

    public MaterialAttributes(MaterialAttributes obj)
    {
        super(11);

        smoothing = obj.smoothing;
        forceMIR = obj.forceMIR;
        cleanZonesOnly = obj.cleanZonesOnly;
        needValidConnectivity = obj.needValidConnectivity;
        algorithm = obj.algorithm;
        iterationEnabled = obj.iterationEnabled;
        numIterations = obj.numIterations;
        iterationDamping = obj.iterationDamping;
        simplifyHeavilyMixedZones = obj.simplifyHeavilyMixedZones;
        maxMaterialsPerZone = obj.maxMaterialsPerZone;
        isoVolumeFraction = obj.isoVolumeFraction;

        SelectAll();
    }

    public boolean equals(MaterialAttributes obj)
    {
        // Create the return value
        return ((smoothing == obj.smoothing) &&
                (forceMIR == obj.forceMIR) &&
                (cleanZonesOnly == obj.cleanZonesOnly) &&
                (needValidConnectivity == obj.needValidConnectivity) &&
                (algorithm == obj.algorithm) &&
                (iterationEnabled == obj.iterationEnabled) &&
                (numIterations == obj.numIterations) &&
                (iterationDamping == obj.iterationDamping) &&
                (simplifyHeavilyMixedZones == obj.simplifyHeavilyMixedZones) &&
                (maxMaterialsPerZone == obj.maxMaterialsPerZone) &&
                (isoVolumeFraction == obj.isoVolumeFraction));
    }

    // Property setting methods
    public void SetSmoothing(boolean smoothing_)
    {
        smoothing = smoothing_;
        Select(0);
    }

    public void SetForceMIR(boolean forceMIR_)
    {
        forceMIR = forceMIR_;
        Select(1);
    }

    public void SetCleanZonesOnly(boolean cleanZonesOnly_)
    {
        cleanZonesOnly = cleanZonesOnly_;
        Select(2);
    }

    public void SetNeedValidConnectivity(boolean needValidConnectivity_)
    {
        needValidConnectivity = needValidConnectivity_;
        Select(3);
    }

    public void SetAlgorithm(int algorithm_)
    {
        algorithm = algorithm_;
        Select(4);
    }

    public void SetIterationEnabled(boolean iterationEnabled_)
    {
        iterationEnabled = iterationEnabled_;
        Select(5);
    }

    public void SetNumIterations(int numIterations_)
    {
        numIterations = numIterations_;
        Select(6);
    }

    public void SetIterationDamping(float iterationDamping_)
    {
        iterationDamping = iterationDamping_;
        Select(7);
    }

    public void SetSimplifyHeavilyMixedZones(boolean simplifyHeavilyMixedZones_)
    {
        simplifyHeavilyMixedZones = simplifyHeavilyMixedZones_;
        Select(8);
    }

    public void SetMaxMaterialsPerZone(int maxMaterialsPerZone_)
    {
        maxMaterialsPerZone = maxMaterialsPerZone_;
        Select(9);
    }

    public void SetIsoVolumeFraction(float isoVolumeFraction_)
    {
        isoVolumeFraction = isoVolumeFraction_;
        Select(10);
    }

    // Property getting methods
    public boolean GetSmoothing() { return smoothing; }
    public boolean GetForceMIR() { return forceMIR; }
    public boolean GetCleanZonesOnly() { return cleanZonesOnly; }
    public boolean GetNeedValidConnectivity() { return needValidConnectivity; }
    public int     GetAlgorithm() { return algorithm; }
    public boolean GetIterationEnabled() { return iterationEnabled; }
    public int     GetNumIterations() { return numIterations; }
    public float   GetIterationDamping() { return iterationDamping; }
    public boolean GetSimplifyHeavilyMixedZones() { return simplifyHeavilyMixedZones; }
    public int     GetMaxMaterialsPerZone() { return maxMaterialsPerZone; }
    public float   GetIsoVolumeFraction() { return isoVolumeFraction; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteBool(smoothing);
        if(WriteSelect(1, buf))
            buf.WriteBool(forceMIR);
        if(WriteSelect(2, buf))
            buf.WriteBool(cleanZonesOnly);
        if(WriteSelect(3, buf))
            buf.WriteBool(needValidConnectivity);
        if(WriteSelect(4, buf))
            buf.WriteInt(algorithm);
        if(WriteSelect(5, buf))
            buf.WriteBool(iterationEnabled);
        if(WriteSelect(6, buf))
            buf.WriteInt(numIterations);
        if(WriteSelect(7, buf))
            buf.WriteFloat(iterationDamping);
        if(WriteSelect(8, buf))
            buf.WriteBool(simplifyHeavilyMixedZones);
        if(WriteSelect(9, buf))
            buf.WriteInt(maxMaterialsPerZone);
        if(WriteSelect(10, buf))
            buf.WriteFloat(isoVolumeFraction);
    }

    public void ReadAtts(int n, CommunicationBuffer buf)
    {
        for(int i = 0; i < n; ++i)
        {
            int index = (int)buf.ReadByte();
            switch(index)
            {
            case 0:
                SetSmoothing(buf.ReadBool());
                break;
            case 1:
                SetForceMIR(buf.ReadBool());
                break;
            case 2:
                SetCleanZonesOnly(buf.ReadBool());
                break;
            case 3:
                SetNeedValidConnectivity(buf.ReadBool());
                break;
            case 4:
                SetAlgorithm(buf.ReadInt());
                break;
            case 5:
                SetIterationEnabled(buf.ReadBool());
                break;
            case 6:
                SetNumIterations(buf.ReadInt());
                break;
            case 7:
                SetIterationDamping(buf.ReadFloat());
                break;
            case 8:
                SetSimplifyHeavilyMixedZones(buf.ReadBool());
                break;
            case 9:
                SetMaxMaterialsPerZone(buf.ReadInt());
                break;
            case 10:
                SetIsoVolumeFraction(buf.ReadFloat());
                break;
            }
        }
    }

    public String toString(String indent)
    {
        String str = new String();
        str = str + boolToString("smoothing", smoothing, indent) + "\n";
        str = str + boolToString("forceMIR", forceMIR, indent) + "\n";
        str = str + boolToString("cleanZonesOnly", cleanZonesOnly, indent) + "\n";
        str = str + boolToString("needValidConnectivity", needValidConnectivity, indent) + "\n";
        str = str + indent + "algorithm = ";
        if(algorithm == ALGORITHM_TETRAHEDRAL)
            str = str + "ALGORITHM_TETRAHEDRAL";
        if(algorithm == ALGORITHM_ZOOCLIPPING)
            str = str + "ALGORITHM_ZOOCLIPPING";
        if(algorithm == ALGORITHM_ISOVOLUME)
            str = str + "ALGORITHM_ISOVOLUME";
        str = str + "\n";
        str = str + boolToString("iterationEnabled", iterationEnabled, indent) + "\n";
        str = str + intToString("numIterations", numIterations, indent) + "\n";
        str = str + floatToString("iterationDamping", iterationDamping, indent) + "\n";
        str = str + boolToString("simplifyHeavilyMixedZones", simplifyHeavilyMixedZones, indent) + "\n";
        str = str + intToString("maxMaterialsPerZone", maxMaterialsPerZone, indent) + "\n";
        str = str + floatToString("isoVolumeFraction", isoVolumeFraction, indent) + "\n";
        return str;
    }


    // Attributes
    private boolean smoothing;
    private boolean forceMIR;
    private boolean cleanZonesOnly;
    private boolean needValidConnectivity;
    private int     algorithm;
    private boolean iterationEnabled;
    private int     numIterations;
    private float   iterationDamping;
    private boolean simplifyHeavilyMixedZones;
    private int     maxMaterialsPerZone;
    private float   isoVolumeFraction;
}

