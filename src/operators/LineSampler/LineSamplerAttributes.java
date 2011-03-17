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

package llnl.visit.operators;

import llnl.visit.AttributeSubject;
import llnl.visit.CommunicationBuffer;
import llnl.visit.Plugin;

// ****************************************************************************
// Class: LineSamplerAttributes
//
// Purpose:
//    This class contains attributes for the line sampler operator.
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

public class LineSamplerAttributes extends AttributeSubject implements Plugin
{
    private static int LineSamplerAttributes_numAdditionalAtts = 19;

    // Enum values
    public final static int COORDINATESYSTEM_CARTESIAN = 0;
    public final static int COORDINATESYSTEM_CYLINDRICAL = 1;

    public final static int BEAMPROJECTION_PARALLEL = 0;
    public final static int BEAMPROJECTION_DIVERGENT = 1;

    public final static int BEAMAXIS_R = 0;
    public final static int BEAMAXIS_Z = 1;

    public final static int BEAMSHAPE_LINE = 0;
    public final static int BEAMSHAPE_CYLINDER = 1;
    public final static int BEAMSHAPE_CONE = 2;

    public final static int VIEWDIMENSION_ONE = 0;
    public final static int VIEWDIMENSION_TWO = 1;
    public final static int VIEWDIMENSION_THREE = 2;

    public final static int BEAMTYPE_TOPHAT = 0;
    public final static int BEAMTYPE_GAUSSIAN = 1;


    public LineSamplerAttributes()
    {
        super(LineSamplerAttributes_numAdditionalAtts);

        coordinateSystem = COORDINATESYSTEM_CYLINDRICAL;
        beamShape = BEAMSHAPE_LINE;
        radius = 0.1;
        divergence = 1;
        beamProjection = BEAMPROJECTION_PARALLEL;
        nBeams = 5;
        offset = 0.1;
        angle = 5;
        origin = new double[3];
        origin[0] = 0;
        origin[1] = 0;
        origin[2] = 0;
        beamAxis = BEAMAXIS_Z;
        poloialAngle = 0;
        poloialRTilt = 0;
        poloialZTilt = 0;
        toroialAngle = 0;
        viewDimension = VIEWDIMENSION_THREE;
        beamType = BEAMTYPE_TOPHAT;
        standardDeviation = 1;
        sampleDistance = 0.1;
        sampleArc = 10;
    }

    public LineSamplerAttributes(int nMoreFields)
    {
        super(LineSamplerAttributes_numAdditionalAtts + nMoreFields);

        coordinateSystem = COORDINATESYSTEM_CYLINDRICAL;
        beamShape = BEAMSHAPE_LINE;
        radius = 0.1;
        divergence = 1;
        beamProjection = BEAMPROJECTION_PARALLEL;
        nBeams = 5;
        offset = 0.1;
        angle = 5;
        origin = new double[3];
        origin[0] = 0;
        origin[1] = 0;
        origin[2] = 0;
        beamAxis = BEAMAXIS_Z;
        poloialAngle = 0;
        poloialRTilt = 0;
        poloialZTilt = 0;
        toroialAngle = 0;
        viewDimension = VIEWDIMENSION_THREE;
        beamType = BEAMTYPE_TOPHAT;
        standardDeviation = 1;
        sampleDistance = 0.1;
        sampleArc = 10;
    }

    public LineSamplerAttributes(LineSamplerAttributes obj)
    {
        super(LineSamplerAttributes_numAdditionalAtts);

        int i;

        coordinateSystem = obj.coordinateSystem;
        beamShape = obj.beamShape;
        radius = obj.radius;
        divergence = obj.divergence;
        beamProjection = obj.beamProjection;
        nBeams = obj.nBeams;
        offset = obj.offset;
        angle = obj.angle;
        origin = new double[3];
        origin[0] = obj.origin[0];
        origin[1] = obj.origin[1];
        origin[2] = obj.origin[2];

        beamAxis = obj.beamAxis;
        poloialAngle = obj.poloialAngle;
        poloialRTilt = obj.poloialRTilt;
        poloialZTilt = obj.poloialZTilt;
        toroialAngle = obj.toroialAngle;
        viewDimension = obj.viewDimension;
        beamType = obj.beamType;
        standardDeviation = obj.standardDeviation;
        sampleDistance = obj.sampleDistance;
        sampleArc = obj.sampleArc;

        SelectAll();
    }

    public int Offset()
    {
        return super.Offset() + super.GetNumAdditionalAttributes();
    }

    public int GetNumAdditionalAttributes()
    {
        return LineSamplerAttributes_numAdditionalAtts;
    }

    public boolean equals(LineSamplerAttributes obj)
    {
        int i;

        // Compare the origin arrays.
        boolean origin_equal = true;
        for(i = 0; i < 3 && origin_equal; ++i)
            origin_equal = (origin[i] == obj.origin[i]);

        // Create the return value
        return ((coordinateSystem == obj.coordinateSystem) &&
                (beamShape == obj.beamShape) &&
                (radius == obj.radius) &&
                (divergence == obj.divergence) &&
                (beamProjection == obj.beamProjection) &&
                (nBeams == obj.nBeams) &&
                (offset == obj.offset) &&
                (angle == obj.angle) &&
                origin_equal &&
                (beamAxis == obj.beamAxis) &&
                (poloialAngle == obj.poloialAngle) &&
                (poloialRTilt == obj.poloialRTilt) &&
                (poloialZTilt == obj.poloialZTilt) &&
                (toroialAngle == obj.toroialAngle) &&
                (viewDimension == obj.viewDimension) &&
                (beamType == obj.beamType) &&
                (standardDeviation == obj.standardDeviation) &&
                (sampleDistance == obj.sampleDistance) &&
                (sampleArc == obj.sampleArc));
    }

    public String GetName() { return "LineSampler"; }
    public String GetVersion() { return "1.0"; }

    // Property setting methods
    public void SetCoordinateSystem(int coordinateSystem_)
    {
        coordinateSystem = coordinateSystem_;
        Select(0);
    }

    public void SetBeamShape(int beamShape_)
    {
        beamShape = beamShape_;
        Select(1);
    }

    public void SetRadius(double radius_)
    {
        radius = radius_;
        Select(2);
    }

    public void SetDivergence(double divergence_)
    {
        divergence = divergence_;
        Select(3);
    }

    public void SetBeamProjection(int beamProjection_)
    {
        beamProjection = beamProjection_;
        Select(4);
    }

    public void SetNBeams(int nBeams_)
    {
        nBeams = nBeams_;
        Select(5);
    }

    public void SetOffset(double offset_)
    {
        offset = offset_;
        Select(6);
    }

    public void SetAngle(double angle_)
    {
        angle = angle_;
        Select(7);
    }

    public void SetOrigin(double[] origin_)
    {
        origin[0] = origin_[0];
        origin[1] = origin_[1];
        origin[2] = origin_[2];
        Select(8);
    }

    public void SetOrigin(double e0, double e1, double e2)
    {
        origin[0] = e0;
        origin[1] = e1;
        origin[2] = e2;
        Select(8);
    }

    public void SetBeamAxis(int beamAxis_)
    {
        beamAxis = beamAxis_;
        Select(9);
    }

    public void SetPoloialAngle(double poloialAngle_)
    {
        poloialAngle = poloialAngle_;
        Select(10);
    }

    public void SetPoloialRTilt(double poloialRTilt_)
    {
        poloialRTilt = poloialRTilt_;
        Select(11);
    }

    public void SetPoloialZTilt(double poloialZTilt_)
    {
        poloialZTilt = poloialZTilt_;
        Select(12);
    }

    public void SetToroialAngle(double toroialAngle_)
    {
        toroialAngle = toroialAngle_;
        Select(13);
    }

    public void SetViewDimension(int viewDimension_)
    {
        viewDimension = viewDimension_;
        Select(14);
    }

    public void SetBeamType(int beamType_)
    {
        beamType = beamType_;
        Select(15);
    }

    public void SetStandardDeviation(double standardDeviation_)
    {
        standardDeviation = standardDeviation_;
        Select(16);
    }

    public void SetSampleDistance(double sampleDistance_)
    {
        sampleDistance = sampleDistance_;
        Select(17);
    }

    public void SetSampleArc(double sampleArc_)
    {
        sampleArc = sampleArc_;
        Select(18);
    }

    // Property getting methods
    public int      GetCoordinateSystem() { return coordinateSystem; }
    public int      GetBeamShape() { return beamShape; }
    public double   GetRadius() { return radius; }
    public double   GetDivergence() { return divergence; }
    public int      GetBeamProjection() { return beamProjection; }
    public int      GetNBeams() { return nBeams; }
    public double   GetOffset() { return offset; }
    public double   GetAngle() { return angle; }
    public double[] GetOrigin() { return origin; }
    public int      GetBeamAxis() { return beamAxis; }
    public double   GetPoloialAngle() { return poloialAngle; }
    public double   GetPoloialRTilt() { return poloialRTilt; }
    public double   GetPoloialZTilt() { return poloialZTilt; }
    public double   GetToroialAngle() { return toroialAngle; }
    public int      GetViewDimension() { return viewDimension; }
    public int      GetBeamType() { return beamType; }
    public double   GetStandardDeviation() { return standardDeviation; }
    public double   GetSampleDistance() { return sampleDistance; }
    public double   GetSampleArc() { return sampleArc; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteInt(coordinateSystem);
        if(WriteSelect(1, buf))
            buf.WriteInt(beamShape);
        if(WriteSelect(2, buf))
            buf.WriteDouble(radius);
        if(WriteSelect(3, buf))
            buf.WriteDouble(divergence);
        if(WriteSelect(4, buf))
            buf.WriteInt(beamProjection);
        if(WriteSelect(5, buf))
            buf.WriteInt(nBeams);
        if(WriteSelect(6, buf))
            buf.WriteDouble(offset);
        if(WriteSelect(7, buf))
            buf.WriteDouble(angle);
        if(WriteSelect(8, buf))
            buf.WriteDoubleArray(origin);
        if(WriteSelect(9, buf))
            buf.WriteInt(beamAxis);
        if(WriteSelect(10, buf))
            buf.WriteDouble(poloialAngle);
        if(WriteSelect(11, buf))
            buf.WriteDouble(poloialRTilt);
        if(WriteSelect(12, buf))
            buf.WriteDouble(poloialZTilt);
        if(WriteSelect(13, buf))
            buf.WriteDouble(toroialAngle);
        if(WriteSelect(14, buf))
            buf.WriteInt(viewDimension);
        if(WriteSelect(15, buf))
            buf.WriteInt(beamType);
        if(WriteSelect(16, buf))
            buf.WriteDouble(standardDeviation);
        if(WriteSelect(17, buf))
            buf.WriteDouble(sampleDistance);
        if(WriteSelect(18, buf))
            buf.WriteDouble(sampleArc);
    }

    public void ReadAtts(int index, CommunicationBuffer buf)
    {
        switch(index)
        {
        case 0:
            SetCoordinateSystem(buf.ReadInt());
            break;
        case 1:
            SetBeamShape(buf.ReadInt());
            break;
        case 2:
            SetRadius(buf.ReadDouble());
            break;
        case 3:
            SetDivergence(buf.ReadDouble());
            break;
        case 4:
            SetBeamProjection(buf.ReadInt());
            break;
        case 5:
            SetNBeams(buf.ReadInt());
            break;
        case 6:
            SetOffset(buf.ReadDouble());
            break;
        case 7:
            SetAngle(buf.ReadDouble());
            break;
        case 8:
            SetOrigin(buf.ReadDoubleArray());
            break;
        case 9:
            SetBeamAxis(buf.ReadInt());
            break;
        case 10:
            SetPoloialAngle(buf.ReadDouble());
            break;
        case 11:
            SetPoloialRTilt(buf.ReadDouble());
            break;
        case 12:
            SetPoloialZTilt(buf.ReadDouble());
            break;
        case 13:
            SetToroialAngle(buf.ReadDouble());
            break;
        case 14:
            SetViewDimension(buf.ReadInt());
            break;
        case 15:
            SetBeamType(buf.ReadInt());
            break;
        case 16:
            SetStandardDeviation(buf.ReadDouble());
            break;
        case 17:
            SetSampleDistance(buf.ReadDouble());
            break;
        case 18:
            SetSampleArc(buf.ReadDouble());
            break;
        }
    }

    public String toString(String indent)
    {
        String str = new String();
        str = str + indent + "coordinateSystem = ";
        if(coordinateSystem == COORDINATESYSTEM_CARTESIAN)
            str = str + "COORDINATESYSTEM_CARTESIAN";
        if(coordinateSystem == COORDINATESYSTEM_CYLINDRICAL)
            str = str + "COORDINATESYSTEM_CYLINDRICAL";
        str = str + "\n";
        str = str + indent + "beamShape = ";
        if(beamShape == BEAMSHAPE_LINE)
            str = str + "BEAMSHAPE_LINE";
        if(beamShape == BEAMSHAPE_CYLINDER)
            str = str + "BEAMSHAPE_CYLINDER";
        if(beamShape == BEAMSHAPE_CONE)
            str = str + "BEAMSHAPE_CONE";
        str = str + "\n";
        str = str + doubleToString("radius", radius, indent) + "\n";
        str = str + doubleToString("divergence", divergence, indent) + "\n";
        str = str + indent + "beamProjection = ";
        if(beamProjection == BEAMPROJECTION_PARALLEL)
            str = str + "BEAMPROJECTION_PARALLEL";
        if(beamProjection == BEAMPROJECTION_DIVERGENT)
            str = str + "BEAMPROJECTION_DIVERGENT";
        str = str + "\n";
        str = str + intToString("nBeams", nBeams, indent) + "\n";
        str = str + doubleToString("offset", offset, indent) + "\n";
        str = str + doubleToString("angle", angle, indent) + "\n";
        str = str + doubleArrayToString("origin", origin, indent) + "\n";
        str = str + indent + "beamAxis = ";
        if(beamAxis == BEAMAXIS_R)
            str = str + "BEAMAXIS_R";
        if(beamAxis == BEAMAXIS_Z)
            str = str + "BEAMAXIS_Z";
        str = str + "\n";
        str = str + doubleToString("poloialAngle", poloialAngle, indent) + "\n";
        str = str + doubleToString("poloialRTilt", poloialRTilt, indent) + "\n";
        str = str + doubleToString("poloialZTilt", poloialZTilt, indent) + "\n";
        str = str + doubleToString("toroialAngle", toroialAngle, indent) + "\n";
        str = str + indent + "viewDimension = ";
        if(viewDimension == VIEWDIMENSION_ONE)
            str = str + "VIEWDIMENSION_ONE";
        if(viewDimension == VIEWDIMENSION_TWO)
            str = str + "VIEWDIMENSION_TWO";
        if(viewDimension == VIEWDIMENSION_THREE)
            str = str + "VIEWDIMENSION_THREE";
        str = str + "\n";
        str = str + indent + "beamType = ";
        if(beamType == BEAMTYPE_TOPHAT)
            str = str + "BEAMTYPE_TOPHAT";
        if(beamType == BEAMTYPE_GAUSSIAN)
            str = str + "BEAMTYPE_GAUSSIAN";
        str = str + "\n";
        str = str + doubleToString("standardDeviation", standardDeviation, indent) + "\n";
        str = str + doubleToString("sampleDistance", sampleDistance, indent) + "\n";
        str = str + doubleToString("sampleArc", sampleArc, indent) + "\n";
        return str;
    }


    // Attributes
    private int      coordinateSystem;
    private int      beamShape;
    private double   radius;
    private double   divergence;
    private int      beamProjection;
    private int      nBeams;
    private double   offset;
    private double   angle;
    private double[] origin;
    private int      beamAxis;
    private double   poloialAngle;
    private double   poloialRTilt;
    private double   poloialZTilt;
    private double   toroialAngle;
    private int      viewDimension;
    private int      beamType;
    private double   standardDeviation;
    private double   sampleDistance;
    private double   sampleArc;
}

