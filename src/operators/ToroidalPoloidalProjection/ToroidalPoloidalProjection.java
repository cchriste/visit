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
// Class: ToroidalPoloidalProjection
//
// Purpose:
//    Projects Exterior of Torus from 3D to ToroidalPoloidal mapping in 2D
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

public class ToroidalPoloidalProjection extends AttributeSubject implements Plugin
{
    private static int ToroidalPoloidalProjection_numAdditionalAtts = 5;

    // Enum values
    public final static int CENTROIDSOURCE_MANUAL = 0;
    public final static int CENTROIDSOURCE_AUTO = 1;


    public ToroidalPoloidalProjection()
    {
        super(ToroidalPoloidalProjection_numAdditionalAtts);

        R0 = 1.5;
        r = 0.2;
        centroidSource = CENTROIDSOURCE_MANUAL;
        centroid = new double[3];
        centroid[0] = 0;
        centroid[1] = 0;
        centroid[2] = 0;
        project2D = true;
    }

    public ToroidalPoloidalProjection(int nMoreFields)
    {
        super(ToroidalPoloidalProjection_numAdditionalAtts + nMoreFields);

        R0 = 1.5;
        r = 0.2;
        centroidSource = CENTROIDSOURCE_MANUAL;
        centroid = new double[3];
        centroid[0] = 0;
        centroid[1] = 0;
        centroid[2] = 0;
        project2D = true;
    }

    public ToroidalPoloidalProjection(ToroidalPoloidalProjection obj)
    {
        super(ToroidalPoloidalProjection_numAdditionalAtts);

        int i;

        R0 = obj.R0;
        r = obj.r;
        centroidSource = obj.centroidSource;
        centroid = new double[3];
        centroid[0] = obj.centroid[0];
        centroid[1] = obj.centroid[1];
        centroid[2] = obj.centroid[2];

        project2D = obj.project2D;

        SelectAll();
    }

    public int Offset()
    {
        return super.Offset() + super.GetNumAdditionalAttributes();
    }

    public int GetNumAdditionalAttributes()
    {
        return ToroidalPoloidalProjection_numAdditionalAtts;
    }

    public boolean equals(ToroidalPoloidalProjection obj)
    {
        int i;

        // Compare the centroid arrays.
        boolean centroid_equal = true;
        for(i = 0; i < 3 && centroid_equal; ++i)
            centroid_equal = (centroid[i] == obj.centroid[i]);

        // Create the return value
        return ((R0 == obj.R0) &&
                (r == obj.r) &&
                (centroidSource == obj.centroidSource) &&
                centroid_equal &&
                (project2D == obj.project2D));
    }

    public String GetName() { return "ToroidalPoloidalProjection"; }
    public String GetVersion() { return "1.0"; }

    // Property setting methods
    public void SetR0(double R0_)
    {
        R0 = R0_;
        Select(0);
    }

    public void SetR(double r_)
    {
        r = r_;
        Select(1);
    }

    public void SetCentroidSource(int centroidSource_)
    {
        centroidSource = centroidSource_;
        Select(2);
    }

    public void SetCentroid(double[] centroid_)
    {
        centroid[0] = centroid_[0];
        centroid[1] = centroid_[1];
        centroid[2] = centroid_[2];
        Select(3);
    }

    public void SetCentroid(double e0, double e1, double e2)
    {
        centroid[0] = e0;
        centroid[1] = e1;
        centroid[2] = e2;
        Select(3);
    }

    public void SetProject2D(boolean project2D_)
    {
        project2D = project2D_;
        Select(4);
    }

    // Property getting methods
    public double   GetR0() { return R0; }
    public double   GetR() { return r; }
    public int      GetCentroidSource() { return centroidSource; }
    public double[] GetCentroid() { return centroid; }
    public boolean  GetProject2D() { return project2D; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteDouble(R0);
        if(WriteSelect(1, buf))
            buf.WriteDouble(r);
        if(WriteSelect(2, buf))
            buf.WriteInt(centroidSource);
        if(WriteSelect(3, buf))
            buf.WriteDoubleArray(centroid);
        if(WriteSelect(4, buf))
            buf.WriteBool(project2D);
    }

    public void ReadAtts(int index, CommunicationBuffer buf)
    {
        switch(index)
        {
        case 0:
            SetR0(buf.ReadDouble());
            break;
        case 1:
            SetR(buf.ReadDouble());
            break;
        case 2:
            SetCentroidSource(buf.ReadInt());
            break;
        case 3:
            SetCentroid(buf.ReadDoubleArray());
            break;
        case 4:
            SetProject2D(buf.ReadBool());
            break;
        }
    }

    public String toString(String indent)
    {
        String str = new String();
        str = str + doubleToString("R0", R0, indent) + "\n";
        str = str + doubleToString("r", r, indent) + "\n";
        str = str + indent + "centroidSource = ";
        if(centroidSource == CENTROIDSOURCE_MANUAL)
            str = str + "CENTROIDSOURCE_MANUAL";
        if(centroidSource == CENTROIDSOURCE_AUTO)
            str = str + "CENTROIDSOURCE_AUTO";
        str = str + "\n";
        str = str + doubleArrayToString("centroid", centroid, indent) + "\n";
        str = str + boolToString("project2D", project2D, indent) + "\n";
        return str;
    }


    // Attributes
    private double   R0;
    private double   r;
    private int      centroidSource;
    private double[] centroid;
    private boolean  project2D;
}

