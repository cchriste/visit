// ***************************************************************************
//
// Copyright (c) 2000 - 2017, Lawrence Livermore National Security, LLC
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
// Class: EllipsoidSliceAttributes
//
// Purpose:
//    EllipsoidSliceAttributes
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

public class EllipsoidSliceAttributes extends AttributeSubject implements Plugin
{
    private static int EllipsoidSliceAttributes_numAdditionalAtts = 3;

    public EllipsoidSliceAttributes()
    {
        super(EllipsoidSliceAttributes_numAdditionalAtts);

        origin = new double[3];
        origin[0] = 0;
        origin[1] = 0;
        origin[2] = 0;
        radii = new double[3];
        radii[0] = 1;
        radii[1] = 1;
        radii[2] = 1;
        rotationAngle = new double[3];
        rotationAngle[0] = 0;
        rotationAngle[1] = 0;
        rotationAngle[2] = 0;
    }

    public EllipsoidSliceAttributes(int nMoreFields)
    {
        super(EllipsoidSliceAttributes_numAdditionalAtts + nMoreFields);

        origin = new double[3];
        origin[0] = 0;
        origin[1] = 0;
        origin[2] = 0;
        radii = new double[3];
        radii[0] = 1;
        radii[1] = 1;
        radii[2] = 1;
        rotationAngle = new double[3];
        rotationAngle[0] = 0;
        rotationAngle[1] = 0;
        rotationAngle[2] = 0;
    }

    public EllipsoidSliceAttributes(EllipsoidSliceAttributes obj)
    {
        super(obj);

        int i;

        origin = new double[3];
        origin[0] = obj.origin[0];
        origin[1] = obj.origin[1];
        origin[2] = obj.origin[2];

        radii = new double[3];
        radii[0] = obj.radii[0];
        radii[1] = obj.radii[1];
        radii[2] = obj.radii[2];

        rotationAngle = new double[3];
        rotationAngle[0] = obj.rotationAngle[0];
        rotationAngle[1] = obj.rotationAngle[1];
        rotationAngle[2] = obj.rotationAngle[2];


        SelectAll();
    }

    public int Offset()
    {
        return super.Offset() + super.GetNumAdditionalAttributes();
    }

    public int GetNumAdditionalAttributes()
    {
        return EllipsoidSliceAttributes_numAdditionalAtts;
    }

    public boolean equals(EllipsoidSliceAttributes obj)
    {
        int i;

        // Compare the origin arrays.
        boolean origin_equal = true;
        for(i = 0; i < 3 && origin_equal; ++i)
            origin_equal = (origin[i] == obj.origin[i]);

        // Compare the radii arrays.
        boolean radii_equal = true;
        for(i = 0; i < 3 && radii_equal; ++i)
            radii_equal = (radii[i] == obj.radii[i]);

        // Compare the rotationAngle arrays.
        boolean rotationAngle_equal = true;
        for(i = 0; i < 3 && rotationAngle_equal; ++i)
            rotationAngle_equal = (rotationAngle[i] == obj.rotationAngle[i]);

        // Create the return value
        return (origin_equal &&
                radii_equal &&
                rotationAngle_equal);
    }

    public String GetName() { return "EllipsoidSlice"; }
    public String GetVersion() { return "1.0"; }

    // Property setting methods
    public void SetOrigin(double[] origin_)
    {
        origin[0] = origin_[0];
        origin[1] = origin_[1];
        origin[2] = origin_[2];
        Select(0);
    }

    public void SetOrigin(double e0, double e1, double e2)
    {
        origin[0] = e0;
        origin[1] = e1;
        origin[2] = e2;
        Select(0);
    }

    public void SetRadii(double[] radii_)
    {
        radii[0] = radii_[0];
        radii[1] = radii_[1];
        radii[2] = radii_[2];
        Select(1);
    }

    public void SetRadii(double e0, double e1, double e2)
    {
        radii[0] = e0;
        radii[1] = e1;
        radii[2] = e2;
        Select(1);
    }

    public void SetRotationAngle(double[] rotationAngle_)
    {
        rotationAngle[0] = rotationAngle_[0];
        rotationAngle[1] = rotationAngle_[1];
        rotationAngle[2] = rotationAngle_[2];
        Select(2);
    }

    public void SetRotationAngle(double e0, double e1, double e2)
    {
        rotationAngle[0] = e0;
        rotationAngle[1] = e1;
        rotationAngle[2] = e2;
        Select(2);
    }

    // Property getting methods
    public double[] GetOrigin() { return origin; }
    public double[] GetRadii() { return radii; }
    public double[] GetRotationAngle() { return rotationAngle; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteDoubleArray(origin);
        if(WriteSelect(1, buf))
            buf.WriteDoubleArray(radii);
        if(WriteSelect(2, buf))
            buf.WriteDoubleArray(rotationAngle);
    }

    public void ReadAtts(int index, CommunicationBuffer buf)
    {
        switch(index)
        {
        case 0:
            SetOrigin(buf.ReadDoubleArray());
            break;
        case 1:
            SetRadii(buf.ReadDoubleArray());
            break;
        case 2:
            SetRotationAngle(buf.ReadDoubleArray());
            break;
        }
    }

    public String toString(String indent)
    {
        String str = new String();
        str = str + doubleArrayToString("origin", origin, indent) + "\n";
        str = str + doubleArrayToString("radii", radii, indent) + "\n";
        str = str + doubleArrayToString("rotationAngle", rotationAngle, indent) + "\n";
        return str;
    }


    // Attributes
    private double[] origin;
    private double[] radii;
    private double[] rotationAngle;
}

