// ****************************************************************************
//
// Copyright (c) 2000 - 2007, The Regents of the University of California
// Produced at the Lawrence Livermore National Laboratory
// All rights reserved.
//
// This file is part of VisIt. For details, see http://www.llnl.gov/visit/. The
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
//    documentation and/or materials provided with the distribution.
//  - Neither the name of the UC/LLNL nor  the names of its contributors may be
//    used to  endorse or  promote products derived from  this software without
//    specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT  HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR  IMPLIED WARRANTIES, INCLUDING,  BUT NOT  LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND  FITNESS FOR A PARTICULAR  PURPOSE
// ARE  DISCLAIMED.  IN  NO  EVENT  SHALL  THE  REGENTS  OF  THE  UNIVERSITY OF
// CALIFORNIA, THE U.S.  DEPARTMENT  OF  ENERGY OR CONTRIBUTORS BE  LIABLE  FOR
// ANY  DIRECT,  INDIRECT,  INCIDENTAL,  SPECIAL,  EXEMPLARY,  OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT  LIMITED TO, PROCUREMENT OF  SUBSTITUTE GOODS OR
// SERVICES; LOSS OF  USE, DATA, OR PROFITS; OR  BUSINESS INTERRUPTION) HOWEVER
// CAUSED  AND  ON  ANY  THEORY  OF  LIABILITY,  WHETHER  IN  CONTRACT,  STRICT
// LIABILITY, OR TORT  (INCLUDING NEGLIGENCE OR OTHERWISE)  ARISING IN ANY  WAY
// OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ****************************************************************************

package llnl.visit.operators;

import llnl.visit.AttributeSubject;
import llnl.visit.CommunicationBuffer;
import llnl.visit.Plugin;

// ****************************************************************************
// Class: SphereSliceAttributes
//
// Purpose:
//    This class contains attributes for the spherical slice.
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   Thu Jul 31 16:09:28 PST 2003
//
// Modifications:
//   
// ****************************************************************************

public class SphereSliceAttributes extends AttributeSubject implements Plugin
{
    public SphereSliceAttributes()
    {
        super(2);

        origin = new double[3];
        origin[0] = 0;
        origin[1] = 0;
        origin[2] = 0;
        radius = 1;
    }

    public SphereSliceAttributes(SphereSliceAttributes obj)
    {
        super(2);

        int i;

        origin = new double[3];
        origin[0] = obj.origin[0];
        origin[1] = obj.origin[1];
        origin[2] = obj.origin[2];

        radius = obj.radius;

        SelectAll();
    }

    public boolean equals(SphereSliceAttributes obj)
    {
        int i;

        // Compare the origin arrays.
        boolean origin_equal = true;
        for(i = 0; i < 3 && origin_equal; ++i)
            origin_equal = (origin[i] == obj.origin[i]);

        // Create the return value
        return (origin_equal &&
                (radius == obj.radius));
    }

    public String GetName() { return "SphereSlice"; }
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

    public void SetRadius(double radius_)
    {
        radius = radius_;
        Select(1);
    }

    // Property getting methods
    public double[] GetOrigin() { return origin; }
    public double   GetRadius() { return radius; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteDoubleArray(origin);
        if(WriteSelect(1, buf))
            buf.WriteDouble(radius);
    }

    public void ReadAtts(int n, CommunicationBuffer buf)
    {
        for(int i = 0; i < n; ++i)
        {
            int index = (int)buf.ReadByte();
            switch(index)
            {
            case 0:
                SetOrigin(buf.ReadDoubleArray());
                break;
            case 1:
                SetRadius(buf.ReadDouble());
                break;
            }
        }
    }


    // Attributes
    private double[] origin;
    private double   radius;
}

