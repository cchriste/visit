// ***************************************************************************
//
// Copyright (c) 2000 - 2016, Lawrence Livermore National Security, LLC
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
// Class: GaussianControlPoint
//
// Purpose:
//    This class contains the information for a gaussian in the opacity bar.
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

public class GaussianControlPoint extends AttributeSubject
{
    private static int GaussianControlPoint_numAdditionalAtts = 5;

    public GaussianControlPoint()
    {
        super(GaussianControlPoint_numAdditionalAtts);

        x = 0f;
        height = 0f;
        width = 0.001f;
        xBias = 0f;
        yBias = 0f;
    }

    public GaussianControlPoint(int nMoreFields)
    {
        super(GaussianControlPoint_numAdditionalAtts + nMoreFields);

        x = 0f;
        height = 0f;
        width = 0.001f;
        xBias = 0f;
        yBias = 0f;
    }

    public GaussianControlPoint(GaussianControlPoint obj)
    {
        super(GaussianControlPoint_numAdditionalAtts);

        x = obj.x;
        height = obj.height;
        width = obj.width;
        xBias = obj.xBias;
        yBias = obj.yBias;

        SelectAll();
    }

    public int Offset()
    {
        return super.Offset() + super.GetNumAdditionalAttributes();
    }

    public int GetNumAdditionalAttributes()
    {
        return GaussianControlPoint_numAdditionalAtts;
    }

    public boolean equals(GaussianControlPoint obj)
    {
        // Create the return value
        return ((x == obj.x) &&
                (height == obj.height) &&
                (width == obj.width) &&
                (xBias == obj.xBias) &&
                (yBias == obj.yBias));
    }

    // Property setting methods
    public void SetX(float x_)
    {
        x = x_;
        Select(0);
    }

    public void SetHeight(float height_)
    {
        height = height_;
        Select(1);
    }

    public void SetWidth(float width_)
    {
        width = width_;
        Select(2);
    }

    public void SetXBias(float xBias_)
    {
        xBias = xBias_;
        Select(3);
    }

    public void SetYBias(float yBias_)
    {
        yBias = yBias_;
        Select(4);
    }

    // Property getting methods
    public float GetX() { return x; }
    public float GetHeight() { return height; }
    public float GetWidth() { return width; }
    public float GetXBias() { return xBias; }
    public float GetYBias() { return yBias; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteFloat(x);
        if(WriteSelect(1, buf))
            buf.WriteFloat(height);
        if(WriteSelect(2, buf))
            buf.WriteFloat(width);
        if(WriteSelect(3, buf))
            buf.WriteFloat(xBias);
        if(WriteSelect(4, buf))
            buf.WriteFloat(yBias);
    }

    public void ReadAtts(int index, CommunicationBuffer buf)
    {
        switch(index)
        {
        case 0:
            SetX(buf.ReadFloat());
            break;
        case 1:
            SetHeight(buf.ReadFloat());
            break;
        case 2:
            SetWidth(buf.ReadFloat());
            break;
        case 3:
            SetXBias(buf.ReadFloat());
            break;
        case 4:
            SetYBias(buf.ReadFloat());
            break;
        }
    }

    public String toString(String indent)
    {
        String str = new String();
        str = str + floatToString("x", x, indent) + "\n";
        str = str + floatToString("height", height, indent) + "\n";
        str = str + floatToString("width", width, indent) + "\n";
        str = str + floatToString("xBias", xBias, indent) + "\n";
        str = str + floatToString("yBias", yBias, indent) + "\n";
        return str;
    }


    // Attributes
    private float x;
    private float height;
    private float width;
    private float xBias;
    private float yBias;
}

