// ***************************************************************************
//
// Copyright (c) 2000 - 2015, Lawrence Livermore National Security, LLC
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
// Class: ColorAttribute
//
// Purpose:
//    This class contains RGBA color information
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

public class ColorAttribute extends AttributeSubject
{
    private static int ColorAttribute_numAdditionalAtts = 1;


    public ColorAttribute()
    {
        super(ColorAttribute_numAdditionalAtts);

        color = new byte[4];
        color[0] = (byte)0;
        color[1] = (byte)0;
        color[2] = (byte)0;
        color[3] = (byte)255;
    }

    public ColorAttribute(int nMoreFields)
    {
        super(ColorAttribute_numAdditionalAtts + nMoreFields);

        color = new byte[4];
        color[0] = (byte)0;
        color[1] = (byte)0;
        color[2] = (byte)0;
        color[3] = (byte)255;
    }

    public ColorAttribute(ColorAttribute obj)
    {
        super(ColorAttribute_numAdditionalAtts);

        int i;

        color = new byte[4];
        for(i = 0; i < obj.color.length; ++i)
            color[i] = obj.color[i];


        SelectAll();
    }

    public int Offset()
    {
        return super.Offset() + super.GetNumAdditionalAttributes();
    }

    public int GetNumAdditionalAttributes()
    {
        return ColorAttribute_numAdditionalAtts;
    }

    public boolean equals(ColorAttribute obj)
    {
        int i;

        // Compare the color arrays.
        boolean color_equal = true;
        for(i = 0; i < 4 && color_equal; ++i)
            color_equal = (color[i] == obj.color[i]);

        // Create the return value
        return (color_equal);
    }

    // Property setting methods
    public void SetColor(byte[] color_)
    {
        color[0] = color_[0];
        color[1] = color_[1];
        color[2] = color_[2];
        color[3] = color_[3];
        Select(0);
    }

    public void SetColor(byte e0, byte e1, byte e2, byte e3)
    {
        color[0] = e0;
        color[1] = e1;
        color[2] = e2;
        color[3] = e3;
        Select(0);
    }

    // Property getting methods
    public byte[] GetColor() { return color; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteByteArray(color, true);
    }

    public void ReadAtts(int index, CommunicationBuffer buf)
    {
        SetColor(buf.ReadByteArray());
    }

    public String toString(String indent)
    {
        String str = new String();
        str = str + ucharArrayToString("color", color, indent) + "\n";
        return str;
    }


    public ColorAttribute(int r, int g, int b)
    {
        super(1);
        color = new byte[4];
        color[0] = (byte)(0xff & r);
        color[1] = (byte)(0xff & g);
        color[2] = (byte)(0xff & b);
        color[3] = (byte)(0xff & 255);
    }

    public ColorAttribute(int r, int g, int b, int a)
    {
        super(1);
        color = new byte[4];
        color[0] = (byte)(0xff & r);
        color[1] = (byte)(0xff & g);
        color[2] = (byte)(0xff & b);
        color[3] = (byte)(0xff & a);
    }

    public int Red()
    { 
       int mask = 0xff;
       return mask & color[0];
    }

    public int Green()
    { 
       int mask = 0xff;
       return mask & color[1];
    }

    public int Blue()
    { 
       int mask = 0xff;
       return mask & color[2];
    }

    public int Alpha()
    { 
       int mask = 0xff;
       return mask & color[3];
    }

    // Attributes
    private byte[] color;
}

