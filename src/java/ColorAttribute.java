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
// Creation:   Fri Aug 9 13:36:19 PST 2002
//
// Modifications:
//   Brad Whitlock, Fri Aug 9 13:51:19 PST 2002
//   I changed the class so it uses ints internally and reads and writes
//   itself as bytes.
//
// ****************************************************************************

public class ColorAttribute extends AttributeSubject
{
    public ColorAttribute()
    {
        super(1);
        color = new int[4];
        color[0] = 0;
        color[1] = 0;
        color[2] = 0;
        color[3] = 255;
    }

    public ColorAttribute(ColorAttribute obj)
    {
        super(1);
        int i;

        color = new int[4];
        for(i = 0; i < 4; ++i)
            color[i] = obj.color[i];

        SelectAll();
    }

    public ColorAttribute(int r, int g, int b)
    {
        super(1);
        color = new int[4];
        color[0] = r;
        color[1] = g;
        color[2] = b;
        color[3] = 255;
    }

    public ColorAttribute(int r, int g, int b, int a)
    {
        super(1);
        color = new int[4];
        color[0] = r;
        color[1] = g;
        color[2] = b;
        color[3] = a;
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
    public void SetColor(int[] color_)
    {
        color[0] = color_[0];
        color[1] = color_[1];
        color[2] = color_[2];
        color[3] = color_[3];
        Select(0);
    }

    // Property getting methods
    public int[] GetColor() { return color; }
    public int Red() { return color[0]; }
    public int Green() { return color[1]; }
    public int Blue() { return color[2]; }
    public int Alpha() { return color[3]; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
        {
            byte b[] = new byte[4];
            for(int i = 0; i < 4; ++i)
                b[i] = (byte)color[i];
            buf.WriteByteArray(b, true);
        }
    }

    public void ReadAtts(int n, CommunicationBuffer buf)
    {
        buf.ReadByte();
        byte[] b = buf.ReadByteArray();
        for(int i = 0; i < 4; ++i)
        {
            color[i] = (int)b[i];
            color[i] += ((color[i] < 0) ? 256 : 0);
        }
        Select(0);    
    }

    // Attributes
    private int[] color;
}

