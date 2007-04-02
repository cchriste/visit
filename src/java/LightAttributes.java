// ***************************************************************************
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
// ***************************************************************************

package llnl.visit;


// ****************************************************************************
// Class: LightAttributes
//
// Purpose:
//    This class is a light in a light list.
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   Wed Mar 14 17:54:58 PST 2007
//
// Modifications:
//   
// ****************************************************************************

public class LightAttributes extends AttributeSubject
{
    // Enum values
    public final static int LIGHTTYPE_AMBIENT = 0;
    public final static int LIGHTTYPE_OBJECT = 1;
    public final static int LIGHTTYPE_CAMERA = 2;


    public LightAttributes()
    {
        super(6);

        enabledFlagCanBeToggled = true;
        enabledFlag = true;
        type = LIGHTTYPE_CAMERA;
        direction = new double[3];
        direction[0] = 0;
        direction[1] = 0;
        direction[2] = -1;
        color = new ColorAttribute(255, 255, 255);
        brightness = 1;
    }

    public LightAttributes(LightAttributes obj)
    {
        super(6);

        int i;

        enabledFlagCanBeToggled = obj.enabledFlagCanBeToggled;
        enabledFlag = obj.enabledFlag;
        type = obj.type;
        direction = new double[3];
        direction[0] = obj.direction[0];
        direction[1] = obj.direction[1];
        direction[2] = obj.direction[2];

        color = new ColorAttribute(obj.color);
        brightness = obj.brightness;

        SelectAll();
    }

    public boolean equals(LightAttributes obj)
    {
        int i;

        // Compare the direction arrays.
        boolean direction_equal = true;
        for(i = 0; i < 3 && direction_equal; ++i)
            direction_equal = (direction[i] == obj.direction[i]);

        // Create the return value
        return ((enabledFlagCanBeToggled == obj.enabledFlagCanBeToggled) &&
                (enabledFlag == obj.enabledFlag) &&
                (type == obj.type) &&
                direction_equal &&
                (color == obj.color) &&
                (brightness == obj.brightness));
    }

    // Property setting methods
    public void SetEnabledFlagCanBeToggled(boolean enabledFlagCanBeToggled_)
    {
        enabledFlagCanBeToggled = enabledFlagCanBeToggled_;
        Select(0);
    }

    public void SetEnabledFlag(boolean enabledFlag_)
    {
        enabledFlag = enabledFlag_;
        Select(1);
    }

    public void SetType(int type_)
    {
        type = type_;
        Select(2);
    }

    public void SetDirection(double[] direction_)
    {
        direction[0] = direction_[0];
        direction[1] = direction_[1];
        direction[2] = direction_[2];
        Select(3);
    }

    public void SetDirection(double e0, double e1, double e2)
    {
        direction[0] = e0;
        direction[1] = e1;
        direction[2] = e2;
        Select(3);
    }

    public void SetColor(ColorAttribute color_)
    {
        color = color_;
        Select(4);
    }

    public void SetBrightness(double brightness_)
    {
        brightness = brightness_;
        Select(5);
    }

    // Property getting methods
    public boolean        GetEnabledFlagCanBeToggled() { return enabledFlagCanBeToggled; }
    public boolean        GetEnabledFlag() { return enabledFlag; }
    public int            GetType() { return type; }
    public double[]       GetDirection() { return direction; }
    public ColorAttribute GetColor() { return color; }
    public double         GetBrightness() { return brightness; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteBool(enabledFlagCanBeToggled);
        if(WriteSelect(1, buf))
            buf.WriteBool(enabledFlag);
        if(WriteSelect(2, buf))
            buf.WriteInt(type);
        if(WriteSelect(3, buf))
            buf.WriteDoubleArray(direction);
        if(WriteSelect(4, buf))
            color.Write(buf);
        if(WriteSelect(5, buf))
            buf.WriteDouble(brightness);
    }

    public void ReadAtts(int n, CommunicationBuffer buf)
    {
        for(int i = 0; i < n; ++i)
        {
            int index = (int)buf.ReadByte();
            switch(index)
            {
            case 0:
                SetEnabledFlagCanBeToggled(buf.ReadBool());
                break;
            case 1:
                SetEnabledFlag(buf.ReadBool());
                break;
            case 2:
                SetType(buf.ReadInt());
                break;
            case 3:
                SetDirection(buf.ReadDoubleArray());
                break;
            case 4:
                color.Read(buf);
                Select(4);
                break;
            case 5:
                SetBrightness(buf.ReadDouble());
                break;
            }
        }
    }


    // Attributes
    private boolean        enabledFlagCanBeToggled;
    private boolean        enabledFlag;
    private int            type;
    private double[]       direction;
    private ColorAttribute color;
    private double         brightness;
}

