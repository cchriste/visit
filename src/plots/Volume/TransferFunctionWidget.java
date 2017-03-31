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

package llnl.visit;


// ****************************************************************************
// Class: TransferFunctionWidget
//
// Purpose:
//    Widget for 2D transfer function
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

public class TransferFunctionWidget extends AttributeSubject
{
    private static int TransferFunctionWidget_numAdditionalAtts = 4;

    // Enum values
    public final static int WIDGETTYPE_RECTANGLE = 0;
    public final static int WIDGETTYPE_TRIANGLE = 1;
    public final static int WIDGETTYPE_PARABOLOID = 2;
    public final static int WIDGETTYPE_ELLIPSOID = 3;


    public TransferFunctionWidget()
    {
        super(TransferFunctionWidget_numAdditionalAtts);

        Type = WIDGETTYPE_RECTANGLE;
        Name = new String("unnamed");
        BaseColor = new float[4];
        BaseColor[0] = 1f;
        BaseColor[1] = 1f;
        BaseColor[2] = 1f;
        BaseColor[3] = 1f;
        Position = new float[8];
        Position[0] = 0f;
        Position[1] = 0f;
        Position[2] = 0f;
        Position[3] = 0f;
        Position[4] = 0f;
        Position[5] = 0f;
        Position[6] = 0f;
        Position[7] = 0f;
    }

    public TransferFunctionWidget(int nMoreFields)
    {
        super(TransferFunctionWidget_numAdditionalAtts + nMoreFields);

        Type = WIDGETTYPE_RECTANGLE;
        Name = new String("unnamed");
        BaseColor = new float[4];
        BaseColor[0] = 1f;
        BaseColor[1] = 1f;
        BaseColor[2] = 1f;
        BaseColor[3] = 1f;
        Position = new float[8];
        Position[0] = 0f;
        Position[1] = 0f;
        Position[2] = 0f;
        Position[3] = 0f;
        Position[4] = 0f;
        Position[5] = 0f;
        Position[6] = 0f;
        Position[7] = 0f;
    }

    public TransferFunctionWidget(TransferFunctionWidget obj)
    {
        super(obj);

        int i;

        Type = obj.Type;
        Name = new String(obj.Name);
        BaseColor = new float[4];
        for(i = 0; i < obj.BaseColor.length; ++i)
            BaseColor[i] = obj.BaseColor[i];

        Position = new float[8];
        for(i = 0; i < obj.Position.length; ++i)
            Position[i] = obj.Position[i];


        SelectAll();
    }

    public int Offset()
    {
        return super.Offset() + super.GetNumAdditionalAttributes();
    }

    public int GetNumAdditionalAttributes()
    {
        return TransferFunctionWidget_numAdditionalAtts;
    }

    public boolean equals(TransferFunctionWidget obj)
    {
        int i;

        // Compare the BaseColor arrays.
        boolean BaseColor_equal = true;
        for(i = 0; i < 4 && BaseColor_equal; ++i)
            BaseColor_equal = (BaseColor[i] == obj.BaseColor[i]);

        // Compare the Position arrays.
        boolean Position_equal = true;
        for(i = 0; i < 8 && Position_equal; ++i)
            Position_equal = (Position[i] == obj.Position[i]);

        // Create the return value
        return ((Type == obj.Type) &&
                (Name.equals(obj.Name)) &&
                BaseColor_equal &&
                Position_equal);
    }

    // Property setting methods
    public void SetType(int Type_)
    {
        Type = Type_;
        Select(0);
    }

    public void SetName(String Name_)
    {
        Name = Name_;
        Select(1);
    }

    public void SetBaseColor(float[] BaseColor_)
    {
        BaseColor[0] = BaseColor_[0];
        BaseColor[1] = BaseColor_[1];
        BaseColor[2] = BaseColor_[2];
        BaseColor[3] = BaseColor_[3];
        Select(2);
    }

    public void SetBaseColor(float e0, float e1, float e2, float e3)
    {
        BaseColor[0] = e0;
        BaseColor[1] = e1;
        BaseColor[2] = e2;
        BaseColor[3] = e3;
        Select(2);
    }

    public void SetPosition(float[] Position_)
    {
        for(int i = 0; i < 8; ++i)
             Position[i] = Position_[i];
        Select(3);
    }

    // Property getting methods
    public int     GetType() { return Type; }
    public String  GetName() { return Name; }
    public float[] GetBaseColor() { return BaseColor; }
    public float[] GetPosition() { return Position; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteInt(Type);
        if(WriteSelect(1, buf))
            buf.WriteString(Name);
        if(WriteSelect(2, buf))
            buf.WriteFloatArray(BaseColor);
        if(WriteSelect(3, buf))
            buf.WriteFloatArray(Position);
    }

    public void ReadAtts(int index, CommunicationBuffer buf)
    {
        switch(index)
        {
        case 0:
            SetType(buf.ReadInt());
            break;
        case 1:
            SetName(buf.ReadString());
            break;
        case 2:
            SetBaseColor(buf.ReadFloatArray());
            break;
        case 3:
            SetPosition(buf.ReadFloatArray());
            break;
        }
    }

    public String toString(String indent)
    {
        String str = new String();
        str = str + indent + "Type = ";
        if(Type == WIDGETTYPE_RECTANGLE)
            str = str + "WIDGETTYPE_RECTANGLE";
        if(Type == WIDGETTYPE_TRIANGLE)
            str = str + "WIDGETTYPE_TRIANGLE";
        if(Type == WIDGETTYPE_PARABOLOID)
            str = str + "WIDGETTYPE_PARABOLOID";
        if(Type == WIDGETTYPE_ELLIPSOID)
            str = str + "WIDGETTYPE_ELLIPSOID";
        str = str + "\n";
        str = str + stringToString("Name", Name, indent) + "\n";
        str = str + floatArrayToString("BaseColor", BaseColor, indent) + "\n";
        str = str + floatArrayToString("Position", Position, indent) + "\n";
        return str;
    }


    // Attributes
    private int     Type;
    private String  Name;
    private float[] BaseColor;
    private float[] Position;
}

