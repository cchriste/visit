// ***************************************************************************
//
// Copyright (c) 2000 - 2018, Lawrence Livermore National Security, LLC
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

import java.util.Vector;

// ****************************************************************************
// Class: AnnotationObject
//
// Purpose:
//    This class defines a general set of attributes that are used to set the attributes for all annotation objects.
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

public class AnnotationObject extends AttributeSubject
{
    private static int AnnotationObject_numAdditionalAtts = 18;

    // Enum values
    public final static int ANNOTATIONTYPE_TEXT2D = 0;
    public final static int ANNOTATIONTYPE_TEXT3D = 1;
    public final static int ANNOTATIONTYPE_TIMESLIDER = 2;
    public final static int ANNOTATIONTYPE_LINE2D = 3;
    public final static int ANNOTATIONTYPE_LINE3D = 4;
    public final static int ANNOTATIONTYPE_ARROW2D = 5;
    public final static int ANNOTATIONTYPE_ARROW3D = 6;
    public final static int ANNOTATIONTYPE_BOX = 7;
    public final static int ANNOTATIONTYPE_IMAGE = 8;
    public final static int ANNOTATIONTYPE_LEGENDATTRIBUTES = 9;
    public final static int ANNOTATIONTYPE_MAXANNOTATIONTYPE = 10;

    public final static int FONTFAMILY_ARIAL = 0;
    public final static int FONTFAMILY_COURIER = 1;
    public final static int FONTFAMILY_TIMES = 2;


    public AnnotationObject()
    {
        super(AnnotationObject_numAdditionalAtts);

        objectName = new String("");
        objectType = ANNOTATIONTYPE_TEXT2D;
        visible = false;
        active = false;
        position = new double[3];
        position[0] = 0;
        position[1] = 0;
        position[2] = 0;
        position2 = new double[3];
        position2[0] = 0;
        position2[1] = 0;
        position2[2] = 0;
        textColor = new ColorAttribute();
        useForegroundForTextColor = true;
        color1 = new ColorAttribute();
        color2 = new ColorAttribute();
        text = new Vector();
        fontFamily = FONTFAMILY_ARIAL;
        fontBold = false;
        fontItalic = false;
        fontShadow = false;
        doubleAttribute1 = 0;
        intAttribute1 = 0;
        options = new MapNode();
    }

    public AnnotationObject(int nMoreFields)
    {
        super(AnnotationObject_numAdditionalAtts + nMoreFields);

        objectName = new String("");
        objectType = ANNOTATIONTYPE_TEXT2D;
        visible = false;
        active = false;
        position = new double[3];
        position[0] = 0;
        position[1] = 0;
        position[2] = 0;
        position2 = new double[3];
        position2[0] = 0;
        position2[1] = 0;
        position2[2] = 0;
        textColor = new ColorAttribute();
        useForegroundForTextColor = true;
        color1 = new ColorAttribute();
        color2 = new ColorAttribute();
        text = new Vector();
        fontFamily = FONTFAMILY_ARIAL;
        fontBold = false;
        fontItalic = false;
        fontShadow = false;
        doubleAttribute1 = 0;
        intAttribute1 = 0;
        options = new MapNode();
    }

    public AnnotationObject(AnnotationObject obj)
    {
        super(obj);

        int i;

        objectName = new String(obj.objectName);
        objectType = obj.objectType;
        visible = obj.visible;
        active = obj.active;
        position = new double[3];
        position[0] = obj.position[0];
        position[1] = obj.position[1];
        position[2] = obj.position[2];

        position2 = new double[3];
        position2[0] = obj.position2[0];
        position2[1] = obj.position2[1];
        position2[2] = obj.position2[2];

        textColor = new ColorAttribute(obj.textColor);
        useForegroundForTextColor = obj.useForegroundForTextColor;
        color1 = new ColorAttribute(obj.color1);
        color2 = new ColorAttribute(obj.color2);
        text = new Vector(obj.text.size());
        for(i = 0; i < obj.text.size(); ++i)
            text.addElement(new String((String)obj.text.elementAt(i)));

        fontFamily = obj.fontFamily;
        fontBold = obj.fontBold;
        fontItalic = obj.fontItalic;
        fontShadow = obj.fontShadow;
        doubleAttribute1 = obj.doubleAttribute1;
        intAttribute1 = obj.intAttribute1;
        options = new MapNode(obj.options);

        SelectAll();
    }

    public int Offset()
    {
        return super.Offset() + super.GetNumAdditionalAttributes();
    }

    public int GetNumAdditionalAttributes()
    {
        return AnnotationObject_numAdditionalAtts;
    }

    public boolean equals(AnnotationObject obj)
    {
        int i;

        // Compare the position arrays.
        boolean position_equal = true;
        for(i = 0; i < 3 && position_equal; ++i)
            position_equal = (position[i] == obj.position[i]);

        // Compare the position2 arrays.
        boolean position2_equal = true;
        for(i = 0; i < 3 && position2_equal; ++i)
            position2_equal = (position2[i] == obj.position2[i]);

        // Compare the elements in the text vector.
        boolean text_equal = (obj.text.size() == text.size());
        for(i = 0; (i < text.size()) && text_equal; ++i)
        {
            // Make references to String from Object.
            String text1 = (String)text.elementAt(i);
            String text2 = (String)obj.text.elementAt(i);
            text_equal = text1.equals(text2);
        }
        // Create the return value
        return ((objectName.equals(obj.objectName)) &&
                (objectType == obj.objectType) &&
                (visible == obj.visible) &&
                (active == obj.active) &&
                position_equal &&
                position2_equal &&
                (textColor == obj.textColor) &&
                (useForegroundForTextColor == obj.useForegroundForTextColor) &&
                (color1 == obj.color1) &&
                (color2 == obj.color2) &&
                text_equal &&
                (fontFamily == obj.fontFamily) &&
                (fontBold == obj.fontBold) &&
                (fontItalic == obj.fontItalic) &&
                (fontShadow == obj.fontShadow) &&
                (doubleAttribute1 == obj.doubleAttribute1) &&
                (intAttribute1 == obj.intAttribute1) &&
                (options.equals(obj.options)));
    }

    // Property setting methods
    public void SetObjectName(String objectName_)
    {
        objectName = objectName_;
        Select(0);
    }

    public void SetObjectType(int objectType_)
    {
        objectType = objectType_;
        Select(1);
    }

    public void SetVisible(boolean visible_)
    {
        visible = visible_;
        Select(2);
    }

    public void SetActive(boolean active_)
    {
        active = active_;
        Select(3);
    }

    public void SetPosition(double[] position_)
    {
        position[0] = position_[0];
        position[1] = position_[1];
        position[2] = position_[2];
        Select(4);
    }

    public void SetPosition(double e0, double e1, double e2)
    {
        position[0] = e0;
        position[1] = e1;
        position[2] = e2;
        Select(4);
    }

    public void SetPosition2(double[] position2_)
    {
        position2[0] = position2_[0];
        position2[1] = position2_[1];
        position2[2] = position2_[2];
        Select(5);
    }

    public void SetPosition2(double e0, double e1, double e2)
    {
        position2[0] = e0;
        position2[1] = e1;
        position2[2] = e2;
        Select(5);
    }

    public void SetTextColor(ColorAttribute textColor_)
    {
        textColor = textColor_;
        Select(6);
    }

    public void SetUseForegroundForTextColor(boolean useForegroundForTextColor_)
    {
        useForegroundForTextColor = useForegroundForTextColor_;
        Select(7);
    }

    public void SetColor1(ColorAttribute color1_)
    {
        color1 = color1_;
        Select(8);
    }

    public void SetColor2(ColorAttribute color2_)
    {
        color2 = color2_;
        Select(9);
    }

    public void SetText(Vector text_)
    {
        text = text_;
        Select(10);
    }

    public void SetFontFamily(int fontFamily_)
    {
        fontFamily = fontFamily_;
        Select(11);
    }

    public void SetFontBold(boolean fontBold_)
    {
        fontBold = fontBold_;
        Select(12);
    }

    public void SetFontItalic(boolean fontItalic_)
    {
        fontItalic = fontItalic_;
        Select(13);
    }

    public void SetFontShadow(boolean fontShadow_)
    {
        fontShadow = fontShadow_;
        Select(14);
    }

    public void SetDoubleAttribute1(double doubleAttribute1_)
    {
        doubleAttribute1 = doubleAttribute1_;
        Select(15);
    }

    public void SetIntAttribute1(int intAttribute1_)
    {
        intAttribute1 = intAttribute1_;
        Select(16);
    }

    public void SetOptions(MapNode options_)
    {
        options = options_;
        Select(17);
    }

    // Property getting methods
    public String         GetObjectName() { return objectName; }
    public int            GetObjectType() { return objectType; }
    public boolean        GetVisible() { return visible; }
    public boolean        GetActive() { return active; }
    public double[]       GetPosition() { return position; }
    public double[]       GetPosition2() { return position2; }
    public ColorAttribute GetTextColor() { return textColor; }
    public boolean        GetUseForegroundForTextColor() { return useForegroundForTextColor; }
    public ColorAttribute GetColor1() { return color1; }
    public ColorAttribute GetColor2() { return color2; }
    public Vector         GetText() { return text; }
    public int            GetFontFamily() { return fontFamily; }
    public boolean        GetFontBold() { return fontBold; }
    public boolean        GetFontItalic() { return fontItalic; }
    public boolean        GetFontShadow() { return fontShadow; }
    public double         GetDoubleAttribute1() { return doubleAttribute1; }
    public int            GetIntAttribute1() { return intAttribute1; }
    public MapNode        GetOptions() { return options; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteString(objectName);
        if(WriteSelect(1, buf))
            buf.WriteInt(objectType);
        if(WriteSelect(2, buf))
            buf.WriteBool(visible);
        if(WriteSelect(3, buf))
            buf.WriteBool(active);
        if(WriteSelect(4, buf))
            buf.WriteDoubleArray(position);
        if(WriteSelect(5, buf))
            buf.WriteDoubleArray(position2);
        if(WriteSelect(6, buf))
            textColor.Write(buf);
        if(WriteSelect(7, buf))
            buf.WriteBool(useForegroundForTextColor);
        if(WriteSelect(8, buf))
            color1.Write(buf);
        if(WriteSelect(9, buf))
            color2.Write(buf);
        if(WriteSelect(10, buf))
            buf.WriteStringVector(text);
        if(WriteSelect(11, buf))
            buf.WriteInt(fontFamily);
        if(WriteSelect(12, buf))
            buf.WriteBool(fontBold);
        if(WriteSelect(13, buf))
            buf.WriteBool(fontItalic);
        if(WriteSelect(14, buf))
            buf.WriteBool(fontShadow);
        if(WriteSelect(15, buf))
            buf.WriteDouble(doubleAttribute1);
        if(WriteSelect(16, buf))
            buf.WriteInt(intAttribute1);
        if(WriteSelect(17, buf))
            options.Write(buf);
    }

    public void ReadAtts(int index, CommunicationBuffer buf)
    {
        switch(index)
        {
        case 0:
            SetObjectName(buf.ReadString());
            break;
        case 1:
            SetObjectType(buf.ReadInt());
            break;
        case 2:
            SetVisible(buf.ReadBool());
            break;
        case 3:
            SetActive(buf.ReadBool());
            break;
        case 4:
            SetPosition(buf.ReadDoubleArray());
            break;
        case 5:
            SetPosition2(buf.ReadDoubleArray());
            break;
        case 6:
            textColor.Read(buf);
            Select(6);
            break;
        case 7:
            SetUseForegroundForTextColor(buf.ReadBool());
            break;
        case 8:
            color1.Read(buf);
            Select(8);
            break;
        case 9:
            color2.Read(buf);
            Select(9);
            break;
        case 10:
            SetText(buf.ReadStringVector());
            break;
        case 11:
            SetFontFamily(buf.ReadInt());
            break;
        case 12:
            SetFontBold(buf.ReadBool());
            break;
        case 13:
            SetFontItalic(buf.ReadBool());
            break;
        case 14:
            SetFontShadow(buf.ReadBool());
            break;
        case 15:
            SetDoubleAttribute1(buf.ReadDouble());
            break;
        case 16:
            SetIntAttribute1(buf.ReadInt());
            break;
        case 17:
            options.Read(buf);
            break;
        }
    }

    public String toString(String indent)
    {
        String str = new String();
        str = str + stringToString("objectName", objectName, indent) + "\n";
        str = str + indent + "objectType = ";
        if(objectType == ANNOTATIONTYPE_TEXT2D)
            str = str + "ANNOTATIONTYPE_TEXT2D";
        if(objectType == ANNOTATIONTYPE_TEXT3D)
            str = str + "ANNOTATIONTYPE_TEXT3D";
        if(objectType == ANNOTATIONTYPE_TIMESLIDER)
            str = str + "ANNOTATIONTYPE_TIMESLIDER";
        if(objectType == ANNOTATIONTYPE_LINE2D)
            str = str + "ANNOTATIONTYPE_LINE2D";
        if(objectType == ANNOTATIONTYPE_LINE3D)
            str = str + "ANNOTATIONTYPE_LINE3D";
        if(objectType == ANNOTATIONTYPE_ARROW2D)
            str = str + "ANNOTATIONTYPE_ARROW2D";
        if(objectType == ANNOTATIONTYPE_ARROW3D)
            str = str + "ANNOTATIONTYPE_ARROW3D";
        if(objectType == ANNOTATIONTYPE_BOX)
            str = str + "ANNOTATIONTYPE_BOX";
        if(objectType == ANNOTATIONTYPE_IMAGE)
            str = str + "ANNOTATIONTYPE_IMAGE";
        if(objectType == ANNOTATIONTYPE_LEGENDATTRIBUTES)
            str = str + "ANNOTATIONTYPE_LEGENDATTRIBUTES";
        if(objectType == ANNOTATIONTYPE_MAXANNOTATIONTYPE)
            str = str + "ANNOTATIONTYPE_MAXANNOTATIONTYPE";
        str = str + "\n";
        str = str + boolToString("visible", visible, indent) + "\n";
        str = str + boolToString("active", active, indent) + "\n";
        str = str + doubleArrayToString("position", position, indent) + "\n";
        str = str + doubleArrayToString("position2", position2, indent) + "\n";
        str = str + indent + "textColor = {" + textColor.Red() + ", " + textColor.Green() + ", " + textColor.Blue() + ", " + textColor.Alpha() + "}\n";
        str = str + boolToString("useForegroundForTextColor", useForegroundForTextColor, indent) + "\n";
        str = str + indent + "color1 = {" + color1.Red() + ", " + color1.Green() + ", " + color1.Blue() + ", " + color1.Alpha() + "}\n";
        str = str + indent + "color2 = {" + color2.Red() + ", " + color2.Green() + ", " + color2.Blue() + ", " + color2.Alpha() + "}\n";
        str = str + stringVectorToString("text", text, indent) + "\n";
        str = str + indent + "fontFamily = ";
        if(fontFamily == FONTFAMILY_ARIAL)
            str = str + "FONTFAMILY_ARIAL";
        if(fontFamily == FONTFAMILY_COURIER)
            str = str + "FONTFAMILY_COURIER";
        if(fontFamily == FONTFAMILY_TIMES)
            str = str + "FONTFAMILY_TIMES";
        str = str + "\n";
        str = str + boolToString("fontBold", fontBold, indent) + "\n";
        str = str + boolToString("fontItalic", fontItalic, indent) + "\n";
        str = str + boolToString("fontShadow", fontShadow, indent) + "\n";
        str = str + doubleToString("doubleAttribute1", doubleAttribute1, indent) + "\n";
        str = str + intToString("intAttribute1", intAttribute1, indent) + "\n";
        str = str + indent + "options = " + options.toString(indent);
        return str;
    }


    // Attributes
    private String         objectName;
    private int            objectType;
    private boolean        visible;
    private boolean        active;
    private double[]       position;
    private double[]       position2;
    private ColorAttribute textColor;
    private boolean        useForegroundForTextColor;
    private ColorAttribute color1;
    private ColorAttribute color2;
    private Vector         text; // vector of String objects
    private int            fontFamily;
    private boolean        fontBold;
    private boolean        fontItalic;
    private boolean        fontShadow;
    private double         doubleAttribute1;
    private int            intAttribute1;
    private MapNode        options;
}

