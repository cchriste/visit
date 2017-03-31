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
// Class: AppearanceAttributes
//
// Purpose:
//    This class contains the GUI/viewer appearance attributes.
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

public class AppearanceAttributes extends AttributeSubject
{
    private static int AppearanceAttributes_numAdditionalAtts = 11;

    public AppearanceAttributes()
    {
        super(AppearanceAttributes_numAdditionalAtts);

        useSystemDefault = true;
        background = new String("#c0c0c0");
        foreground = new String("#000000");
        fontName = new String("Helvetica,12,-1,5,50,0,0,0,0,0");
        style = new String("motif");
        orientation = 0;
        defaultForeground = new String("");
        defaultBackground = new String("");
        defaultFontName = new String("");
        defaultStyle = new String("");
        defaultOrientation = 0;
    }

    public AppearanceAttributes(int nMoreFields)
    {
        super(AppearanceAttributes_numAdditionalAtts + nMoreFields);

        useSystemDefault = true;
        background = new String("#c0c0c0");
        foreground = new String("#000000");
        fontName = new String("Helvetica,12,-1,5,50,0,0,0,0,0");
        style = new String("motif");
        orientation = 0;
        defaultForeground = new String("");
        defaultBackground = new String("");
        defaultFontName = new String("");
        defaultStyle = new String("");
        defaultOrientation = 0;
    }

    public AppearanceAttributes(AppearanceAttributes obj)
    {
        super(obj);

        useSystemDefault = obj.useSystemDefault;
        background = new String(obj.background);
        foreground = new String(obj.foreground);
        fontName = new String(obj.fontName);
        style = new String(obj.style);
        orientation = obj.orientation;
        defaultForeground = new String(obj.defaultForeground);
        defaultBackground = new String(obj.defaultBackground);
        defaultFontName = new String(obj.defaultFontName);
        defaultStyle = new String(obj.defaultStyle);
        defaultOrientation = obj.defaultOrientation;

        SelectAll();
    }

    public int Offset()
    {
        return super.Offset() + super.GetNumAdditionalAttributes();
    }

    public int GetNumAdditionalAttributes()
    {
        return AppearanceAttributes_numAdditionalAtts;
    }

    public boolean equals(AppearanceAttributes obj)
    {
        // Create the return value
        return ((useSystemDefault == obj.useSystemDefault) &&
                (background.equals(obj.background)) &&
                (foreground.equals(obj.foreground)) &&
                (fontName.equals(obj.fontName)) &&
                (style.equals(obj.style)) &&
                (orientation == obj.orientation) &&
                (defaultForeground.equals(obj.defaultForeground)) &&
                (defaultBackground.equals(obj.defaultBackground)) &&
                (defaultFontName.equals(obj.defaultFontName)) &&
                (defaultStyle.equals(obj.defaultStyle)) &&
                (defaultOrientation == obj.defaultOrientation));
    }

    // Property setting methods
    public void SetUseSystemDefault(boolean useSystemDefault_)
    {
        useSystemDefault = useSystemDefault_;
        Select(0);
    }

    public void SetBackground(String background_)
    {
        background = background_;
        Select(1);
    }

    public void SetForeground(String foreground_)
    {
        foreground = foreground_;
        Select(2);
    }

    public void SetFontName(String fontName_)
    {
        fontName = fontName_;
        Select(3);
    }

    public void SetStyle(String style_)
    {
        style = style_;
        Select(4);
    }

    public void SetOrientation(int orientation_)
    {
        orientation = orientation_;
        Select(5);
    }

    public void SetDefaultForeground(String defaultForeground_)
    {
        defaultForeground = defaultForeground_;
        Select(6);
    }

    public void SetDefaultBackground(String defaultBackground_)
    {
        defaultBackground = defaultBackground_;
        Select(7);
    }

    public void SetDefaultFontName(String defaultFontName_)
    {
        defaultFontName = defaultFontName_;
        Select(8);
    }

    public void SetDefaultStyle(String defaultStyle_)
    {
        defaultStyle = defaultStyle_;
        Select(9);
    }

    public void SetDefaultOrientation(int defaultOrientation_)
    {
        defaultOrientation = defaultOrientation_;
        Select(10);
    }

    // Property getting methods
    public boolean GetUseSystemDefault() { return useSystemDefault; }
    public String  GetBackground() { return background; }
    public String  GetForeground() { return foreground; }
    public String  GetFontName() { return fontName; }
    public String  GetStyle() { return style; }
    public int     GetOrientation() { return orientation; }
    public String  GetDefaultForeground() { return defaultForeground; }
    public String  GetDefaultBackground() { return defaultBackground; }
    public String  GetDefaultFontName() { return defaultFontName; }
    public String  GetDefaultStyle() { return defaultStyle; }
    public int     GetDefaultOrientation() { return defaultOrientation; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteBool(useSystemDefault);
        if(WriteSelect(1, buf))
            buf.WriteString(background);
        if(WriteSelect(2, buf))
            buf.WriteString(foreground);
        if(WriteSelect(3, buf))
            buf.WriteString(fontName);
        if(WriteSelect(4, buf))
            buf.WriteString(style);
        if(WriteSelect(5, buf))
            buf.WriteInt(orientation);
        if(WriteSelect(6, buf))
            buf.WriteString(defaultForeground);
        if(WriteSelect(7, buf))
            buf.WriteString(defaultBackground);
        if(WriteSelect(8, buf))
            buf.WriteString(defaultFontName);
        if(WriteSelect(9, buf))
            buf.WriteString(defaultStyle);
        if(WriteSelect(10, buf))
            buf.WriteInt(defaultOrientation);
    }

    public void ReadAtts(int index, CommunicationBuffer buf)
    {
        switch(index)
        {
        case 0:
            SetUseSystemDefault(buf.ReadBool());
            break;
        case 1:
            SetBackground(buf.ReadString());
            break;
        case 2:
            SetForeground(buf.ReadString());
            break;
        case 3:
            SetFontName(buf.ReadString());
            break;
        case 4:
            SetStyle(buf.ReadString());
            break;
        case 5:
            SetOrientation(buf.ReadInt());
            break;
        case 6:
            SetDefaultForeground(buf.ReadString());
            break;
        case 7:
            SetDefaultBackground(buf.ReadString());
            break;
        case 8:
            SetDefaultFontName(buf.ReadString());
            break;
        case 9:
            SetDefaultStyle(buf.ReadString());
            break;
        case 10:
            SetDefaultOrientation(buf.ReadInt());
            break;
        }
    }

    public String toString(String indent)
    {
        String str = new String();
        str = str + boolToString("useSystemDefault", useSystemDefault, indent) + "\n";
        str = str + stringToString("background", background, indent) + "\n";
        str = str + stringToString("foreground", foreground, indent) + "\n";
        str = str + stringToString("fontName", fontName, indent) + "\n";
        str = str + stringToString("style", style, indent) + "\n";
        str = str + intToString("orientation", orientation, indent) + "\n";
        str = str + stringToString("defaultForeground", defaultForeground, indent) + "\n";
        str = str + stringToString("defaultBackground", defaultBackground, indent) + "\n";
        str = str + stringToString("defaultFontName", defaultFontName, indent) + "\n";
        str = str + stringToString("defaultStyle", defaultStyle, indent) + "\n";
        str = str + intToString("defaultOrientation", defaultOrientation, indent) + "\n";
        return str;
    }


    // Attributes
    private boolean useSystemDefault;
    private String  background;
    private String  foreground;
    private String  fontName;
    private String  style;
    private int     orientation;
    private String  defaultForeground;
    private String  defaultBackground;
    private String  defaultFontName;
    private String  defaultStyle;
    private int     defaultOrientation;
}

