// ***************************************************************************
//
// Copyright (c) 2000 - 2010, Lawrence Livermore National Security, LLC
// Produced at the Lawrence Livermore National Laboratory
// LLNL-CODE-400124
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
// Class: AxisTitles
//
// Purpose:
//    Contains the title properties for one axis.
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

public class AxisTitles extends AttributeSubject
{
    private static int AxisTitles_numAdditionalAtts = 6;

    public AxisTitles()
    {
        super(AxisTitles_numAdditionalAtts);

        visible = true;
        font = new FontAttributes();
        userTitle = false;
        userUnits = false;
        title = new String("");
        units = new String("");
    }

    public AxisTitles(int nMoreFields)
    {
        super(AxisTitles_numAdditionalAtts + nMoreFields);

        visible = true;
        font = new FontAttributes();
        userTitle = false;
        userUnits = false;
        title = new String("");
        units = new String("");
    }

    public AxisTitles(AxisTitles obj)
    {
        super(AxisTitles_numAdditionalAtts);

        visible = obj.visible;
        font = new FontAttributes(obj.font);
        userTitle = obj.userTitle;
        userUnits = obj.userUnits;
        title = new String(obj.title);
        units = new String(obj.units);

        SelectAll();
    }

    public int Offset()
    {
        return super.Offset() + super.GetNumAdditionalAttributes();
    }

    public int GetNumAdditionalAttributes()
    {
        return AxisTitles_numAdditionalAtts;
    }

    public boolean equals(AxisTitles obj)
    {
        // Create the return value
        return ((visible == obj.visible) &&
                (font.equals(obj.font)) &&
                (userTitle == obj.userTitle) &&
                (userUnits == obj.userUnits) &&
                (title.equals(obj.title)) &&
                (units.equals(obj.units)));
    }

    // Property setting methods
    public void SetVisible(boolean visible_)
    {
        visible = visible_;
        Select(0);
    }

    public void SetFont(FontAttributes font_)
    {
        font = font_;
        Select(1);
    }

    public void SetUserTitle(boolean userTitle_)
    {
        userTitle = userTitle_;
        Select(2);
    }

    public void SetUserUnits(boolean userUnits_)
    {
        userUnits = userUnits_;
        Select(3);
    }

    public void SetTitle(String title_)
    {
        title = title_;
        Select(4);
    }

    public void SetUnits(String units_)
    {
        units = units_;
        Select(5);
    }

    // Property getting methods
    public boolean        GetVisible() { return visible; }
    public FontAttributes GetFont() { return font; }
    public boolean        GetUserTitle() { return userTitle; }
    public boolean        GetUserUnits() { return userUnits; }
    public String         GetTitle() { return title; }
    public String         GetUnits() { return units; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteBool(visible);
        if(WriteSelect(1, buf))
            font.Write(buf);
        if(WriteSelect(2, buf))
            buf.WriteBool(userTitle);
        if(WriteSelect(3, buf))
            buf.WriteBool(userUnits);
        if(WriteSelect(4, buf))
            buf.WriteString(title);
        if(WriteSelect(5, buf))
            buf.WriteString(units);
    }

    public void ReadAtts(int index, CommunicationBuffer buf)
    {
        switch(index)
        {
        case 0:
            SetVisible(buf.ReadBool());
            break;
        case 1:
            font.Read(buf);
            Select(1);
            break;
        case 2:
            SetUserTitle(buf.ReadBool());
            break;
        case 3:
            SetUserUnits(buf.ReadBool());
            break;
        case 4:
            SetTitle(buf.ReadString());
            break;
        case 5:
            SetUnits(buf.ReadString());
            break;
        }
    }

    public String toString(String indent)
    {
        String str = new String();
        str = str + boolToString("visible", visible, indent) + "\n";
        str = str + indent + "font = {\n" + font.toString(indent + "    ") + indent + "}\n";
        str = str + boolToString("userTitle", userTitle, indent) + "\n";
        str = str + boolToString("userUnits", userUnits, indent) + "\n";
        str = str + stringToString("title", title, indent) + "\n";
        str = str + stringToString("units", units, indent) + "\n";
        return str;
    }


    // Attributes
    private boolean        visible;
    private FontAttributes font;
    private boolean        userTitle;
    private boolean        userUnits;
    private String         title;
    private String         units;
}

