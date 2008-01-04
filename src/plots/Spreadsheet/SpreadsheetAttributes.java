// ***************************************************************************
//
// Copyright (c) 2000 - 2008, Lawrence Livermore National Security, LLC
// Produced at the Lawrence Livermore National Laboratory
// LLNL-CODE-400142
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

package llnl.visit.plots;

import llnl.visit.AttributeSubject;
import llnl.visit.CommunicationBuffer;
import llnl.visit.Plugin;
import llnl.visit.ColorAttribute;
import java.lang.Double;
import java.util.Vector;

// ****************************************************************************
// Class: SpreadsheetAttributes
//
// Purpose:
//    Contains the attributes for the visual spreadsheet.
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   Wed Nov 7 14:49:57 PST 2007
//
// Modifications:
//   
// ****************************************************************************

public class SpreadsheetAttributes extends AttributeSubject implements Plugin
{
    // Enum values
    public final static int NORMALAXIS_X = 0;
    public final static int NORMALAXIS_Y = 1;
    public final static int NORMALAXIS_Z = 2;


    public SpreadsheetAttributes()
    {
        super(15);

        subsetName = new String("Whole");
        formatString = new String("%1.6f");
        useColorTable = false;
        colorTableName = new String("Default");
        showTracerPlane = true;
        tracerColor = new ColorAttribute(255, 0, 0, 150);
        normal = NORMALAXIS_Z;
        sliceIndex = 0;
        currentPick = new double[3];
        for (int i = 0; i < currentPick.length; ++i)
            currentPick[i] = 0.;
        currentPickValid = false;
        pastPicks = new Vector();
        currentPickLetter = new String("");
        pastPickLetters = new Vector();
        spreadsheetFont = new String("Courier,12,-1,5,50,0,0,0,0,0");
        showPatchOutline = true;
    }

    public SpreadsheetAttributes(SpreadsheetAttributes obj)
    {
        super(15);

        int i;

        subsetName = new String(obj.subsetName);
        formatString = new String(obj.formatString);
        useColorTable = obj.useColorTable;
        colorTableName = new String(obj.colorTableName);
        showTracerPlane = obj.showTracerPlane;
        tracerColor = new ColorAttribute(obj.tracerColor);
        normal = obj.normal;
        sliceIndex = obj.sliceIndex;
        currentPick = new double[3];
        currentPick[0] = obj.currentPick[0];
        currentPick[1] = obj.currentPick[1];
        currentPick[2] = obj.currentPick[2];

        currentPickValid = obj.currentPickValid;
        pastPicks = new Vector(obj.pastPicks.size());
        for(i = 0; i < obj.pastPicks.size(); ++i)
        {
            Double dv = (Double)obj.pastPicks.elementAt(i);
            pastPicks.addElement(new Double(dv.doubleValue()));
        }

        currentPickLetter = new String(obj.currentPickLetter);
        pastPickLetters = new Vector(obj.pastPickLetters.size());
        for(i = 0; i < obj.pastPickLetters.size(); ++i)
            pastPickLetters.addElement(new String((String)obj.pastPickLetters.elementAt(i)));

        spreadsheetFont = new String(obj.spreadsheetFont);
        showPatchOutline = obj.showPatchOutline;

        SelectAll();
    }

    public boolean equals(SpreadsheetAttributes obj)
    {
        int i;

        // Compare the currentPick arrays.
        boolean currentPick_equal = true;
        for(i = 0; i < 3 && currentPick_equal; ++i)
            currentPick_equal = (currentPick[i] == obj.currentPick[i]);

        // Create the return value
        return ((subsetName == obj.subsetName) &&
                (formatString == obj.formatString) &&
                (useColorTable == obj.useColorTable) &&
                (colorTableName == obj.colorTableName) &&
                (showTracerPlane == obj.showTracerPlane) &&
                (tracerColor == obj.tracerColor) &&
                (normal == obj.normal) &&
                (sliceIndex == obj.sliceIndex) &&
                currentPick_equal &&
                (currentPickValid == obj.currentPickValid) &&
                (pastPicks == obj.pastPicks) &&
                (currentPickLetter == obj.currentPickLetter) &&
                (pastPickLetters == obj.pastPickLetters) &&
                (spreadsheetFont == obj.spreadsheetFont) &&
                (showPatchOutline == obj.showPatchOutline));
    }

    public String GetName() { return "Spreadsheet"; }
    public String GetVersion() { return "1.0"; }

    // Property setting methods
    public void SetSubsetName(String subsetName_)
    {
        subsetName = subsetName_;
        Select(0);
    }

    public void SetFormatString(String formatString_)
    {
        formatString = formatString_;
        Select(1);
    }

    public void SetUseColorTable(boolean useColorTable_)
    {
        useColorTable = useColorTable_;
        Select(2);
    }

    public void SetColorTableName(String colorTableName_)
    {
        colorTableName = colorTableName_;
        Select(3);
    }

    public void SetShowTracerPlane(boolean showTracerPlane_)
    {
        showTracerPlane = showTracerPlane_;
        Select(4);
    }

    public void SetTracerColor(ColorAttribute tracerColor_)
    {
        tracerColor = tracerColor_;
        Select(5);
    }

    public void SetNormal(int normal_)
    {
        normal = normal_;
        Select(6);
    }

    public void SetSliceIndex(int sliceIndex_)
    {
        sliceIndex = sliceIndex_;
        Select(7);
    }

    public void SetCurrentPick(double[] currentPick_)
    {
        currentPick[0] = currentPick_[0];
        currentPick[1] = currentPick_[1];
        currentPick[2] = currentPick_[2];
        Select(8);
    }

    public void SetCurrentPick(double e0, double e1, double e2)
    {
        currentPick[0] = e0;
        currentPick[1] = e1;
        currentPick[2] = e2;
        Select(8);
    }

    public void SetCurrentPickValid(boolean currentPickValid_)
    {
        currentPickValid = currentPickValid_;
        Select(9);
    }

    public void SetPastPicks(Vector pastPicks_)
    {
        pastPicks = pastPicks_;
        Select(10);
    }

    public void SetCurrentPickLetter(String currentPickLetter_)
    {
        currentPickLetter = currentPickLetter_;
        Select(11);
    }

    public void SetPastPickLetters(Vector pastPickLetters_)
    {
        pastPickLetters = pastPickLetters_;
        Select(12);
    }

    public void SetSpreadsheetFont(String spreadsheetFont_)
    {
        spreadsheetFont = spreadsheetFont_;
        Select(13);
    }

    public void SetShowPatchOutline(boolean showPatchOutline_)
    {
        showPatchOutline = showPatchOutline_;
        Select(14);
    }

    // Property getting methods
    public String         GetSubsetName() { return subsetName; }
    public String         GetFormatString() { return formatString; }
    public boolean        GetUseColorTable() { return useColorTable; }
    public String         GetColorTableName() { return colorTableName; }
    public boolean        GetShowTracerPlane() { return showTracerPlane; }
    public ColorAttribute GetTracerColor() { return tracerColor; }
    public int            GetNormal() { return normal; }
    public int            GetSliceIndex() { return sliceIndex; }
    public double[]       GetCurrentPick() { return currentPick; }
    public boolean        GetCurrentPickValid() { return currentPickValid; }
    public Vector         GetPastPicks() { return pastPicks; }
    public String         GetCurrentPickLetter() { return currentPickLetter; }
    public Vector         GetPastPickLetters() { return pastPickLetters; }
    public String         GetSpreadsheetFont() { return spreadsheetFont; }
    public boolean        GetShowPatchOutline() { return showPatchOutline; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteString(subsetName);
        if(WriteSelect(1, buf))
            buf.WriteString(formatString);
        if(WriteSelect(2, buf))
            buf.WriteBool(useColorTable);
        if(WriteSelect(3, buf))
            buf.WriteString(colorTableName);
        if(WriteSelect(4, buf))
            buf.WriteBool(showTracerPlane);
        if(WriteSelect(5, buf))
            tracerColor.Write(buf);
        if(WriteSelect(6, buf))
            buf.WriteInt(normal);
        if(WriteSelect(7, buf))
            buf.WriteInt(sliceIndex);
        if(WriteSelect(8, buf))
            buf.WriteDoubleArray(currentPick);
        if(WriteSelect(9, buf))
            buf.WriteBool(currentPickValid);
        if(WriteSelect(10, buf))
            buf.WriteDoubleVector(pastPicks);
        if(WriteSelect(11, buf))
            buf.WriteString(currentPickLetter);
        if(WriteSelect(12, buf))
            buf.WriteStringVector(pastPickLetters);
        if(WriteSelect(13, buf))
            buf.WriteString(spreadsheetFont);
        if(WriteSelect(14, buf))
            buf.WriteBool(showPatchOutline);
    }

    public void ReadAtts(int n, CommunicationBuffer buf)
    {
        for(int i = 0; i < n; ++i)
        {
            int index = (int)buf.ReadByte();
            switch(index)
            {
            case 0:
                SetSubsetName(buf.ReadString());
                break;
            case 1:
                SetFormatString(buf.ReadString());
                break;
            case 2:
                SetUseColorTable(buf.ReadBool());
                break;
            case 3:
                SetColorTableName(buf.ReadString());
                break;
            case 4:
                SetShowTracerPlane(buf.ReadBool());
                break;
            case 5:
                tracerColor.Read(buf);
                Select(5);
                break;
            case 6:
                SetNormal(buf.ReadInt());
                break;
            case 7:
                SetSliceIndex(buf.ReadInt());
                break;
            case 8:
                SetCurrentPick(buf.ReadDoubleArray());
                break;
            case 9:
                SetCurrentPickValid(buf.ReadBool());
                break;
            case 10:
                SetPastPicks(buf.ReadDoubleVector());
                break;
            case 11:
                SetCurrentPickLetter(buf.ReadString());
                break;
            case 12:
                SetPastPickLetters(buf.ReadStringVector());
                break;
            case 13:
                SetSpreadsheetFont(buf.ReadString());
                break;
            case 14:
                SetShowPatchOutline(buf.ReadBool());
                break;
            }
        }
    }


    // Attributes
    private String         subsetName;
    private String         formatString;
    private boolean        useColorTable;
    private String         colorTableName;
    private boolean        showTracerPlane;
    private ColorAttribute tracerColor;
    private int            normal;
    private int            sliceIndex;
    private double[]       currentPick;
    private boolean        currentPickValid;
    private Vector         pastPicks; // vector of Double objects
    private String         currentPickLetter;
    private Vector         pastPickLetters; // vector of String objects
    private String         spreadsheetFont;
    private boolean        showPatchOutline;
}

