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

// ****************************************************************************
// Class: LabelAttributes
//
// Purpose:
//    This class contains the fields that we need to set the attributes for the Label plot.
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   Tue Sep 18 16:48:38 PST 2007
//
// Modifications:
//   
// ****************************************************************************

public class LabelAttributes extends AttributeSubject implements Plugin
{
    // Enum values
    public final static int LABELINDEXDISPLAY_NATURAL = 0;
    public final static int LABELINDEXDISPLAY_LOGICALINDEX = 1;
    public final static int LABELINDEXDISPLAY_INDEX = 2;

    public final static int LABELHORIZONTALALIGNMENT_HCENTER = 0;
    public final static int LABELHORIZONTALALIGNMENT_LEFT = 1;
    public final static int LABELHORIZONTALALIGNMENT_RIGHT = 2;

    public final static int LABELVERTICALALIGNMENT_VCENTER = 0;
    public final static int LABELVERTICALALIGNMENT_TOP = 1;
    public final static int LABELVERTICALALIGNMENT_BOTTOM = 2;

    public final static int LABELDRAWFACING_FRONT = 0;
    public final static int LABELDRAWFACING_BACK = 1;
    public final static int LABELDRAWFACING_FRONTANDBACK = 2;

    public final static int VARIABLETYPE_LABEL_VT_MESH = 0;
    public final static int VARIABLETYPE_LABEL_VT_SCALAR_VAR = 1;
    public final static int VARIABLETYPE_LABEL_VT_VECTOR_VAR = 2;
    public final static int VARIABLETYPE_LABEL_VT_TENSOR_VAR = 3;
    public final static int VARIABLETYPE_LABEL_VT_SYMMETRIC_TENSOR_VAR = 4;
    public final static int VARIABLETYPE_LABEL_VT_ARRAY_VAR = 5;
    public final static int VARIABLETYPE_LABEL_VT_LABEL_VAR = 6;
    public final static int VARIABLETYPE_LABEL_VT_MATERIAL = 7;
    public final static int VARIABLETYPE_LABEL_VT_SUBSET = 8;
    public final static int VARIABLETYPE_LABEL_VT_UNKNOWN_TYPE = 9;

    public final static int DEPTHTESTMODE_LABEL_DT_AUTO = 0;
    public final static int DEPTHTESTMODE_LABEL_DT_ALWAYS = 1;
    public final static int DEPTHTESTMODE_LABEL_DT_NEVER = 2;


    public LabelAttributes()
    {
        super(18);

        varType = VARIABLETYPE_LABEL_VT_UNKNOWN_TYPE;
        legendFlag = true;
        showNodes = false;
        showCells = true;
        restrictNumberOfLabels = true;
        drawLabelsFacing = LABELDRAWFACING_FRONT;
        labelDisplayFormat = LABELINDEXDISPLAY_NATURAL;
        numberOfLabels = 200;
        specifyTextColor1 = false;
        textColor1 = new ColorAttribute(255, 0, 0, 0);
        textHeight1 = 0.02f;
        specifyTextColor2 = false;
        textColor2 = new ColorAttribute(0, 0, 255, 0);
        textHeight2 = 0.02f;
        horizontalJustification = LABELHORIZONTALALIGNMENT_HCENTER;
        verticalJustification = LABELVERTICALALIGNMENT_VCENTER;
        depthTestMode = DEPTHTESTMODE_LABEL_DT_AUTO;
        formatTemplate = new String("%g");
    }

    public LabelAttributes(LabelAttributes obj)
    {
        super(18);

        varType = obj.varType;
        legendFlag = obj.legendFlag;
        showNodes = obj.showNodes;
        showCells = obj.showCells;
        restrictNumberOfLabels = obj.restrictNumberOfLabels;
        drawLabelsFacing = obj.drawLabelsFacing;
        labelDisplayFormat = obj.labelDisplayFormat;
        numberOfLabels = obj.numberOfLabels;
        specifyTextColor1 = obj.specifyTextColor1;
        textColor1 = new ColorAttribute(obj.textColor1);
        textHeight1 = obj.textHeight1;
        specifyTextColor2 = obj.specifyTextColor2;
        textColor2 = new ColorAttribute(obj.textColor2);
        textHeight2 = obj.textHeight2;
        horizontalJustification = obj.horizontalJustification;
        verticalJustification = obj.verticalJustification;
        depthTestMode = obj.depthTestMode;
        formatTemplate = new String(obj.formatTemplate);

        SelectAll();
    }

    public boolean equals(LabelAttributes obj)
    {
        // Create the return value
        return (true /* can ignore varType */ &&
                (legendFlag == obj.legendFlag) &&
                (showNodes == obj.showNodes) &&
                (showCells == obj.showCells) &&
                (restrictNumberOfLabels == obj.restrictNumberOfLabels) &&
                (drawLabelsFacing == obj.drawLabelsFacing) &&
                (labelDisplayFormat == obj.labelDisplayFormat) &&
                (numberOfLabels == obj.numberOfLabels) &&
                (specifyTextColor1 == obj.specifyTextColor1) &&
                (textColor1 == obj.textColor1) &&
                (textHeight1 == obj.textHeight1) &&
                (specifyTextColor2 == obj.specifyTextColor2) &&
                (textColor2 == obj.textColor2) &&
                (textHeight2 == obj.textHeight2) &&
                (horizontalJustification == obj.horizontalJustification) &&
                (verticalJustification == obj.verticalJustification) &&
                (depthTestMode == obj.depthTestMode) &&
                (formatTemplate == obj.formatTemplate));
    }

    public String GetName() { return "Label"; }
    public String GetVersion() { return "1.0"; }

    // Property setting methods
    public void SetVarType(int varType_)
    {
        varType = varType_;
        Select(0);
    }

    public void SetLegendFlag(boolean legendFlag_)
    {
        legendFlag = legendFlag_;
        Select(1);
    }

    public void SetShowNodes(boolean showNodes_)
    {
        showNodes = showNodes_;
        Select(2);
    }

    public void SetShowCells(boolean showCells_)
    {
        showCells = showCells_;
        Select(3);
    }

    public void SetRestrictNumberOfLabels(boolean restrictNumberOfLabels_)
    {
        restrictNumberOfLabels = restrictNumberOfLabels_;
        Select(4);
    }

    public void SetDrawLabelsFacing(int drawLabelsFacing_)
    {
        drawLabelsFacing = drawLabelsFacing_;
        Select(5);
    }

    public void SetLabelDisplayFormat(int labelDisplayFormat_)
    {
        labelDisplayFormat = labelDisplayFormat_;
        Select(6);
    }

    public void SetNumberOfLabels(int numberOfLabels_)
    {
        numberOfLabels = numberOfLabels_;
        Select(7);
    }

    public void SetSpecifyTextColor1(boolean specifyTextColor1_)
    {
        specifyTextColor1 = specifyTextColor1_;
        Select(8);
    }

    public void SetTextColor1(ColorAttribute textColor1_)
    {
        textColor1 = textColor1_;
        Select(9);
    }

    public void SetTextHeight1(float textHeight1_)
    {
        textHeight1 = textHeight1_;
        Select(10);
    }

    public void SetSpecifyTextColor2(boolean specifyTextColor2_)
    {
        specifyTextColor2 = specifyTextColor2_;
        Select(11);
    }

    public void SetTextColor2(ColorAttribute textColor2_)
    {
        textColor2 = textColor2_;
        Select(12);
    }

    public void SetTextHeight2(float textHeight2_)
    {
        textHeight2 = textHeight2_;
        Select(13);
    }

    public void SetHorizontalJustification(int horizontalJustification_)
    {
        horizontalJustification = horizontalJustification_;
        Select(14);
    }

    public void SetVerticalJustification(int verticalJustification_)
    {
        verticalJustification = verticalJustification_;
        Select(15);
    }

    public void SetDepthTestMode(int depthTestMode_)
    {
        depthTestMode = depthTestMode_;
        Select(16);
    }

    public void SetFormatTemplate(String formatTemplate_)
    {
        formatTemplate = formatTemplate_;
        Select(17);
    }

    // Property getting methods
    public int            GetVarType() { return varType; }
    public boolean        GetLegendFlag() { return legendFlag; }
    public boolean        GetShowNodes() { return showNodes; }
    public boolean        GetShowCells() { return showCells; }
    public boolean        GetRestrictNumberOfLabels() { return restrictNumberOfLabels; }
    public int            GetDrawLabelsFacing() { return drawLabelsFacing; }
    public int            GetLabelDisplayFormat() { return labelDisplayFormat; }
    public int            GetNumberOfLabels() { return numberOfLabels; }
    public boolean        GetSpecifyTextColor1() { return specifyTextColor1; }
    public ColorAttribute GetTextColor1() { return textColor1; }
    public float          GetTextHeight1() { return textHeight1; }
    public boolean        GetSpecifyTextColor2() { return specifyTextColor2; }
    public ColorAttribute GetTextColor2() { return textColor2; }
    public float          GetTextHeight2() { return textHeight2; }
    public int            GetHorizontalJustification() { return horizontalJustification; }
    public int            GetVerticalJustification() { return verticalJustification; }
    public int            GetDepthTestMode() { return depthTestMode; }
    public String         GetFormatTemplate() { return formatTemplate; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteInt(varType);
        if(WriteSelect(1, buf))
            buf.WriteBool(legendFlag);
        if(WriteSelect(2, buf))
            buf.WriteBool(showNodes);
        if(WriteSelect(3, buf))
            buf.WriteBool(showCells);
        if(WriteSelect(4, buf))
            buf.WriteBool(restrictNumberOfLabels);
        if(WriteSelect(5, buf))
            buf.WriteInt(drawLabelsFacing);
        if(WriteSelect(6, buf))
            buf.WriteInt(labelDisplayFormat);
        if(WriteSelect(7, buf))
            buf.WriteInt(numberOfLabels);
        if(WriteSelect(8, buf))
            buf.WriteBool(specifyTextColor1);
        if(WriteSelect(9, buf))
            textColor1.Write(buf);
        if(WriteSelect(10, buf))
            buf.WriteFloat(textHeight1);
        if(WriteSelect(11, buf))
            buf.WriteBool(specifyTextColor2);
        if(WriteSelect(12, buf))
            textColor2.Write(buf);
        if(WriteSelect(13, buf))
            buf.WriteFloat(textHeight2);
        if(WriteSelect(14, buf))
            buf.WriteInt(horizontalJustification);
        if(WriteSelect(15, buf))
            buf.WriteInt(verticalJustification);
        if(WriteSelect(16, buf))
            buf.WriteInt(depthTestMode);
        if(WriteSelect(17, buf))
            buf.WriteString(formatTemplate);
    }

    public void ReadAtts(int n, CommunicationBuffer buf)
    {
        for(int i = 0; i < n; ++i)
        {
            int index = (int)buf.ReadByte();
            switch(index)
            {
            case 0:
                SetVarType(buf.ReadInt());
                break;
            case 1:
                SetLegendFlag(buf.ReadBool());
                break;
            case 2:
                SetShowNodes(buf.ReadBool());
                break;
            case 3:
                SetShowCells(buf.ReadBool());
                break;
            case 4:
                SetRestrictNumberOfLabels(buf.ReadBool());
                break;
            case 5:
                SetDrawLabelsFacing(buf.ReadInt());
                break;
            case 6:
                SetLabelDisplayFormat(buf.ReadInt());
                break;
            case 7:
                SetNumberOfLabels(buf.ReadInt());
                break;
            case 8:
                SetSpecifyTextColor1(buf.ReadBool());
                break;
            case 9:
                textColor1.Read(buf);
                Select(9);
                break;
            case 10:
                SetTextHeight1(buf.ReadFloat());
                break;
            case 11:
                SetSpecifyTextColor2(buf.ReadBool());
                break;
            case 12:
                textColor2.Read(buf);
                Select(12);
                break;
            case 13:
                SetTextHeight2(buf.ReadFloat());
                break;
            case 14:
                SetHorizontalJustification(buf.ReadInt());
                break;
            case 15:
                SetVerticalJustification(buf.ReadInt());
                break;
            case 16:
                SetDepthTestMode(buf.ReadInt());
                break;
            case 17:
                SetFormatTemplate(buf.ReadString());
                break;
            }
        }
    }


    // Attributes
    private int            varType;
    private boolean        legendFlag;
    private boolean        showNodes;
    private boolean        showCells;
    private boolean        restrictNumberOfLabels;
    private int            drawLabelsFacing;
    private int            labelDisplayFormat;
    private int            numberOfLabels;
    private boolean        specifyTextColor1;
    private ColorAttribute textColor1;
    private float          textHeight1;
    private boolean        specifyTextColor2;
    private ColorAttribute textColor2;
    private float          textHeight2;
    private int            horizontalJustification;
    private int            verticalJustification;
    private int            depthTestMode;
    private String         formatTemplate;
}

