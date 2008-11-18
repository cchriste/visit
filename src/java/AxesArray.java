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

package llnl.visit;


// ****************************************************************************
// Class: AxesArray
//
// Purpose:
//    Contains the properties for the array axes.
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

public class AxesArray extends AttributeSubject
{
    public AxesArray()
    {
        super(6);

        visible = true;
        ticksVisible = true;
        autoSetTicks = true;
        autoSetScaling = true;
        lineWidth = 0;
        axes = new AxisAttributes();
    }

    public AxesArray(AxesArray obj)
    {
        super(6);

        visible = obj.visible;
        ticksVisible = obj.ticksVisible;
        autoSetTicks = obj.autoSetTicks;
        autoSetScaling = obj.autoSetScaling;
        lineWidth = obj.lineWidth;
        axes = new AxisAttributes(obj.axes);

        SelectAll();
    }

    public boolean equals(AxesArray obj)
    {
        // Create the return value
        return ((visible == obj.visible) &&
                (ticksVisible == obj.ticksVisible) &&
                (autoSetTicks == obj.autoSetTicks) &&
                (autoSetScaling == obj.autoSetScaling) &&
                (lineWidth == obj.lineWidth) &&
                (axes.equals(obj.axes)));
    }

    // Property setting methods
    public void SetVisible(boolean visible_)
    {
        visible = visible_;
        Select(0);
    }

    public void SetTicksVisible(boolean ticksVisible_)
    {
        ticksVisible = ticksVisible_;
        Select(1);
    }

    public void SetAutoSetTicks(boolean autoSetTicks_)
    {
        autoSetTicks = autoSetTicks_;
        Select(2);
    }

    public void SetAutoSetScaling(boolean autoSetScaling_)
    {
        autoSetScaling = autoSetScaling_;
        Select(3);
    }

    public void SetLineWidth(int lineWidth_)
    {
        lineWidth = lineWidth_;
        Select(4);
    }

    public void SetAxes(AxisAttributes axes_)
    {
        axes = axes_;
        Select(5);
    }

    // Property getting methods
    public boolean        GetVisible() { return visible; }
    public boolean        GetTicksVisible() { return ticksVisible; }
    public boolean        GetAutoSetTicks() { return autoSetTicks; }
    public boolean        GetAutoSetScaling() { return autoSetScaling; }
    public int            GetLineWidth() { return lineWidth; }
    public AxisAttributes GetAxes() { return axes; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteBool(visible);
        if(WriteSelect(1, buf))
            buf.WriteBool(ticksVisible);
        if(WriteSelect(2, buf))
            buf.WriteBool(autoSetTicks);
        if(WriteSelect(3, buf))
            buf.WriteBool(autoSetScaling);
        if(WriteSelect(4, buf))
            buf.WriteInt(lineWidth);
        if(WriteSelect(5, buf))
            axes.Write(buf);
    }

    public void ReadAtts(int n, CommunicationBuffer buf)
    {
        for(int i = 0; i < n; ++i)
        {
            int index = (int)buf.ReadByte();
            switch(index)
            {
            case 0:
                SetVisible(buf.ReadBool());
                break;
            case 1:
                SetTicksVisible(buf.ReadBool());
                break;
            case 2:
                SetAutoSetTicks(buf.ReadBool());
                break;
            case 3:
                SetAutoSetScaling(buf.ReadBool());
                break;
            case 4:
                SetLineWidth(buf.ReadInt());
                break;
            case 5:
                axes.Read(buf);
                Select(5);
                break;
            }
        }
    }

    public String toString(String indent)
    {
        String str = new String();
        str = str + boolToString("visible", visible, indent) + "\n";
        str = str + boolToString("ticksVisible", ticksVisible, indent) + "\n";
        str = str + boolToString("autoSetTicks", autoSetTicks, indent) + "\n";
        str = str + boolToString("autoSetScaling", autoSetScaling, indent) + "\n";
        str = str + intToString("lineWidth", lineWidth, indent) + "\n";
        str = str + indent + "axes = {\n" + axes.toString(indent + "    ") + indent + "}\n";
        return str;
    }


    // Attributes
    private boolean        visible;
    private boolean        ticksVisible;
    private boolean        autoSetTicks;
    private boolean        autoSetScaling;
    private int            lineWidth;
    private AxisAttributes axes;
}

