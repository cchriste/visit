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

import java.util.Vector;

// ****************************************************************************
// Class: avtDefaultPlotMetaData
//
// Purpose:
//    Contains default plot metadata attributes
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   Wed Mar 14 17:56:04 PST 2007
//
// Modifications:
//   
// ****************************************************************************

public class avtDefaultPlotMetaData extends AttributeSubject
{
    public avtDefaultPlotMetaData()
    {
        super(3);

        pluginID = new String("");
        plotVar = new String("var");
        plotAttributes = new Vector();
    }

    public avtDefaultPlotMetaData(avtDefaultPlotMetaData obj)
    {
        super(3);

        int i;

        pluginID = new String(obj.pluginID);
        plotVar = new String(obj.plotVar);
        plotAttributes = new Vector(obj.plotAttributes.size());
        for(i = 0; i < obj.plotAttributes.size(); ++i)
            plotAttributes.addElement(new String((String)obj.plotAttributes.elementAt(i)));


        SelectAll();
    }

    public boolean equals(avtDefaultPlotMetaData obj)
    {
        int i;

        // Create the return value
        return ((pluginID == obj.pluginID) &&
                (plotVar == obj.plotVar) &&
                (plotAttributes == obj.plotAttributes));
    }

    // Property setting methods
    public void SetPluginID(String pluginID_)
    {
        pluginID = pluginID_;
        Select(0);
    }

    public void SetPlotVar(String plotVar_)
    {
        plotVar = plotVar_;
        Select(1);
    }

    public void SetPlotAttributes(Vector plotAttributes_)
    {
        plotAttributes = plotAttributes_;
        Select(2);
    }

    // Property getting methods
    public String GetPluginID() { return pluginID; }
    public String GetPlotVar() { return plotVar; }
    public Vector GetPlotAttributes() { return plotAttributes; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteString(pluginID);
        if(WriteSelect(1, buf))
            buf.WriteString(plotVar);
        if(WriteSelect(2, buf))
            buf.WriteStringVector(plotAttributes);
    }

    public void ReadAtts(int n, CommunicationBuffer buf)
    {
        for(int i = 0; i < n; ++i)
        {
            int index = (int)buf.ReadByte();
            switch(index)
            {
            case 0:
                SetPluginID(buf.ReadString());
                break;
            case 1:
                SetPlotVar(buf.ReadString());
                break;
            case 2:
                SetPlotAttributes(buf.ReadStringVector());
                break;
            }
        }
    }


    // Attributes
    private String pluginID;
    private String plotVar;
    private Vector plotAttributes; // vector of String objects
}

