// ***************************************************************************
//
// Copyright (c) 2000 - 2010, Lawrence Livermore National Security, LLC
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

package llnl.visit.operators;

import llnl.visit.AttributeSubject;
import llnl.visit.CommunicationBuffer;
import llnl.visit.Plugin;

// ****************************************************************************
// Class: MultiresControlAttributes
//
// Purpose:
//    
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

public class MultiresControlAttributes extends AttributeSubject implements Plugin
{
    private static int MultiresControlAttributes_numAdditionalAtts = 3;

    public MultiresControlAttributes()
    {
        super(MultiresControlAttributes_numAdditionalAtts);

        resolution = 0;
        maxResolution = 1;
        info = new String("");
    }

    public MultiresControlAttributes(int nMoreFields)
    {
        super(MultiresControlAttributes_numAdditionalAtts + nMoreFields);

        resolution = 0;
        maxResolution = 1;
        info = new String("");
    }

    public MultiresControlAttributes(MultiresControlAttributes obj)
    {
        super(MultiresControlAttributes_numAdditionalAtts);

        resolution = obj.resolution;
        maxResolution = obj.maxResolution;
        info = new String(obj.info);

        SelectAll();
    }

    public int Offset()
    {
        return super.Offset() + super.GetNumAdditionalAttributes();
    }

    public int GetNumAdditionalAttributes()
    {
        return MultiresControlAttributes_numAdditionalAtts;
    }

    public boolean equals(MultiresControlAttributes obj)
    {
        // Create the return value
        return ((resolution == obj.resolution) &&
                (maxResolution == obj.maxResolution) &&
                (info.equals(obj.info)));
    }

    public String GetName() { return "MultiresControl"; }
    public String GetVersion() { return "1.0"; }

    // Property setting methods
    public void SetResolution(int resolution_)
    {
        resolution = resolution_;
        Select(0);
    }

    public void SetMaxResolution(int maxResolution_)
    {
        maxResolution = maxResolution_;
        Select(1);
    }

    public void SetInfo(String info_)
    {
        info = info_;
        Select(2);
    }

    // Property getting methods
    public int    GetResolution() { return resolution; }
    public int    GetMaxResolution() { return maxResolution; }
    public String GetInfo() { return info; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteInt(resolution);
        if(WriteSelect(1, buf))
            buf.WriteInt(maxResolution);
        if(WriteSelect(2, buf))
            buf.WriteString(info);
    }

    public void ReadAtts(int index, CommunicationBuffer buf)
    {
        switch(index)
        {
        case 0:
            SetResolution(buf.ReadInt());
            break;
        case 1:
            SetMaxResolution(buf.ReadInt());
            break;
        case 2:
            SetInfo(buf.ReadString());
            break;
        }
    }

    public String toString(String indent)
    {
        String str = new String();
        str = str + intToString("resolution", resolution, indent) + "\n";
        str = str + intToString("maxResolution", maxResolution, indent) + "\n";
        str = str + stringToString("info", info, indent) + "\n";
        return str;
    }


    // Attributes
    private int    resolution;
    private int    maxResolution;
    private String info;
}

