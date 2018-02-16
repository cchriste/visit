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


// ****************************************************************************
// Class: PlotQueryInfo
//
// Purpose:
//    This class is a .
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

public class PlotQueryInfo extends AttributeSubject
{
    private static int PlotQueryInfo_numAdditionalAtts = 3;

    // Enum values
    public final static int CHANGETYPE_NONE = 0;
    public final static int CHANGETYPE_DATABASE = 1;
    public final static int CHANGETYPE_VARNAME = 2;
    public final static int CHANGETYPE_ADDOP = 3;
    public final static int CHANGETYPE_OPATTS = 4;
    public final static int CHANGETYPE_PLOTATTS = 5;
    public final static int CHANGETYPE_MOVEOPERATOR = 6;
    public final static int CHANGETYPE_REMOVEOPERATOR = 7;
    public final static int CHANGETYPE_REMOVEALL = 8;
    public final static int CHANGETYPE_REMOVELAST = 9;
    public final static int CHANGETYPE_CACHEINDEX = 10;


    public PlotQueryInfo()
    {
        super(PlotQueryInfo_numAdditionalAtts);

        changeType = CHANGETYPE_NONE;
        oldFrameIndex = 0;
        newFrameIndex = 0;
    }

    public PlotQueryInfo(int nMoreFields)
    {
        super(PlotQueryInfo_numAdditionalAtts + nMoreFields);

        changeType = CHANGETYPE_NONE;
        oldFrameIndex = 0;
        newFrameIndex = 0;
    }

    public PlotQueryInfo(PlotQueryInfo obj)
    {
        super(obj);

        changeType = obj.changeType;
        oldFrameIndex = obj.oldFrameIndex;
        newFrameIndex = obj.newFrameIndex;

        SelectAll();
    }

    public int Offset()
    {
        return super.Offset() + super.GetNumAdditionalAttributes();
    }

    public int GetNumAdditionalAttributes()
    {
        return PlotQueryInfo_numAdditionalAtts;
    }

    public boolean equals(PlotQueryInfo obj)
    {
        // Create the return value
        return ((changeType == obj.changeType) &&
                (oldFrameIndex == obj.oldFrameIndex) &&
                (newFrameIndex == obj.newFrameIndex));
    }

    // Property setting methods
    public void SetChangeType(int changeType_)
    {
        changeType = changeType_;
        Select(0);
    }

    public void SetOldFrameIndex(int oldFrameIndex_)
    {
        oldFrameIndex = oldFrameIndex_;
        Select(1);
    }

    public void SetNewFrameIndex(int newFrameIndex_)
    {
        newFrameIndex = newFrameIndex_;
        Select(2);
    }

    // Property getting methods
    public int GetChangeType() { return changeType; }
    public int GetOldFrameIndex() { return oldFrameIndex; }
    public int GetNewFrameIndex() { return newFrameIndex; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteInt(changeType);
        if(WriteSelect(1, buf))
            buf.WriteInt(oldFrameIndex);
        if(WriteSelect(2, buf))
            buf.WriteInt(newFrameIndex);
    }

    public void ReadAtts(int index, CommunicationBuffer buf)
    {
        switch(index)
        {
        case 0:
            SetChangeType(buf.ReadInt());
            break;
        case 1:
            SetOldFrameIndex(buf.ReadInt());
            break;
        case 2:
            SetNewFrameIndex(buf.ReadInt());
            break;
        }
    }

    public String toString(String indent)
    {
        String str = new String();
        str = str + indent + "changeType = ";
        if(changeType == CHANGETYPE_NONE)
            str = str + "CHANGETYPE_NONE";
        if(changeType == CHANGETYPE_DATABASE)
            str = str + "CHANGETYPE_DATABASE";
        if(changeType == CHANGETYPE_VARNAME)
            str = str + "CHANGETYPE_VARNAME";
        if(changeType == CHANGETYPE_ADDOP)
            str = str + "CHANGETYPE_ADDOP";
        if(changeType == CHANGETYPE_OPATTS)
            str = str + "CHANGETYPE_OPATTS";
        if(changeType == CHANGETYPE_PLOTATTS)
            str = str + "CHANGETYPE_PLOTATTS";
        if(changeType == CHANGETYPE_MOVEOPERATOR)
            str = str + "CHANGETYPE_MOVEOPERATOR";
        if(changeType == CHANGETYPE_REMOVEOPERATOR)
            str = str + "CHANGETYPE_REMOVEOPERATOR";
        if(changeType == CHANGETYPE_REMOVEALL)
            str = str + "CHANGETYPE_REMOVEALL";
        if(changeType == CHANGETYPE_REMOVELAST)
            str = str + "CHANGETYPE_REMOVELAST";
        if(changeType == CHANGETYPE_CACHEINDEX)
            str = str + "CHANGETYPE_CACHEINDEX";
        str = str + "\n";
        str = str + intToString("oldFrameIndex", oldFrameIndex, indent) + "\n";
        str = str + intToString("newFrameIndex", newFrameIndex, indent) + "\n";
        return str;
    }


    // Attributes
    private int changeType;
    private int oldFrameIndex;
    private int newFrameIndex;
}

