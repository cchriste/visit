// ***************************************************************************
//
// Copyright (c) 2000 - 2015, Lawrence Livermore National Security, LLC
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

import java.lang.Integer;
import java.util.Vector;

// ****************************************************************************
// Class: ViewerClientAttributes
//
// Purpose:
//    This class contains attributes used for viewer client
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

public class ViewerClientAttributes extends AttributeSubject
{
    private static int ViewerClientAttributes_numAdditionalAtts = 8;

    // Enum values
    public final static int RENDERTYPE_NONE = 0;
    public final static int RENDERTYPE_IMAGE = 1;
    public final static int RENDERTYPE_DATA = 2;


    public ViewerClientAttributes()
    {
        super(ViewerClientAttributes_numAdditionalAtts);

        renderingType = RENDERTYPE_NONE;
        id = -1;
        title = new String("client");
        windowIds = new Vector();
        imageWidth = 300;
        imageHeight = 300;
        imageResolutionPcnt = 100;
        externalClient = false;
    }

    public ViewerClientAttributes(int nMoreFields)
    {
        super(ViewerClientAttributes_numAdditionalAtts + nMoreFields);

        renderingType = RENDERTYPE_NONE;
        id = -1;
        title = new String("client");
        windowIds = new Vector();
        imageWidth = 300;
        imageHeight = 300;
        imageResolutionPcnt = 100;
        externalClient = false;
    }

    public ViewerClientAttributes(ViewerClientAttributes obj)
    {
        super(ViewerClientAttributes_numAdditionalAtts);

        int i;

        renderingType = obj.renderingType;
        id = obj.id;
        title = new String(obj.title);
        windowIds = new Vector();
        for(i = 0; i < obj.windowIds.size(); ++i)
        {
            Integer iv = (Integer)obj.windowIds.elementAt(i);
            windowIds.addElement(new Integer(iv.intValue()));
        }
        imageWidth = obj.imageWidth;
        imageHeight = obj.imageHeight;
        imageResolutionPcnt = obj.imageResolutionPcnt;
        externalClient = obj.externalClient;

        SelectAll();
    }

    public int Offset()
    {
        return super.Offset() + super.GetNumAdditionalAttributes();
    }

    public int GetNumAdditionalAttributes()
    {
        return ViewerClientAttributes_numAdditionalAtts;
    }

    public boolean equals(ViewerClientAttributes obj)
    {
        int i;

        // Compare the elements in the windowIds vector.
        boolean windowIds_equal = (obj.windowIds.size() == windowIds.size());
        for(i = 0; (i < windowIds.size()) && windowIds_equal; ++i)
        {
            // Make references to Integer from Object.
            Integer windowIds1 = (Integer)windowIds.elementAt(i);
            Integer windowIds2 = (Integer)obj.windowIds.elementAt(i);
            windowIds_equal = windowIds1.equals(windowIds2);
        }
        // Create the return value
        return ((renderingType == obj.renderingType) &&
                (id == obj.id) &&
                (title.equals(obj.title)) &&
                windowIds_equal &&
                (imageWidth == obj.imageWidth) &&
                (imageHeight == obj.imageHeight) &&
                (imageResolutionPcnt == obj.imageResolutionPcnt) &&
                (externalClient == obj.externalClient));
    }

    // Property setting methods
    public void SetRenderingType(int renderingType_)
    {
        renderingType = renderingType_;
        Select(0);
    }

    public void SetId(int id_)
    {
        id = id_;
        Select(1);
    }

    public void SetTitle(String title_)
    {
        title = title_;
        Select(2);
    }

    public void SetWindowIds(Vector windowIds_)
    {
        windowIds = windowIds_;
        Select(3);
    }

    public void SetImageWidth(int imageWidth_)
    {
        imageWidth = imageWidth_;
        Select(4);
    }

    public void SetImageHeight(int imageHeight_)
    {
        imageHeight = imageHeight_;
        Select(5);
    }

    public void SetImageResolutionPcnt(double imageResolutionPcnt_)
    {
        imageResolutionPcnt = imageResolutionPcnt_;
        Select(6);
    }

    public void SetExternalClient(boolean externalClient_)
    {
        externalClient = externalClient_;
        Select(7);
    }

    // Property getting methods
    public int     GetRenderingType() { return renderingType; }
    public int     GetId() { return id; }
    public String  GetTitle() { return title; }
    public Vector  GetWindowIds() { return windowIds; }
    public int     GetImageWidth() { return imageWidth; }
    public int     GetImageHeight() { return imageHeight; }
    public double  GetImageResolutionPcnt() { return imageResolutionPcnt; }
    public boolean GetExternalClient() { return externalClient; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteInt(renderingType);
        if(WriteSelect(1, buf))
            buf.WriteInt(id);
        if(WriteSelect(2, buf))
            buf.WriteString(title);
        if(WriteSelect(3, buf))
            buf.WriteIntVector(windowIds);
        if(WriteSelect(4, buf))
            buf.WriteInt(imageWidth);
        if(WriteSelect(5, buf))
            buf.WriteInt(imageHeight);
        if(WriteSelect(6, buf))
            buf.WriteDouble(imageResolutionPcnt);
        if(WriteSelect(7, buf))
            buf.WriteBool(externalClient);
    }

    public void ReadAtts(int index, CommunicationBuffer buf)
    {
        switch(index)
        {
        case 0:
            SetRenderingType(buf.ReadInt());
            break;
        case 1:
            SetId(buf.ReadInt());
            break;
        case 2:
            SetTitle(buf.ReadString());
            break;
        case 3:
            SetWindowIds(buf.ReadIntVector());
            break;
        case 4:
            SetImageWidth(buf.ReadInt());
            break;
        case 5:
            SetImageHeight(buf.ReadInt());
            break;
        case 6:
            SetImageResolutionPcnt(buf.ReadDouble());
            break;
        case 7:
            SetExternalClient(buf.ReadBool());
            break;
        }
    }

    public String toString(String indent)
    {
        String str = new String();
        str = str + indent + "renderingType = ";
        if(renderingType == RENDERTYPE_NONE)
            str = str + "RENDERTYPE_NONE";
        if(renderingType == RENDERTYPE_IMAGE)
            str = str + "RENDERTYPE_IMAGE";
        if(renderingType == RENDERTYPE_DATA)
            str = str + "RENDERTYPE_DATA";
        str = str + "\n";
        str = str + intToString("id", id, indent) + "\n";
        str = str + stringToString("title", title, indent) + "\n";
        str = str + intVectorToString("windowIds", windowIds, indent) + "\n";
        str = str + intToString("imageWidth", imageWidth, indent) + "\n";
        str = str + intToString("imageHeight", imageHeight, indent) + "\n";
        str = str + doubleToString("imageResolutionPcnt", imageResolutionPcnt, indent) + "\n";
        str = str + boolToString("externalClient", externalClient, indent) + "\n";
        return str;
    }


    // Attributes
    private int     renderingType;
    private int     id;
    private String  title;
    private Vector  windowIds; // vector of Integer objects
    private int     imageWidth;
    private int     imageHeight;
    private double  imageResolutionPcnt;
    private boolean externalClient;
}

