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
// Class: ViewAttributes
//
// Purpose:
//    This class contains the view attributes.
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

public class ViewAttributes extends AttributeSubject
{
    private static int ViewAttributes_numAdditionalAtts = 14;

    public ViewAttributes()
    {
        super(ViewAttributes_numAdditionalAtts);

        viewNormal = new double[3];
        viewNormal[0] = 0;
        viewNormal[1] = 0;
        viewNormal[2] = 1;
        focus = new double[3];
        focus[0] = 0;
        focus[1] = 0;
        focus[2] = 0;
        viewUp = new double[3];
        viewUp[0] = 0;
        viewUp[1] = 1;
        viewUp[2] = 0;
        viewAngle = 30;
        setScale = false;
        parallelScale = 1;
        nearPlane = 0.001;
        farPlane = 100;
        imagePan = new double[2];
        imagePan[0] = 0;
        imagePan[1] = 0;
        imageZoom = 1;
        perspective = true;
        windowCoords = new double[4];
        windowCoords[0] = 0;
        windowCoords[1] = 0;
        windowCoords[2] = 1;
        windowCoords[3] = 1;
        viewportCoords = new double[4];
        viewportCoords[0] = 0.1;
        viewportCoords[1] = 0.1;
        viewportCoords[2] = 0.9;
        viewportCoords[3] = 0.9;
        eyeAngle = 2;
    }

    public ViewAttributes(int nMoreFields)
    {
        super(ViewAttributes_numAdditionalAtts + nMoreFields);

        viewNormal = new double[3];
        viewNormal[0] = 0;
        viewNormal[1] = 0;
        viewNormal[2] = 1;
        focus = new double[3];
        focus[0] = 0;
        focus[1] = 0;
        focus[2] = 0;
        viewUp = new double[3];
        viewUp[0] = 0;
        viewUp[1] = 1;
        viewUp[2] = 0;
        viewAngle = 30;
        setScale = false;
        parallelScale = 1;
        nearPlane = 0.001;
        farPlane = 100;
        imagePan = new double[2];
        imagePan[0] = 0;
        imagePan[1] = 0;
        imageZoom = 1;
        perspective = true;
        windowCoords = new double[4];
        windowCoords[0] = 0;
        windowCoords[1] = 0;
        windowCoords[2] = 1;
        windowCoords[3] = 1;
        viewportCoords = new double[4];
        viewportCoords[0] = 0.1;
        viewportCoords[1] = 0.1;
        viewportCoords[2] = 0.9;
        viewportCoords[3] = 0.9;
        eyeAngle = 2;
    }

    public ViewAttributes(ViewAttributes obj)
    {
        super(obj);

        int i;

        viewNormal = new double[3];
        viewNormal[0] = obj.viewNormal[0];
        viewNormal[1] = obj.viewNormal[1];
        viewNormal[2] = obj.viewNormal[2];

        focus = new double[3];
        focus[0] = obj.focus[0];
        focus[1] = obj.focus[1];
        focus[2] = obj.focus[2];

        viewUp = new double[3];
        viewUp[0] = obj.viewUp[0];
        viewUp[1] = obj.viewUp[1];
        viewUp[2] = obj.viewUp[2];

        viewAngle = obj.viewAngle;
        setScale = obj.setScale;
        parallelScale = obj.parallelScale;
        nearPlane = obj.nearPlane;
        farPlane = obj.farPlane;
        imagePan = new double[2];
        imagePan[0] = obj.imagePan[0];
        imagePan[1] = obj.imagePan[1];

        imageZoom = obj.imageZoom;
        perspective = obj.perspective;
        windowCoords = new double[4];
        for(i = 0; i < obj.windowCoords.length; ++i)
            windowCoords[i] = obj.windowCoords[i];

        viewportCoords = new double[4];
        for(i = 0; i < obj.viewportCoords.length; ++i)
            viewportCoords[i] = obj.viewportCoords[i];

        eyeAngle = obj.eyeAngle;

        SelectAll();
    }

    public int Offset()
    {
        return super.Offset() + super.GetNumAdditionalAttributes();
    }

    public int GetNumAdditionalAttributes()
    {
        return ViewAttributes_numAdditionalAtts;
    }

    public boolean equals(ViewAttributes obj)
    {
        int i;

        // Compare the viewNormal arrays.
        boolean viewNormal_equal = true;
        for(i = 0; i < 3 && viewNormal_equal; ++i)
            viewNormal_equal = (viewNormal[i] == obj.viewNormal[i]);

        // Compare the focus arrays.
        boolean focus_equal = true;
        for(i = 0; i < 3 && focus_equal; ++i)
            focus_equal = (focus[i] == obj.focus[i]);

        // Compare the viewUp arrays.
        boolean viewUp_equal = true;
        for(i = 0; i < 3 && viewUp_equal; ++i)
            viewUp_equal = (viewUp[i] == obj.viewUp[i]);

        // Compare the imagePan arrays.
        boolean imagePan_equal = true;
        for(i = 0; i < 2 && imagePan_equal; ++i)
            imagePan_equal = (imagePan[i] == obj.imagePan[i]);

        // Compare the windowCoords arrays.
        boolean windowCoords_equal = true;
        for(i = 0; i < 4 && windowCoords_equal; ++i)
            windowCoords_equal = (windowCoords[i] == obj.windowCoords[i]);

        // Compare the viewportCoords arrays.
        boolean viewportCoords_equal = true;
        for(i = 0; i < 4 && viewportCoords_equal; ++i)
            viewportCoords_equal = (viewportCoords[i] == obj.viewportCoords[i]);

        // Create the return value
        return (viewNormal_equal &&
                focus_equal &&
                viewUp_equal &&
                (viewAngle == obj.viewAngle) &&
                (setScale == obj.setScale) &&
                (parallelScale == obj.parallelScale) &&
                (nearPlane == obj.nearPlane) &&
                (farPlane == obj.farPlane) &&
                imagePan_equal &&
                (imageZoom == obj.imageZoom) &&
                (perspective == obj.perspective) &&
                windowCoords_equal &&
                viewportCoords_equal &&
                (eyeAngle == obj.eyeAngle));
    }

    // Property setting methods
    public void SetViewNormal(double[] viewNormal_)
    {
        viewNormal[0] = viewNormal_[0];
        viewNormal[1] = viewNormal_[1];
        viewNormal[2] = viewNormal_[2];
        Select(0);
    }

    public void SetViewNormal(double e0, double e1, double e2)
    {
        viewNormal[0] = e0;
        viewNormal[1] = e1;
        viewNormal[2] = e2;
        Select(0);
    }

    public void SetFocus(double[] focus_)
    {
        focus[0] = focus_[0];
        focus[1] = focus_[1];
        focus[2] = focus_[2];
        Select(1);
    }

    public void SetFocus(double e0, double e1, double e2)
    {
        focus[0] = e0;
        focus[1] = e1;
        focus[2] = e2;
        Select(1);
    }

    public void SetViewUp(double[] viewUp_)
    {
        viewUp[0] = viewUp_[0];
        viewUp[1] = viewUp_[1];
        viewUp[2] = viewUp_[2];
        Select(2);
    }

    public void SetViewUp(double e0, double e1, double e2)
    {
        viewUp[0] = e0;
        viewUp[1] = e1;
        viewUp[2] = e2;
        Select(2);
    }

    public void SetViewAngle(double viewAngle_)
    {
        viewAngle = viewAngle_;
        Select(3);
    }

    public void SetSetScale(boolean setScale_)
    {
        setScale = setScale_;
        Select(4);
    }

    public void SetParallelScale(double parallelScale_)
    {
        parallelScale = parallelScale_;
        Select(5);
    }

    public void SetNearPlane(double nearPlane_)
    {
        nearPlane = nearPlane_;
        Select(6);
    }

    public void SetFarPlane(double farPlane_)
    {
        farPlane = farPlane_;
        Select(7);
    }

    public void SetImagePan(double[] imagePan_)
    {
        imagePan[0] = imagePan_[0];
        imagePan[1] = imagePan_[1];
        Select(8);
    }

    public void SetImagePan(double e0, double e1)
    {
        imagePan[0] = e0;
        imagePan[1] = e1;
        Select(8);
    }

    public void SetImageZoom(double imageZoom_)
    {
        imageZoom = imageZoom_;
        Select(9);
    }

    public void SetPerspective(boolean perspective_)
    {
        perspective = perspective_;
        Select(10);
    }

    public void SetWindowCoords(double[] windowCoords_)
    {
        windowCoords[0] = windowCoords_[0];
        windowCoords[1] = windowCoords_[1];
        windowCoords[2] = windowCoords_[2];
        windowCoords[3] = windowCoords_[3];
        Select(11);
    }

    public void SetWindowCoords(double e0, double e1, double e2, double e3)
    {
        windowCoords[0] = e0;
        windowCoords[1] = e1;
        windowCoords[2] = e2;
        windowCoords[3] = e3;
        Select(11);
    }

    public void SetViewportCoords(double[] viewportCoords_)
    {
        viewportCoords[0] = viewportCoords_[0];
        viewportCoords[1] = viewportCoords_[1];
        viewportCoords[2] = viewportCoords_[2];
        viewportCoords[3] = viewportCoords_[3];
        Select(12);
    }

    public void SetViewportCoords(double e0, double e1, double e2, double e3)
    {
        viewportCoords[0] = e0;
        viewportCoords[1] = e1;
        viewportCoords[2] = e2;
        viewportCoords[3] = e3;
        Select(12);
    }

    public void SetEyeAngle(double eyeAngle_)
    {
        eyeAngle = eyeAngle_;
        Select(13);
    }

    // Property getting methods
    public double[] GetViewNormal() { return viewNormal; }
    public double[] GetFocus() { return focus; }
    public double[] GetViewUp() { return viewUp; }
    public double   GetViewAngle() { return viewAngle; }
    public boolean  GetSetScale() { return setScale; }
    public double   GetParallelScale() { return parallelScale; }
    public double   GetNearPlane() { return nearPlane; }
    public double   GetFarPlane() { return farPlane; }
    public double[] GetImagePan() { return imagePan; }
    public double   GetImageZoom() { return imageZoom; }
    public boolean  GetPerspective() { return perspective; }
    public double[] GetWindowCoords() { return windowCoords; }
    public double[] GetViewportCoords() { return viewportCoords; }
    public double   GetEyeAngle() { return eyeAngle; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteDoubleArray(viewNormal);
        if(WriteSelect(1, buf))
            buf.WriteDoubleArray(focus);
        if(WriteSelect(2, buf))
            buf.WriteDoubleArray(viewUp);
        if(WriteSelect(3, buf))
            buf.WriteDouble(viewAngle);
        if(WriteSelect(4, buf))
            buf.WriteBool(setScale);
        if(WriteSelect(5, buf))
            buf.WriteDouble(parallelScale);
        if(WriteSelect(6, buf))
            buf.WriteDouble(nearPlane);
        if(WriteSelect(7, buf))
            buf.WriteDouble(farPlane);
        if(WriteSelect(8, buf))
            buf.WriteDoubleArray(imagePan);
        if(WriteSelect(9, buf))
            buf.WriteDouble(imageZoom);
        if(WriteSelect(10, buf))
            buf.WriteBool(perspective);
        if(WriteSelect(11, buf))
            buf.WriteDoubleArray(windowCoords);
        if(WriteSelect(12, buf))
            buf.WriteDoubleArray(viewportCoords);
        if(WriteSelect(13, buf))
            buf.WriteDouble(eyeAngle);
    }

    public void ReadAtts(int index, CommunicationBuffer buf)
    {
        switch(index)
        {
        case 0:
            SetViewNormal(buf.ReadDoubleArray());
            break;
        case 1:
            SetFocus(buf.ReadDoubleArray());
            break;
        case 2:
            SetViewUp(buf.ReadDoubleArray());
            break;
        case 3:
            SetViewAngle(buf.ReadDouble());
            break;
        case 4:
            SetSetScale(buf.ReadBool());
            break;
        case 5:
            SetParallelScale(buf.ReadDouble());
            break;
        case 6:
            SetNearPlane(buf.ReadDouble());
            break;
        case 7:
            SetFarPlane(buf.ReadDouble());
            break;
        case 8:
            SetImagePan(buf.ReadDoubleArray());
            break;
        case 9:
            SetImageZoom(buf.ReadDouble());
            break;
        case 10:
            SetPerspective(buf.ReadBool());
            break;
        case 11:
            SetWindowCoords(buf.ReadDoubleArray());
            break;
        case 12:
            SetViewportCoords(buf.ReadDoubleArray());
            break;
        case 13:
            SetEyeAngle(buf.ReadDouble());
            break;
        }
    }

    public String toString(String indent)
    {
        String str = new String();
        str = str + doubleArrayToString("viewNormal", viewNormal, indent) + "\n";
        str = str + doubleArrayToString("focus", focus, indent) + "\n";
        str = str + doubleArrayToString("viewUp", viewUp, indent) + "\n";
        str = str + doubleToString("viewAngle", viewAngle, indent) + "\n";
        str = str + boolToString("setScale", setScale, indent) + "\n";
        str = str + doubleToString("parallelScale", parallelScale, indent) + "\n";
        str = str + doubleToString("nearPlane", nearPlane, indent) + "\n";
        str = str + doubleToString("farPlane", farPlane, indent) + "\n";
        str = str + doubleArrayToString("imagePan", imagePan, indent) + "\n";
        str = str + doubleToString("imageZoom", imageZoom, indent) + "\n";
        str = str + boolToString("perspective", perspective, indent) + "\n";
        str = str + doubleArrayToString("windowCoords", windowCoords, indent) + "\n";
        str = str + doubleArrayToString("viewportCoords", viewportCoords, indent) + "\n";
        str = str + doubleToString("eyeAngle", eyeAngle, indent) + "\n";
        return str;
    }


    // Attributes
    private double[] viewNormal;
    private double[] focus;
    private double[] viewUp;
    private double   viewAngle;
    private boolean  setScale;
    private double   parallelScale;
    private double   nearPlane;
    private double   farPlane;
    private double[] imagePan;
    private double   imageZoom;
    private boolean  perspective;
    private double[] windowCoords;
    private double[] viewportCoords;
    private double   eyeAngle;
}

