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
// Class: View3DAttributes
//
// Purpose:
//    This class contains the 3d view attributes.
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   Mon Feb 25 15:14:56 PST 2008
//
// Modifications:
//   
// ****************************************************************************

public class View3DAttributes extends AttributeSubject
{
    public View3DAttributes()
    {
        super(13);

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
        parallelScale = 1;
        nearPlane = 0.001;
        farPlane = 100;
        imagePan = new double[2];
        imagePan[0] = 0;
        imagePan[1] = 0;
        imageZoom = 1;
        perspective = true;
        eyeAngle = 2;
        centerOfRotationSet = false;
        centerOfRotation = new double[3];
        centerOfRotation[0] = 0;
        centerOfRotation[1] = 0;
        centerOfRotation[2] = 0;
    }

    public View3DAttributes(View3DAttributes obj)
    {
        super(13);

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
        parallelScale = obj.parallelScale;
        nearPlane = obj.nearPlane;
        farPlane = obj.farPlane;
        imagePan = new double[2];
        imagePan[0] = obj.imagePan[0];
        imagePan[1] = obj.imagePan[1];

        imageZoom = obj.imageZoom;
        perspective = obj.perspective;
        eyeAngle = obj.eyeAngle;
        centerOfRotationSet = obj.centerOfRotationSet;
        centerOfRotation = new double[3];
        centerOfRotation[0] = obj.centerOfRotation[0];
        centerOfRotation[1] = obj.centerOfRotation[1];
        centerOfRotation[2] = obj.centerOfRotation[2];


        SelectAll();
    }

    public boolean equals(View3DAttributes obj)
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

        // Compare the centerOfRotation arrays.
        boolean centerOfRotation_equal = true;
        for(i = 0; i < 3 && centerOfRotation_equal; ++i)
            centerOfRotation_equal = (centerOfRotation[i] == obj.centerOfRotation[i]);

        // Create the return value
        return (viewNormal_equal &&
                focus_equal &&
                viewUp_equal &&
                (viewAngle == obj.viewAngle) &&
                (parallelScale == obj.parallelScale) &&
                (nearPlane == obj.nearPlane) &&
                (farPlane == obj.farPlane) &&
                imagePan_equal &&
                (imageZoom == obj.imageZoom) &&
                (perspective == obj.perspective) &&
                (eyeAngle == obj.eyeAngle) &&
                (centerOfRotationSet == obj.centerOfRotationSet) &&
                centerOfRotation_equal);
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

    public void SetParallelScale(double parallelScale_)
    {
        parallelScale = parallelScale_;
        Select(4);
    }

    public void SetNearPlane(double nearPlane_)
    {
        nearPlane = nearPlane_;
        Select(5);
    }

    public void SetFarPlane(double farPlane_)
    {
        farPlane = farPlane_;
        Select(6);
    }

    public void SetImagePan(double[] imagePan_)
    {
        imagePan[0] = imagePan_[0];
        imagePan[1] = imagePan_[1];
        Select(7);
    }

    public void SetImagePan(double e0, double e1)
    {
        imagePan[0] = e0;
        imagePan[1] = e1;
        Select(7);
    }

    public void SetImageZoom(double imageZoom_)
    {
        imageZoom = imageZoom_;
        Select(8);
    }

    public void SetPerspective(boolean perspective_)
    {
        perspective = perspective_;
        Select(9);
    }

    public void SetEyeAngle(double eyeAngle_)
    {
        eyeAngle = eyeAngle_;
        Select(10);
    }

    public void SetCenterOfRotationSet(boolean centerOfRotationSet_)
    {
        centerOfRotationSet = centerOfRotationSet_;
        Select(11);
    }

    public void SetCenterOfRotation(double[] centerOfRotation_)
    {
        centerOfRotation[0] = centerOfRotation_[0];
        centerOfRotation[1] = centerOfRotation_[1];
        centerOfRotation[2] = centerOfRotation_[2];
        Select(12);
    }

    public void SetCenterOfRotation(double e0, double e1, double e2)
    {
        centerOfRotation[0] = e0;
        centerOfRotation[1] = e1;
        centerOfRotation[2] = e2;
        Select(12);
    }

    // Property getting methods
    public double[] GetViewNormal() { return viewNormal; }
    public double[] GetFocus() { return focus; }
    public double[] GetViewUp() { return viewUp; }
    public double   GetViewAngle() { return viewAngle; }
    public double   GetParallelScale() { return parallelScale; }
    public double   GetNearPlane() { return nearPlane; }
    public double   GetFarPlane() { return farPlane; }
    public double[] GetImagePan() { return imagePan; }
    public double   GetImageZoom() { return imageZoom; }
    public boolean  GetPerspective() { return perspective; }
    public double   GetEyeAngle() { return eyeAngle; }
    public boolean  GetCenterOfRotationSet() { return centerOfRotationSet; }
    public double[] GetCenterOfRotation() { return centerOfRotation; }

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
            buf.WriteDouble(parallelScale);
        if(WriteSelect(5, buf))
            buf.WriteDouble(nearPlane);
        if(WriteSelect(6, buf))
            buf.WriteDouble(farPlane);
        if(WriteSelect(7, buf))
            buf.WriteDoubleArray(imagePan);
        if(WriteSelect(8, buf))
            buf.WriteDouble(imageZoom);
        if(WriteSelect(9, buf))
            buf.WriteBool(perspective);
        if(WriteSelect(10, buf))
            buf.WriteDouble(eyeAngle);
        if(WriteSelect(11, buf))
            buf.WriteBool(centerOfRotationSet);
        if(WriteSelect(12, buf))
            buf.WriteDoubleArray(centerOfRotation);
    }

    public void ReadAtts(int n, CommunicationBuffer buf)
    {
        for(int i = 0; i < n; ++i)
        {
            int index = (int)buf.ReadByte();
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
                SetParallelScale(buf.ReadDouble());
                break;
            case 5:
                SetNearPlane(buf.ReadDouble());
                break;
            case 6:
                SetFarPlane(buf.ReadDouble());
                break;
            case 7:
                SetImagePan(buf.ReadDoubleArray());
                break;
            case 8:
                SetImageZoom(buf.ReadDouble());
                break;
            case 9:
                SetPerspective(buf.ReadBool());
                break;
            case 10:
                SetEyeAngle(buf.ReadDouble());
                break;
            case 11:
                SetCenterOfRotationSet(buf.ReadBool());
                break;
            case 12:
                SetCenterOfRotation(buf.ReadDoubleArray());
                break;
            }
        }
    }

    public String toString(String indent)
    {
        String str = new String();
        str = str + doubleArrayToString("viewNormal", viewNormal, indent) + "\n";
        str = str + doubleArrayToString("focus", focus, indent) + "\n";
        str = str + doubleArrayToString("viewUp", viewUp, indent) + "\n";
        str = str + doubleToString("viewAngle", viewAngle, indent) + "\n";
        str = str + doubleToString("parallelScale", parallelScale, indent) + "\n";
        str = str + doubleToString("nearPlane", nearPlane, indent) + "\n";
        str = str + doubleToString("farPlane", farPlane, indent) + "\n";
        str = str + doubleArrayToString("imagePan", imagePan, indent) + "\n";
        str = str + doubleToString("imageZoom", imageZoom, indent) + "\n";
        str = str + boolToString("perspective", perspective, indent) + "\n";
        str = str + doubleToString("eyeAngle", eyeAngle, indent) + "\n";
        str = str + boolToString("centerOfRotationSet", centerOfRotationSet, indent) + "\n";
        str = str + doubleArrayToString("centerOfRotation", centerOfRotation, indent) + "\n";
        return str;
    }


    // Attributes
    private double[] viewNormal;
    private double[] focus;
    private double[] viewUp;
    private double   viewAngle;
    private double   parallelScale;
    private double   nearPlane;
    private double   farPlane;
    private double[] imagePan;
    private double   imageZoom;
    private boolean  perspective;
    private double   eyeAngle;
    private boolean  centerOfRotationSet;
    private double[] centerOfRotation;
}

