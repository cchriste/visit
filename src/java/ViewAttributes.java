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
// Creation:   Thu Aug 15 12:24:43 PDT 2002
//
// Modifications:
//   
// ****************************************************************************

public class ViewAttributes extends AttributeSubject
{
    public ViewAttributes()
    {
        super(11);

        camera = new double[3];
        camera[0] = 0;
        camera[1] = 0;
        camera[2] = -1;
        focus = new double[3];
        focus[0] = 0;
        focus[1] = 0;
        focus[2] = 0;
        viewUp = new double[3];
        viewUp[0] = 0;
        viewUp[1] = 0;
        viewUp[2] = 0;
        viewAngle = 30;
        setScale = false;
        parallelScale = 1;
        nearPlane = 0.001;
        farPlane = 100;
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
    }

    public ViewAttributes(ViewAttributes obj)
    {
        super(11);

        int i;

        camera = new double[3];
        camera[0] = obj.camera[0];
        camera[1] = obj.camera[1];
        camera[2] = obj.camera[2];

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
        perspective = obj.perspective;
        windowCoords = new double[4];
        for(i = 0; i < obj.windowCoords.length; ++i)
            windowCoords[i] = obj.windowCoords[i];

        viewportCoords = new double[4];
        for(i = 0; i < obj.viewportCoords.length; ++i)
            viewportCoords[i] = obj.viewportCoords[i];


        SelectAll();
    }

    public boolean equals(ViewAttributes obj)
    {
        int i;

        // Compare the camera arrays.
        boolean camera_equal = true;
        for(i = 0; i < 3 && camera_equal; ++i)
            camera_equal = (camera[i] == obj.camera[i]);

        // Compare the focus arrays.
        boolean focus_equal = true;
        for(i = 0; i < 3 && focus_equal; ++i)
            focus_equal = (focus[i] == obj.focus[i]);

        // Compare the viewUp arrays.
        boolean viewUp_equal = true;
        for(i = 0; i < 3 && viewUp_equal; ++i)
            viewUp_equal = (viewUp[i] == obj.viewUp[i]);

        // Compare the windowCoords arrays.
        boolean windowCoords_equal = true;
        for(i = 0; i < 4 && windowCoords_equal; ++i)
            windowCoords_equal = (windowCoords[i] == obj.windowCoords[i]);

        // Compare the viewportCoords arrays.
        boolean viewportCoords_equal = true;
        for(i = 0; i < 4 && viewportCoords_equal; ++i)
            viewportCoords_equal = (viewportCoords[i] == obj.viewportCoords[i]);

        // Create the return value
        return (camera_equal &&
                focus_equal &&
                viewUp_equal &&
                (viewAngle == obj.viewAngle) &&
                (setScale == obj.setScale) &&
                (parallelScale == obj.parallelScale) &&
                (nearPlane == obj.nearPlane) &&
                (farPlane == obj.farPlane) &&
                (perspective == obj.perspective) &&
                windowCoords_equal &&
                viewportCoords_equal);
    }

    // Property setting methods
    public void SetCamera(double[] camera_)
    {
        camera[0] = camera_[0];
        camera[1] = camera_[1];
        camera[2] = camera_[2];
        Select(0);
    }

    public void SetCamera(double e0, double e1, double e2)
    {
        camera[0] = e0;
        camera[1] = e1;
        camera[2] = e2;
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

    public void SetPerspective(boolean perspective_)
    {
        perspective = perspective_;
        Select(8);
    }

    public void SetWindowCoords(double[] windowCoords_)
    {
        windowCoords[0] = windowCoords_[0];
        windowCoords[1] = windowCoords_[1];
        windowCoords[2] = windowCoords_[2];
        windowCoords[3] = windowCoords_[3];
        Select(9);
    }

    public void SetWindowCoords(double e0, double e1, double e2, double e3)
    {
        windowCoords[0] = e0;
        windowCoords[1] = e1;
        windowCoords[2] = e2;
        windowCoords[3] = e3;
        Select(9);
    }

    public void SetViewportCoords(double[] viewportCoords_)
    {
        viewportCoords[0] = viewportCoords_[0];
        viewportCoords[1] = viewportCoords_[1];
        viewportCoords[2] = viewportCoords_[2];
        viewportCoords[3] = viewportCoords_[3];
        Select(10);
    }

    public void SetViewportCoords(double e0, double e1, double e2, double e3)
    {
        viewportCoords[0] = e0;
        viewportCoords[1] = e1;
        viewportCoords[2] = e2;
        viewportCoords[3] = e3;
        Select(10);
    }

    // Property getting methods
    public double[] GetCamera() { return camera; }
    public double[] GetFocus() { return focus; }
    public double[] GetViewUp() { return viewUp; }
    public double   GetViewAngle() { return viewAngle; }
    public boolean  GetSetScale() { return setScale; }
    public double   GetParallelScale() { return parallelScale; }
    public double   GetNearPlane() { return nearPlane; }
    public double   GetFarPlane() { return farPlane; }
    public boolean  GetPerspective() { return perspective; }
    public double[] GetWindowCoords() { return windowCoords; }
    public double[] GetViewportCoords() { return viewportCoords; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteDoubleArray(camera);
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
            buf.WriteBool(perspective);
        if(WriteSelect(9, buf))
            buf.WriteDoubleArray(windowCoords);
        if(WriteSelect(10, buf))
            buf.WriteDoubleArray(viewportCoords);
    }

    public void ReadAtts(int n, CommunicationBuffer buf)
    {
        for(int i = 0; i < n; ++i)
        {
            int index = (int)buf.ReadByte();
            switch(index)
            {
            case 0:
                SetCamera(buf.ReadDoubleArray());
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
                SetPerspective(buf.ReadBool());
                break;
            case 9:
                SetWindowCoords(buf.ReadDoubleArray());
                break;
            case 10:
                SetViewportCoords(buf.ReadDoubleArray());
                break;
            }
        }
    }


    // Attributes
    private double[] camera;
    private double[] focus;
    private double[] viewUp;
    private double   viewAngle;
    private boolean  setScale;
    private double   parallelScale;
    private double   nearPlane;
    private double   farPlane;
    private boolean  perspective;
    private double[] windowCoords;
    private double[] viewportCoords;
}

