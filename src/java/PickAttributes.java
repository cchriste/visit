package llnl.visit;

import java.util.Vector;
import java.lang.Integer;

// ****************************************************************************
// Class: PickAttributes
//
// Purpose:
//    This class contains attributes used for pick.
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   Wed Dec 17 16:17:03 PST 2003
//
// Modifications:
//   
// ****************************************************************************

public class PickAttributes extends AttributeSubject
{
    // Constants
    public final static int PICKTYPE_ZONE = 0;
    public final static int PICKTYPE_NODE = 1;
    public final static int PICKTYPE_CURVEZONE = 2;
    public final static int PICKTYPE_CURVENODE = 3;


    public PickAttributes()
    {
        super(36);

        variables = new Vector();
        variables.addElement(new String("default"));
        displayIncidentElements = true;
        showNodeId = true;
        showNodeDomainLogicalCoords = false;
        showNodeBlockLogicalCoords = false;
        showNodePhysicalCoords = false;
        showZoneId = true;
        showZoneDomainLogicalCoords = false;
        showZoneBlockLogicalCoords = false;
        clearWindow = false;
        pickLetter = new String("(null)");
        fulfilled = false;
        pickType = PICKTYPE_ZONE;
        domain = -1;
        elementNumber = -1;
        incidentElements = new Vector();
        timeStep = -1;
        dimension = -1;
        databaseName = new String("(null)");
        activeVariable = new String("(null)");
        pickPoint = new float[3];
        pickPoint[0] = 0f;
        pickPoint[1] = 0f;
        pickPoint[2] = 0f;
        cellPoint = new float[3];
        cellPoint[0] = 0f;
        cellPoint[1] = 0f;
        cellPoint[2] = 0f;
        nodePoint = new float[3];
        nodePoint[0] = 0f;
        nodePoint[1] = 0f;
        nodePoint[2] = 0f;
        plotBounds = new float[6];
        plotBounds[0] = 0f;
        plotBounds[1] = 0f;
        plotBounds[2] = 0f;
        plotBounds[3] = 0f;
        plotBounds[4] = 0f;
        plotBounds[5] = 0f;
        rayPoint1 = new float[3];
        rayPoint1[0] = 0f;
        rayPoint1[1] = 0f;
        rayPoint1[2] = 0f;
        rayPoint2 = new float[3];
        rayPoint2[0] = 0f;
        rayPoint2[1] = 0f;
        rayPoint2[2] = 0f;
        meshInfo = new String("(null)");
        realElementNumber = -1;
        realIncidentElements = new Vector();
        pnodeCoords = new Vector();
        dnodeCoords = new Vector();
        bnodeCoords = new Vector();
        dzoneCoords = new Vector();
        bzoneCoords = new Vector();
        needTransformMessage = false;
        varInfo = new Vector();
    }

    public PickAttributes(PickAttributes obj)
    {
        super(36);

        int i;

        variables = new Vector(obj.variables.size());
        for(i = 0; i < obj.variables.size(); ++i)
            variables.addElement(new String((String)obj.variables.elementAt(i)));

        displayIncidentElements = obj.displayIncidentElements;
        showNodeId = obj.showNodeId;
        showNodeDomainLogicalCoords = obj.showNodeDomainLogicalCoords;
        showNodeBlockLogicalCoords = obj.showNodeBlockLogicalCoords;
        showNodePhysicalCoords = obj.showNodePhysicalCoords;
        showZoneId = obj.showZoneId;
        showZoneDomainLogicalCoords = obj.showZoneDomainLogicalCoords;
        showZoneBlockLogicalCoords = obj.showZoneBlockLogicalCoords;
        clearWindow = obj.clearWindow;
        pickLetter = new String(obj.pickLetter);
        fulfilled = obj.fulfilled;
        pickType = obj.pickType;
        domain = obj.domain;
        elementNumber = obj.elementNumber;
        incidentElements = new Vector();
        for(i = 0; i < obj.incidentElements.size(); ++i)
        {
            Integer iv = (Integer)obj.incidentElements.elementAt(i);
            incidentElements.addElement(new Integer(iv.intValue()));
        }
        timeStep = obj.timeStep;
        dimension = obj.dimension;
        databaseName = new String(obj.databaseName);
        activeVariable = new String(obj.activeVariable);
        pickPoint = new float[3];
        pickPoint[0] = obj.pickPoint[0];
        pickPoint[1] = obj.pickPoint[1];
        pickPoint[2] = obj.pickPoint[2];

        cellPoint = new float[3];
        cellPoint[0] = obj.cellPoint[0];
        cellPoint[1] = obj.cellPoint[1];
        cellPoint[2] = obj.cellPoint[2];

        nodePoint = new float[3];
        nodePoint[0] = obj.nodePoint[0];
        nodePoint[1] = obj.nodePoint[1];
        nodePoint[2] = obj.nodePoint[2];

        plotBounds = new float[6];
        for(i = 0; i < obj.plotBounds.length; ++i)
            plotBounds[i] = obj.plotBounds[i];

        rayPoint1 = new float[3];
        rayPoint1[0] = obj.rayPoint1[0];
        rayPoint1[1] = obj.rayPoint1[1];
        rayPoint1[2] = obj.rayPoint1[2];

        rayPoint2 = new float[3];
        rayPoint2[0] = obj.rayPoint2[0];
        rayPoint2[1] = obj.rayPoint2[1];
        rayPoint2[2] = obj.rayPoint2[2];

        meshInfo = new String(obj.meshInfo);
        realElementNumber = obj.realElementNumber;
        realIncidentElements = new Vector();
        for(i = 0; i < obj.realIncidentElements.size(); ++i)
        {
            Integer iv = (Integer)obj.realIncidentElements.elementAt(i);
            realIncidentElements.addElement(new Integer(iv.intValue()));
        }
        pnodeCoords = new Vector(obj.pnodeCoords.size());
        for(i = 0; i < obj.pnodeCoords.size(); ++i)
            pnodeCoords.addElement(new String((String)obj.pnodeCoords.elementAt(i)));

        dnodeCoords = new Vector(obj.dnodeCoords.size());
        for(i = 0; i < obj.dnodeCoords.size(); ++i)
            dnodeCoords.addElement(new String((String)obj.dnodeCoords.elementAt(i)));

        bnodeCoords = new Vector(obj.bnodeCoords.size());
        for(i = 0; i < obj.bnodeCoords.size(); ++i)
            bnodeCoords.addElement(new String((String)obj.bnodeCoords.elementAt(i)));

        dzoneCoords = new Vector(obj.dzoneCoords.size());
        for(i = 0; i < obj.dzoneCoords.size(); ++i)
            dzoneCoords.addElement(new String((String)obj.dzoneCoords.elementAt(i)));

        bzoneCoords = new Vector(obj.bzoneCoords.size());
        for(i = 0; i < obj.bzoneCoords.size(); ++i)
            bzoneCoords.addElement(new String((String)obj.bzoneCoords.elementAt(i)));

        needTransformMessage = obj.needTransformMessage;
        // *** Copy the varInfo field ***
        varInfo = new Vector(obj.varInfo.size());
        for(i = 0; i < obj.varInfo.size(); ++i)
        {
            PickVarInfo newObj = (PickVarInfo)varInfo.elementAt(i);
            varInfo.addElement(new PickVarInfo(newObj));
        }


        SelectAll();
    }

    public boolean equals(PickAttributes obj)
    {
        int i;

        // Compare the pickPoint arrays.
        boolean pickPoint_equal = true;
        for(i = 0; i < 3 && pickPoint_equal; ++i)
            pickPoint_equal = (pickPoint[i] == obj.pickPoint[i]);

        // Compare the cellPoint arrays.
        boolean cellPoint_equal = true;
        for(i = 0; i < 3 && cellPoint_equal; ++i)
            cellPoint_equal = (cellPoint[i] == obj.cellPoint[i]);

        // Compare the nodePoint arrays.
        boolean nodePoint_equal = true;
        for(i = 0; i < 3 && nodePoint_equal; ++i)
            nodePoint_equal = (nodePoint[i] == obj.nodePoint[i]);

        // Compare the plotBounds arrays.
        boolean plotBounds_equal = true;
        for(i = 0; i < 6 && plotBounds_equal; ++i)
            plotBounds_equal = (plotBounds[i] == obj.plotBounds[i]);

        // Compare the rayPoint1 arrays.
        boolean rayPoint1_equal = true;
        for(i = 0; i < 3 && rayPoint1_equal; ++i)
            rayPoint1_equal = (rayPoint1[i] == obj.rayPoint1[i]);

        // Compare the rayPoint2 arrays.
        boolean rayPoint2_equal = true;
        for(i = 0; i < 3 && rayPoint2_equal; ++i)
            rayPoint2_equal = (rayPoint2[i] == obj.rayPoint2[i]);

        boolean varInfo_equal = (obj.varInfo.size() == varInfo.size());
        for(i = 0; (i < varInfo.size()) && varInfo_equal; ++i)
        {
            // Make references to PickVarInfo from Object.
            PickVarInfo varInfo1 = (PickVarInfo)varInfo.elementAt(i);
            PickVarInfo varInfo2 = (PickVarInfo)obj.varInfo.elementAt(i);
            varInfo_equal = varInfo1.equals(varInfo2);
        }

        // Create the return value
        return ((variables == obj.variables) &&
                (displayIncidentElements == obj.displayIncidentElements) &&
                (showNodeId == obj.showNodeId) &&
                (showNodeDomainLogicalCoords == obj.showNodeDomainLogicalCoords) &&
                (showNodeBlockLogicalCoords == obj.showNodeBlockLogicalCoords) &&
                (showNodePhysicalCoords == obj.showNodePhysicalCoords) &&
                (showZoneId == obj.showZoneId) &&
                (showZoneDomainLogicalCoords == obj.showZoneDomainLogicalCoords) &&
                (showZoneBlockLogicalCoords == obj.showZoneBlockLogicalCoords) &&
                (clearWindow == obj.clearWindow) &&
                (pickLetter == obj.pickLetter) &&
                (fulfilled == obj.fulfilled) &&
                (pickType == obj.pickType) &&
                (domain == obj.domain) &&
                (elementNumber == obj.elementNumber) &&
                (incidentElements == obj.incidentElements) &&
                (timeStep == obj.timeStep) &&
                (dimension == obj.dimension) &&
                (databaseName == obj.databaseName) &&
                (activeVariable == obj.activeVariable) &&
                pickPoint_equal &&
                cellPoint_equal &&
                nodePoint_equal &&
                plotBounds_equal &&
                rayPoint1_equal &&
                rayPoint2_equal &&
                (meshInfo == obj.meshInfo) &&
                (realElementNumber == obj.realElementNumber) &&
                (realIncidentElements == obj.realIncidentElements) &&
                (pnodeCoords == obj.pnodeCoords) &&
                (dnodeCoords == obj.dnodeCoords) &&
                (bnodeCoords == obj.bnodeCoords) &&
                (dzoneCoords == obj.dzoneCoords) &&
                (bzoneCoords == obj.bzoneCoords) &&
                (needTransformMessage == obj.needTransformMessage) &&
                varInfo_equal);
    }

    // Property setting methods
    public void SetVariables(Vector variables_)
    {
        variables = variables_;
        Select(0);
    }

    public void SetDisplayIncidentElements(boolean displayIncidentElements_)
    {
        displayIncidentElements = displayIncidentElements_;
        Select(1);
    }

    public void SetShowNodeId(boolean showNodeId_)
    {
        showNodeId = showNodeId_;
        Select(2);
    }

    public void SetShowNodeDomainLogicalCoords(boolean showNodeDomainLogicalCoords_)
    {
        showNodeDomainLogicalCoords = showNodeDomainLogicalCoords_;
        Select(3);
    }

    public void SetShowNodeBlockLogicalCoords(boolean showNodeBlockLogicalCoords_)
    {
        showNodeBlockLogicalCoords = showNodeBlockLogicalCoords_;
        Select(4);
    }

    public void SetShowNodePhysicalCoords(boolean showNodePhysicalCoords_)
    {
        showNodePhysicalCoords = showNodePhysicalCoords_;
        Select(5);
    }

    public void SetShowZoneId(boolean showZoneId_)
    {
        showZoneId = showZoneId_;
        Select(6);
    }

    public void SetShowZoneDomainLogicalCoords(boolean showZoneDomainLogicalCoords_)
    {
        showZoneDomainLogicalCoords = showZoneDomainLogicalCoords_;
        Select(7);
    }

    public void SetShowZoneBlockLogicalCoords(boolean showZoneBlockLogicalCoords_)
    {
        showZoneBlockLogicalCoords = showZoneBlockLogicalCoords_;
        Select(8);
    }

    public void SetClearWindow(boolean clearWindow_)
    {
        clearWindow = clearWindow_;
        Select(9);
    }

    public void SetPickLetter(String pickLetter_)
    {
        pickLetter = pickLetter_;
        Select(10);
    }

    public void SetFulfilled(boolean fulfilled_)
    {
        fulfilled = fulfilled_;
        Select(11);
    }

    public void SetPickType(int pickType_)
    {
        pickType = pickType_;
        Select(12);
    }

    public void SetDomain(int domain_)
    {
        domain = domain_;
        Select(13);
    }

    public void SetElementNumber(int elementNumber_)
    {
        elementNumber = elementNumber_;
        Select(14);
    }

    public void SetIncidentElements(Vector incidentElements_)
    {
        incidentElements = incidentElements_;
        Select(15);
    }

    public void SetTimeStep(int timeStep_)
    {
        timeStep = timeStep_;
        Select(16);
    }

    public void SetDimension(int dimension_)
    {
        dimension = dimension_;
        Select(17);
    }

    public void SetDatabaseName(String databaseName_)
    {
        databaseName = databaseName_;
        Select(18);
    }

    public void SetActiveVariable(String activeVariable_)
    {
        activeVariable = activeVariable_;
        Select(19);
    }

    public void SetPickPoint(float[] pickPoint_)
    {
        pickPoint[0] = pickPoint_[0];
        pickPoint[1] = pickPoint_[1];
        pickPoint[2] = pickPoint_[2];
        Select(20);
    }

    public void SetPickPoint(float e0, float e1, float e2)
    {
        pickPoint[0] = e0;
        pickPoint[1] = e1;
        pickPoint[2] = e2;
        Select(20);
    }

    public void SetCellPoint(float[] cellPoint_)
    {
        cellPoint[0] = cellPoint_[0];
        cellPoint[1] = cellPoint_[1];
        cellPoint[2] = cellPoint_[2];
        Select(21);
    }

    public void SetCellPoint(float e0, float e1, float e2)
    {
        cellPoint[0] = e0;
        cellPoint[1] = e1;
        cellPoint[2] = e2;
        Select(21);
    }

    public void SetNodePoint(float[] nodePoint_)
    {
        nodePoint[0] = nodePoint_[0];
        nodePoint[1] = nodePoint_[1];
        nodePoint[2] = nodePoint_[2];
        Select(22);
    }

    public void SetNodePoint(float e0, float e1, float e2)
    {
        nodePoint[0] = e0;
        nodePoint[1] = e1;
        nodePoint[2] = e2;
        Select(22);
    }

    public void SetPlotBounds(float[] plotBounds_)
    {
        for(int i = 0; i < 6; ++i)
             plotBounds[i] = plotBounds_[i];
        Select(23);
    }

    public void SetRayPoint1(float[] rayPoint1_)
    {
        rayPoint1[0] = rayPoint1_[0];
        rayPoint1[1] = rayPoint1_[1];
        rayPoint1[2] = rayPoint1_[2];
        Select(24);
    }

    public void SetRayPoint1(float e0, float e1, float e2)
    {
        rayPoint1[0] = e0;
        rayPoint1[1] = e1;
        rayPoint1[2] = e2;
        Select(24);
    }

    public void SetRayPoint2(float[] rayPoint2_)
    {
        rayPoint2[0] = rayPoint2_[0];
        rayPoint2[1] = rayPoint2_[1];
        rayPoint2[2] = rayPoint2_[2];
        Select(25);
    }

    public void SetRayPoint2(float e0, float e1, float e2)
    {
        rayPoint2[0] = e0;
        rayPoint2[1] = e1;
        rayPoint2[2] = e2;
        Select(25);
    }

    public void SetMeshInfo(String meshInfo_)
    {
        meshInfo = meshInfo_;
        Select(26);
    }

    public void SetRealElementNumber(int realElementNumber_)
    {
        realElementNumber = realElementNumber_;
        Select(27);
    }

    public void SetRealIncidentElements(Vector realIncidentElements_)
    {
        realIncidentElements = realIncidentElements_;
        Select(28);
    }

    public void SetPnodeCoords(Vector pnodeCoords_)
    {
        pnodeCoords = pnodeCoords_;
        Select(29);
    }

    public void SetDnodeCoords(Vector dnodeCoords_)
    {
        dnodeCoords = dnodeCoords_;
        Select(30);
    }

    public void SetBnodeCoords(Vector bnodeCoords_)
    {
        bnodeCoords = bnodeCoords_;
        Select(31);
    }

    public void SetDzoneCoords(Vector dzoneCoords_)
    {
        dzoneCoords = dzoneCoords_;
        Select(32);
    }

    public void SetBzoneCoords(Vector bzoneCoords_)
    {
        bzoneCoords = bzoneCoords_;
        Select(33);
    }

    public void SetNeedTransformMessage(boolean needTransformMessage_)
    {
        needTransformMessage = needTransformMessage_;
        Select(34);
    }

    // Property getting methods
    public Vector  GetVariables() { return variables; }
    public boolean GetDisplayIncidentElements() { return displayIncidentElements; }
    public boolean GetShowNodeId() { return showNodeId; }
    public boolean GetShowNodeDomainLogicalCoords() { return showNodeDomainLogicalCoords; }
    public boolean GetShowNodeBlockLogicalCoords() { return showNodeBlockLogicalCoords; }
    public boolean GetShowNodePhysicalCoords() { return showNodePhysicalCoords; }
    public boolean GetShowZoneId() { return showZoneId; }
    public boolean GetShowZoneDomainLogicalCoords() { return showZoneDomainLogicalCoords; }
    public boolean GetShowZoneBlockLogicalCoords() { return showZoneBlockLogicalCoords; }
    public boolean GetClearWindow() { return clearWindow; }
    public String  GetPickLetter() { return pickLetter; }
    public boolean GetFulfilled() { return fulfilled; }
    public int     GetPickType() { return pickType; }
    public int     GetDomain() { return domain; }
    public int     GetElementNumber() { return elementNumber; }
    public Vector  GetIncidentElements() { return incidentElements; }
    public int     GetTimeStep() { return timeStep; }
    public int     GetDimension() { return dimension; }
    public String  GetDatabaseName() { return databaseName; }
    public String  GetActiveVariable() { return activeVariable; }
    public float[] GetPickPoint() { return pickPoint; }
    public float[] GetCellPoint() { return cellPoint; }
    public float[] GetNodePoint() { return nodePoint; }
    public float[] GetPlotBounds() { return plotBounds; }
    public float[] GetRayPoint1() { return rayPoint1; }
    public float[] GetRayPoint2() { return rayPoint2; }
    public String  GetMeshInfo() { return meshInfo; }
    public int     GetRealElementNumber() { return realElementNumber; }
    public Vector  GetRealIncidentElements() { return realIncidentElements; }
    public Vector  GetPnodeCoords() { return pnodeCoords; }
    public Vector  GetDnodeCoords() { return dnodeCoords; }
    public Vector  GetBnodeCoords() { return bnodeCoords; }
    public Vector  GetDzoneCoords() { return dzoneCoords; }
    public Vector  GetBzoneCoords() { return bzoneCoords; }
    public boolean GetNeedTransformMessage() { return needTransformMessage; }
    public Vector  GetVarInfo() { return varInfo; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteStringVector(variables);
        if(WriteSelect(1, buf))
            buf.WriteBool(displayIncidentElements);
        if(WriteSelect(2, buf))
            buf.WriteBool(showNodeId);
        if(WriteSelect(3, buf))
            buf.WriteBool(showNodeDomainLogicalCoords);
        if(WriteSelect(4, buf))
            buf.WriteBool(showNodeBlockLogicalCoords);
        if(WriteSelect(5, buf))
            buf.WriteBool(showNodePhysicalCoords);
        if(WriteSelect(6, buf))
            buf.WriteBool(showZoneId);
        if(WriteSelect(7, buf))
            buf.WriteBool(showZoneDomainLogicalCoords);
        if(WriteSelect(8, buf))
            buf.WriteBool(showZoneBlockLogicalCoords);
        if(WriteSelect(9, buf))
            buf.WriteBool(clearWindow);
        if(WriteSelect(10, buf))
            buf.WriteString(pickLetter);
        if(WriteSelect(11, buf))
            buf.WriteBool(fulfilled);
        if(WriteSelect(12, buf))
            buf.WriteInt(pickType);
        if(WriteSelect(13, buf))
            buf.WriteInt(domain);
        if(WriteSelect(14, buf))
            buf.WriteInt(elementNumber);
        if(WriteSelect(15, buf))
            buf.WriteIntVector(incidentElements);
        if(WriteSelect(16, buf))
            buf.WriteInt(timeStep);
        if(WriteSelect(17, buf))
            buf.WriteInt(dimension);
        if(WriteSelect(18, buf))
            buf.WriteString(databaseName);
        if(WriteSelect(19, buf))
            buf.WriteString(activeVariable);
        if(WriteSelect(20, buf))
            buf.WriteFloatArray(pickPoint);
        if(WriteSelect(21, buf))
            buf.WriteFloatArray(cellPoint);
        if(WriteSelect(22, buf))
            buf.WriteFloatArray(nodePoint);
        if(WriteSelect(23, buf))
            buf.WriteFloatArray(plotBounds);
        if(WriteSelect(24, buf))
            buf.WriteFloatArray(rayPoint1);
        if(WriteSelect(25, buf))
            buf.WriteFloatArray(rayPoint2);
        if(WriteSelect(26, buf))
            buf.WriteString(meshInfo);
        if(WriteSelect(27, buf))
            buf.WriteInt(realElementNumber);
        if(WriteSelect(28, buf))
            buf.WriteIntVector(realIncidentElements);
        if(WriteSelect(29, buf))
            buf.WriteStringVector(pnodeCoords);
        if(WriteSelect(30, buf))
            buf.WriteStringVector(dnodeCoords);
        if(WriteSelect(31, buf))
            buf.WriteStringVector(bnodeCoords);
        if(WriteSelect(32, buf))
            buf.WriteStringVector(dzoneCoords);
        if(WriteSelect(33, buf))
            buf.WriteStringVector(bzoneCoords);
        if(WriteSelect(34, buf))
            buf.WriteBool(needTransformMessage);
        if(WriteSelect(35, buf))
        {
            buf.WriteInt(varInfo.size());
            for(int i = 0; i < varInfo.size(); ++i)
            {
                PickVarInfo tmp = (PickVarInfo)varInfo.elementAt(i);
                tmp.Write(buf);
            }
        }
    }

    public void ReadAtts(int n, CommunicationBuffer buf)
    {
        for(int i = 0; i < n; ++i)
        {
            int index = (int)buf.ReadByte();
            switch(index)
            {
            case 0:
                SetVariables(buf.ReadStringVector());
                break;
            case 1:
                SetDisplayIncidentElements(buf.ReadBool());
                break;
            case 2:
                SetShowNodeId(buf.ReadBool());
                break;
            case 3:
                SetShowNodeDomainLogicalCoords(buf.ReadBool());
                break;
            case 4:
                SetShowNodeBlockLogicalCoords(buf.ReadBool());
                break;
            case 5:
                SetShowNodePhysicalCoords(buf.ReadBool());
                break;
            case 6:
                SetShowZoneId(buf.ReadBool());
                break;
            case 7:
                SetShowZoneDomainLogicalCoords(buf.ReadBool());
                break;
            case 8:
                SetShowZoneBlockLogicalCoords(buf.ReadBool());
                break;
            case 9:
                SetClearWindow(buf.ReadBool());
                break;
            case 10:
                SetPickLetter(buf.ReadString());
                break;
            case 11:
                SetFulfilled(buf.ReadBool());
                break;
            case 12:
                SetPickType(buf.ReadInt());
                break;
            case 13:
                SetDomain(buf.ReadInt());
                break;
            case 14:
                SetElementNumber(buf.ReadInt());
                break;
            case 15:
                SetIncidentElements(buf.ReadIntVector());
                break;
            case 16:
                SetTimeStep(buf.ReadInt());
                break;
            case 17:
                SetDimension(buf.ReadInt());
                break;
            case 18:
                SetDatabaseName(buf.ReadString());
                break;
            case 19:
                SetActiveVariable(buf.ReadString());
                break;
            case 20:
                SetPickPoint(buf.ReadFloatArray());
                break;
            case 21:
                SetCellPoint(buf.ReadFloatArray());
                break;
            case 22:
                SetNodePoint(buf.ReadFloatArray());
                break;
            case 23:
                SetPlotBounds(buf.ReadFloatArray());
                break;
            case 24:
                SetRayPoint1(buf.ReadFloatArray());
                break;
            case 25:
                SetRayPoint2(buf.ReadFloatArray());
                break;
            case 26:
                SetMeshInfo(buf.ReadString());
                break;
            case 27:
                SetRealElementNumber(buf.ReadInt());
                break;
            case 28:
                SetRealIncidentElements(buf.ReadIntVector());
                break;
            case 29:
                SetPnodeCoords(buf.ReadStringVector());
                break;
            case 30:
                SetDnodeCoords(buf.ReadStringVector());
                break;
            case 31:
                SetBnodeCoords(buf.ReadStringVector());
                break;
            case 32:
                SetDzoneCoords(buf.ReadStringVector());
                break;
            case 33:
                SetBzoneCoords(buf.ReadStringVector());
                break;
            case 34:
                SetNeedTransformMessage(buf.ReadBool());
                break;
            case 35:
                {
                    int len = buf.ReadInt();
                    varInfo.clear();
                    for(int j = 0; j < len; ++j)
                    {
                        PickVarInfo tmp = new PickVarInfo();
                        tmp.Read(buf);
                        varInfo.addElement(tmp);
                    }
                }
                Select(35);
                break;
            }
        }
    }

    // Attributegroup convenience methods
    public void AddPickVarInfo(PickVarInfo obj)
    {
        varInfo.addElement(new PickVarInfo(obj));
        Select(35);
    }

    public void ClearPickVarInfos()
    {
        varInfo.clear();
        Select(35);
    }

    public void RemovePickVarInfo(int index)
    {
        if(index >= 0 && index < varInfo.size())
        {
            varInfo.remove(index);
            Select(35);
        }
    }

    public int GetNumPickVarInfos()
    {
        return varInfo.size();
    }

    public PickVarInfo GetPickVarInfo(int i)
    {
        PickVarInfo tmp = (PickVarInfo)varInfo.elementAt(i);
        return tmp;
    }


    // Attributes
    private Vector  variables; // vector of String objects
    private boolean displayIncidentElements;
    private boolean showNodeId;
    private boolean showNodeDomainLogicalCoords;
    private boolean showNodeBlockLogicalCoords;
    private boolean showNodePhysicalCoords;
    private boolean showZoneId;
    private boolean showZoneDomainLogicalCoords;
    private boolean showZoneBlockLogicalCoords;
    private boolean clearWindow;
    private String  pickLetter;
    private boolean fulfilled;
    private int     pickType;
    private int     domain;
    private int     elementNumber;
    private Vector  incidentElements; // vector of Integer objects
    private int     timeStep;
    private int     dimension;
    private String  databaseName;
    private String  activeVariable;
    private float[] pickPoint;
    private float[] cellPoint;
    private float[] nodePoint;
    private float[] plotBounds;
    private float[] rayPoint1;
    private float[] rayPoint2;
    private String  meshInfo;
    private int     realElementNumber;
    private Vector  realIncidentElements; // vector of Integer objects
    private Vector  pnodeCoords; // vector of String objects
    private Vector  dnodeCoords; // vector of String objects
    private Vector  bnodeCoords; // vector of String objects
    private Vector  dzoneCoords; // vector of String objects
    private Vector  bzoneCoords; // vector of String objects
    private boolean needTransformMessage;
    private Vector  varInfo; // vector of PickVarInfo objects
}

