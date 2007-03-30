import java.lang.ArrayIndexOutOfBoundsException;
import java.util.Vector;
import llnl.visit.AttributeSubject;
import llnl.visit.SimpleObserver;
import llnl.visit.QueryAttributes;

// ****************************************************************************
// Class: TryQuery
//
// Purpose:
//   This example program does a plot and queries some values in it.
//
// Notes:      
//
// Programmer: Brad Whitlock
// Creation:   Tue Oct 1 12:51:29 PDT 2002
//
// Modifications:
//   Brad Whitlock, Thu Dec 12 10:44:31 PDT 2002
//   Updated because of changse to color table methods.
//
//   Brad Whitlock, Thu Jan 2 16:05:48 PST 2003
//   Changed because of Lineout method interface change.
//
// ****************************************************************************

public class TryQuery extends RunViewer implements SimpleObserver
{
    public TryQuery()
    {
        super();
        doUpdate = true;

        // Make this object observe the light attributes.
        viewer.GetQueryAttributes().Attach(this);
    }

    protected void work(String[] args)
    {
        // Try and open a database
        if(viewer.OpenDatabase("localhost:/usr/gapps/visit/data/curv2d.silo"))
        {
            viewer.AddPlot("Mesh", "curvmesh2d");
            viewer.AddPlot("Pseudocolor", "d");
            viewer.DrawPlots();

            // Set the colortable to one that has white at the bottom values.
            viewer.SetActiveContinuousColorTable("calewhite");

            // Create the variable list.
            Vector vars = new Vector();
            vars.addElement(new String("default"));

            // Do some picks.
            viewer.Pick(300, 300, vars);
            viewer.Pick(450, 350, vars);
            viewer.Pick(600, 400, vars);

            // Do some lineouts.
            viewer.Lineout(-4.01261, 1.91818, 2.52975, 3.78323, vars);
            viewer.SetActiveWindow(1);
            viewer.Lineout(-3.89903, 1.79309, 2.91593, 3.40794, vars);

            // Change the window layout.
            viewer.SetWindowLayout(2);
        }
        else
            System.out.println("Could not open the database!");
    }

    public void Update(AttributeSubject s)
    {
        QueryAttributes q = (QueryAttributes)s;
        printResults(q);
    }

    public void SetUpdate(boolean val) { doUpdate = val; }
    public boolean GetUpdate() { return doUpdate; }

    protected void printResults(QueryAttributes q)
    {
        System.out.println("Query name="+q.GetName()+
            " resultString="+q.GetResultsMessage());
    }

    public static void main(String args[])
    {
        TryQuery r = new TryQuery();
        r.run(args);
    }

    private boolean doUpdate;
}
