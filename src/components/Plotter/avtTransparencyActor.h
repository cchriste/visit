// ************************************************************************* //
//                          avtTransparencyActor.h                           //
// ************************************************************************* //

#ifndef AVT_TRANSPARENCY_ACTOR_H
#define AVT_TRANSPARENCY_ACTOR_H
#include <plotter_exports.h>
#include <vector>


class     vtkActor;
class     vtkAppendPolyData;
class     vtkAxisDepthSort;
class     vtkCamera;
class     vtkDataSet;
class     vtkDataSetMapper;
class     vtkDepthSortPolyData;
class     vtkMatrix4x4;
class     vtkPolyData;
class     vtkPolyDataMapper;
class     vtkRenderer;


// ****************************************************************************
//  Class: avtTransparencyActor
//
//  Purpose:
//      Rendering transparent polygons correctly requires sorting them in front
//      to back order.  Since the AVT model is to keep groups of polygons
//      separated, this sorting would be impossible.  This class collects the
//      inputs from all of the plots and keeps them together so they can be
//      sorted.
//
//  Programmer: Hank Childs
//  Creation:   July 3, 2002
//
//  Modifications:
//
//    Hank Childs, Thu Jul 11 16:02:45 PDT 2002
//    Added the concept of visibility.  This is "different" from useActor
//    because the VisWindow has both concepts and they don't always agree.  At
//    this point it was easier to adopt its concept than re-write that code.
//
//    Hank Childs, Sun Jul 14 15:49:58 PDT 2002
//    Make use of one VTK module to do all the sorting.
//
//    Jeremy Meredith, Fri Jul 26 14:28:09 PDT 2002
//    Made perfect sorting a permanent mode instead of a one-frame mode.
//
//    Kathleen Bonnell, Wed Jul 16 16:39:02 PDT 2003
//    Added ScaleByVector method. 
//
//    Kathleen Bonnell, Wed Dec  3 16:42:38 PST 2003 
//    Added TransparenciesExist method.
//
//    Chris Wojtan, Fri Jun 25 15:15 PDT 2004
//    Added is2Dimensional bool and functions to get and set its value
//
//    Hank Childs, Wed Sep  8 17:55:34 PDT 2004
//    No longer inline is2Dimensional.
//
// ****************************************************************************

class PLOTTER_API avtTransparencyActor
{
  public:
                                     avtTransparencyActor();
    virtual                         ~avtTransparencyActor();

    int                              AddInput(std::vector<vtkDataSet *> &,
                                              std::vector<vtkDataSetMapper *>&,
                                              std::vector<vtkActor *> &);
    void                             ReplaceInput(int,
                                              std::vector<vtkDataSet *> &,
                                              std::vector<vtkDataSetMapper *>&,
                                              std::vector<vtkActor *> &);

    void                             TurnOffInput(int);
    void                             TurnOnInput(int);
    void                             RemoveInput(int);
    void                             InputWasModified(int);
    void                             SetVisibility(int, bool);
    void                             VisibilityOff(void);
    void                             VisibilityOn(void);

    void                             PrepareForRender(vtkCamera *);
    bool                             UsePerfectSort(bool);

    bool                             TransparenciesExist(void);

    void                             AddToRenderer(vtkRenderer *);
    void                             RemoveFromRenderer(vtkRenderer *);
    void                             ScaleByVector(const float vec[3]);

    bool                             GetIs2Dimensional()
                                                    { return is2Dimensional; };
    void                             SetIs2Dimensional(bool val);

  protected:
    std::vector<std::vector <vtkDataSet *> >         datasets;
    std::vector<std::vector <vtkDataSetMapper *> >   mappers;
    std::vector<std::vector <vtkActor *> >           actors;

    std::vector<std::vector <vtkPolyData *> >        preparedDataset;

    std::vector<bool>                                useActor;
    std::vector<bool>                                visibility;
    std::vector<bool>                                lastExecutionActorList;
    bool                                             inputModified;

    vtkAppendPolyData                               *appender;
    vtkActor                                        *myActor;
    vtkPolyDataMapper                               *myMapper;

    vtkAxisDepthSort                                *axisSort;
    vtkDepthSortPolyData                            *perfectSort;
    bool                                             usePerfectSort;
    vtkMatrix4x4                                    *lastCamera;
    bool                                             is2Dimensional;

    void                                             SetUpActor(void);
    void                                             PrepareDataset(int, int);
};


#endif

    
