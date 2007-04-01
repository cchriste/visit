// ************************************************************************* //
//                            avtPickQuery.h                                 //
// ************************************************************************* //

#ifndef AVT_PICK_QUERY_H
#define AVT_PICK_QUERY_H
#include <query_exports.h>

#include <avtDatasetQuery.h>

#include <string>
#include <PickAttributes.h>
#include <avtTypes.h>

class avtExpressionEvaluatorFilter;
class avtMatrix;
class vtkDataSet;


// ****************************************************************************
//  Class: avtPickQuery
//
//  Purpose:
//    This query gathers information about a cell picked by the user. 
//
//  Programmer: Kathleen Bonnell 
//  Creation:   November 15, 2002 
//
//  Modifications:
//    Kathleen Bonnell, Thu Apr 17 09:39:19 PDT 2003  
//    Added member invTransform.
//
//    Kathleen Bonnell, Wed May  7 13:24:37 PDT 2003 
//    Renamed 'LocateCell' as 'LocatorFindCell'. Added 'RGridFindCell'. 
//
//    Kathleen Bonnell, Wed Jun 18 18:07:55 PDT 2003   
//    Add member ghostType.
//
//    Kathleen Bonnell, Fri Jun 27 17:06:27 PDT 2003  
//    Add methods RetrieveNodes, RetrieveZones, DeterminePickedNode,
//    SetRealIds, RetrieveVarInfo -- to support Node Pick and make code
//    more legible.
//
//    Kathleen Bonnell, Thu Nov 13 09:05:14 PST 2003 
//    Remove RgridFindCell, LocatorFindCell -- functionality now resides
//    in vtkVisItUtility. 
//    
//    Hank Childs, Thu Feb  5 17:11:06 PST 2004
//    Moved inlined destructor definition to .C file because certain compilers
//    have problems with them.
//
//    Kathleen Bonnell, Mon Mar  8 15:39:15 PST 2004 
//    Added SetNeedTransform / SetTransform and setUseSet. 
//    
// ****************************************************************************

class QUERY_API avtPickQuery : public avtDatasetQuery
{
  public:
                                    avtPickQuery();
    virtual                        ~avtPickQuery();


    virtual const char             *GetType(void)
                                             { return "avtPickQuery"; };
    virtual const char             *GetDescription(void)
                                             { return "Picking"; };

    void                            SetPickAtts(const PickAttributes *pa);
    const PickAttributes *          GetPickAtts(void);
    void                            SetTransform(const avtMatrix *m);
    void                            SetNeedTransform(const bool b)
                                        { needTransform = b; };
    void                            SetUseSet(const unsignedCharVector &u)
                                        { useSet = u; };

  protected:
    PickAttributes                  pickAtts;
    int                             cellOrigin;
    int                             blockOrigin;
    avtGhostType                    ghostType;
    const avtMatrix                *invTransform;
    bool                            singleDomain;
    bool                            needTransform;
    avtQueryableSource             *src;
    unsignedCharVector              useSet;
    avtExpressionEvaluatorFilter   *eef;

    virtual void                    Execute(vtkDataSet *, const int);
    virtual void                    PreExecute(void);
    virtual void                    PostExecute(void);
    virtual avtDataObject_p         ApplyFilters(avtDataObject_p);   
    bool                            RetrieveNodes(vtkDataSet *, int);
    bool                            RetrieveZones(vtkDataSet *, int);
    bool                            DeterminePickedNode(vtkDataSet *, int &);
    void                            SetRealIds(vtkDataSet *);
    void                            RetrieveVarInfo(vtkDataSet *);
};


#endif


