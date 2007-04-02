#ifndef __vtkLineLegend_h
#define __vtkLineLegend_h
#include <visit_vtk_exports.h>

#include <vtkActor2D.h>
#include <vtkTextMapper.h>

class vtkPolyData;
class vtkPolyDataMapper2D;
class vtkTransform;
class vtkTransformPolyDataFilter;


class VISIT_VTK_API vtkLineLegend : public vtkActor2D
{
public:
  vtkTypeMacro(vtkLineLegend,vtkActor2D);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Instantiate object. 
  static vtkLineLegend *New();
  
  // Description:
  // Access the Position instance variable. Reimplemented from base
  // class to ensure normalized viewport coordinates
  // This variable controls the lower left corner of the legend. 
  void SetPosition(float,float);
  void SetPosition(float x[2]);

  // Description:
  // Access the Position2 instance variable. This variable controls
  // the upper right corner of the legend. It is by default
  // relative to Position1 and in Normalized Viewport coordinates.
  void SetPosition2(float,float);
  void SetPosition2(float x[2]);
  vtkCoordinate *GetPosition2Coordinate();
  float *GetPosition2();
 
  // Description:
  // Draw the legend and annotation text to the screen.
  int RenderOpaqueGeometry(vtkViewport* viewport);
  int RenderTranslucentGeometry(vtkViewport*) { return 0; };
  virtual int RenderOverlay(vtkViewport* viewport);

  // Description:
  // Release any graphics resources that are being consumed by this actor.
  // The parameter window could be used to determine which graphic
  // resources to release.
  virtual void ReleaseGraphicsResources(vtkWindow *);

  // Description:
  // Enable/Disable bolding annotation text.
  vtkSetMacro(Bold, int);
  vtkGetMacro(Bold, int);
  vtkBooleanMacro(Bold, int);

  // Description:
  // Enable/Disable italicizing annotation text.
  vtkSetMacro(Italic, int);
  vtkGetMacro(Italic, int);
  vtkBooleanMacro(Italic, int);

  // Description:
  // Enable/Disable creating shadows on the annotation text. Shadows make 
  // the text easier to read.
  vtkSetMacro(Shadow, int);
  vtkGetMacro(Shadow, int);
  vtkBooleanMacro(Shadow, int);

  // Description:
  // Set/Get the font family for the annotation text. Three font types 
  // are available: Arial (VTK_ARIAL), Courier (VTK_COURIER), and 
  // Times (VTK_TIMES).
  vtkSetMacro(FontFamily, int);
  vtkGetMacro(FontFamily, int);
  void SetFontFamilyToArial() {this->SetFontFamily(VTK_ARIAL);};
  void SetFontFamilyToCourier() {this->SetFontFamily(VTK_COURIER);};
  void SetFontFamilyToTimes() {this->SetFontFamily(VTK_TIMES);};

  // Description:
  // Set/Get the font height for the annotation text.
  vtkSetClampMacro(FontHeight, float, 0, 0.2);
  vtkGetMacro(FontHeight, float);

  // Description:
  // Set/Get the title of the scalar bar actor,
  vtkSetStringMacro(Title);
  vtkGetStringMacro(Title);

  // Description:
  // Set/Get the visibility of the title annotation text. 
  vtkSetMacro(TitleVisibility, int);
  vtkGetMacro(TitleVisibility, int);
  vtkBooleanMacro(TitleVisibility, int);

  vtkTextProperty * GetTitleProperty() 
     { return this->TitleMapper->GetTextProperty();};

  vtkProperty2D * GetLineProperty() 
     { return this->LineActor->GetProperty();};

  // Description:
  // Set/Get the scalar bar width.
  vtkSetClampMacro(BarWidth,float, 0.0, 0.5);
  vtkGetMacro(BarWidth,float);

  // Shallow copy of a scalar bar actor. Overloads the virtual vtkProp method.
  void ShallowCopy(vtkProp *prop);

protected:
  vtkLineLegend();
  virtual ~vtkLineLegend();

  void BuildTitle(vtkViewport *);
  void BuildLine(vtkViewport *);

  int   Bold;
  int   Italic;
  int   Shadow;
  int   FontFamily;
  float FontHeight;

  vtkCoordinate *Position2Coordinate;

  float BarWidth;
  
  char          *Title;
  vtkTextMapper *TitleMapper;
  vtkActor2D    *TitleActor;
  int            TitleVisibility;
  int            TitleOkayToDraw;

  vtkTransform               *Transform;
  vtkTransformPolyDataFilter *TransformFilter;
  vtkPolyDataMapper2D        *LineMapper;
  vtkActor2D                 *LineActor;

  vtkTimeStamp  BuildTime;
  int LastSize[2];
  int LastOrigin[2];

private:
  vtkLineLegend(const vtkLineLegend&);
  void operator=(const vtkLineLegend&);
};


#endif

