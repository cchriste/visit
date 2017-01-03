/*****************************************************************************
*
* Copyright (c) 2000 - 2017, Lawrence Livermore National Security, LLC
* Produced at the Lawrence Livermore National Laboratory
* LLNL-CODE-442911
* All rights reserved.
*
* This file is  part of VisIt. For  details, see https://visit.llnl.gov/.  The
* full copyright notice is contained in the file COPYRIGHT located at the root
* of the VisIt distribution or at http://www.llnl.gov/visit/copyright.html.
*
* Redistribution  and  use  in  source  and  binary  forms,  with  or  without
* modification, are permitted provided that the following conditions are met:
*
*  - Redistributions of  source code must  retain the above  copyright notice,
*    this list of conditions and the disclaimer below.
*  - Redistributions in binary form must reproduce the above copyright notice,
*    this  list of  conditions  and  the  disclaimer (as noted below)  in  the
*    documentation and/or other materials provided with the distribution.
*  - Neither the name of  the LLNS/LLNL nor the names of  its contributors may
*    be used to endorse or promote products derived from this software without
*    specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT  HOLDERS AND CONTRIBUTORS "AS IS"
* AND ANY EXPRESS OR  IMPLIED WARRANTIES, INCLUDING,  BUT NOT  LIMITED TO, THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND  FITNESS FOR A PARTICULAR  PURPOSE
* ARE  DISCLAIMED. IN  NO EVENT  SHALL LAWRENCE  LIVERMORE NATIONAL  SECURITY,
* LLC, THE  U.S.  DEPARTMENT OF  ENERGY  OR  CONTRIBUTORS BE  LIABLE  FOR  ANY
* DIRECT,  INDIRECT,   INCIDENTAL,   SPECIAL,   EXEMPLARY,  OR   CONSEQUENTIAL
* DAMAGES (INCLUDING, BUT NOT  LIMITED TO, PROCUREMENT OF  SUBSTITUTE GOODS OR
* SERVICES; LOSS OF  USE, DATA, OR PROFITS; OR  BUSINESS INTERRUPTION) HOWEVER
* CAUSED  AND  ON  ANY  THEORY  OF  LIABILITY,  WHETHER  IN  CONTRACT,  STRICT
* LIABILITY, OR TORT  (INCLUDING NEGLIGENCE OR OTHERWISE)  ARISING IN ANY  WAY
* OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
* DAMAGE.
*
*****************************************************************************/


// ************************************************************************* //
//                            avtunvFileFormat.h                           //
// ************************************************************************* //

#ifndef AVT_unv_FILE_FORMAT_H
#define AVT_unv_FILE_FORMAT_H

#include <avtSTSDFileFormat.h>

#include <string>
#include <vector>
#include <set>
#include <vtk_zlib.h>

// Define my classes for the mesh
class UnvRange { // Element class
public:
  double trange[6] ;
} ;
class UnvElement { // Element class
public:
  int number;
  int label;
  int typelt ;
  mutable int matid ; // mutable
  int* nodes ;
  struct compare_UnvElement
  {
    bool operator () (const UnvElement& e1, const UnvElement& e2) const
    {
      return (e1.label < e2.label);
    };
  };
};

class UnvNode { // Node class
public:
  int number; // Numbered nodes
  int label; // Global node label
  double x,y,z ;
  mutable std::vector<int> nod2elts; // Reverse mesh connectivity
  struct compare_UnvNode
  {
    bool operator () (const UnvNode& n1, const UnvNode& n2) const
    {
      return (n1.label < n2.label);
    };
  };
};
class UnvFace { // Face class, can be a boundary face or internal one
public:
  int number;
  int element ;
  int facloc ;
  int matid ;
  double pressure ;
};
class UnvFacePressure { // a set of faces
public:
  int number;
  int label ;
  std::string name ;
  std::vector<UnvFace> faces;
};
// ****************************************************************************
//  Class: avtunvFileFormat
//
//  Purpose:
//      Reads in unv files as a plugin to VisIt.
//
//  Programmer: cessenat -- generated by xml2avt
//  Creation:   Thu Sep 27 17:38:20 PST 2012
//
// ****************************************************************************

class avtunvFileFormat : public avtSTSDFileFormat
{
  public:
                       avtunvFileFormat(const char *filename);
    virtual           ~avtunvFileFormat() {;};

    //
    // This is used to return unconvention data -- ranging from material
    // information to information about block connectivity.
    //
    virtual void      *GetAuxiliaryData(const char *var, const char *type,
                                      void *args, DestructorFunction &);
    //

    //
    // These are used to declare what the current time and cycle are for the
    // file.  These should only be defined if the file format knows what the
    // time and/or cycle is.
    //
    // virtual int       GetCycle(void);
    // virtual double    GetTime(void);
    //

    virtual const char    *GetType(void)   { return "unv"; };
    virtual void           FreeUpResources(void); 

    virtual vtkDataSet    *GetMesh(const char *);
    virtual vtkDataArray  *GetVar(const char *);
    virtual vtkDataArray  *GetVectorVar(const char *);
 
  protected:
    // DATA MEMBERS

    virtual void           PopulateDatabaseMetaData(avtDatabaseMetaData *);
    virtual void           ReadFile();
    virtual int is3DKnownElt(int ); // Provides 3D element value in connectivity nodefac array
    virtual int is2DKnownElt(int ); // Provides 2D element value in connectivity nodefac array
    virtual int is1DKnownElt(int ); // Provides 1D element value in connectivity nodefac array
    virtual int isKnownElt(int ); // Provides element value in counter array nbletsptyp
    virtual int getNbnodes(int ); // Provides 3D element number of nodes
    virtual int getNbfaces(int ); // Provides 3D element number of faces
    virtual int getEltDim(int ); // Provides element topological dimension

    virtual int getNeighbour3D (int iel, int l1, int l2, int l3) ; // Provides a neighbour for iel with 3 nodes.
    virtual int getNeighbour2D (int iel, int l1, int l2) ; // Provides a neighbour for iel with 2 nodes.

    virtual int getfastNeighbour3D (int iel, int l1, int l2, int l3, std::set<UnvNode, UnvNode::compare_UnvNode>::iterator * itrgs) ; // Provides a neighbour for iel with 3 nodes.
    virtual int getfastNeighbour2D (int iel, int l1, int l2, std::set<UnvNode, UnvNode::compare_UnvNode>::iterator * itrgs) ; // Provides a neighbour for iel with 2 nodes.

    virtual int getNbnolsv( ); // Provides the number of nodes to build the Face Pressure Load-Set nbnolsv
    virtual int getNbfaextv( ); // Provides the number of free faces nbfaextv
    virtual int getfastNbfaextv( ); // Provides the number of free faces nbfaextv
    virtual int getNbnodesFreeFaces( ); // Provides the number of nodes to build the Face Faces nbnff
    virtual void getNormal3D (float *one_entry, std::set<UnvElement, UnvElement::compare_UnvElement>::iterator itre, int iflo, int facloc); // Provides the normal to a face
    virtual void getNormal2D (float *one_entry, std::set<UnvElement, UnvElement::compare_UnvElement>::iterator itre, int facloc); // Provides the normal to a face

    virtual int getNbfreeSets(); // Gets the number of boundaries, i.e. free connected faces

    FILE* handle; // File handle for unv file
    gzFile gzhandle ; // File handle for unv.gz file
    int nbnodes ; // Total number of nodes in the mesh
    int maxnodl ; // Maximum node label in the mesh
    int nb3dmats ; // Highest material numbre for 3D elements
    int nb2dmats ; // Highest material numbre for 3D elements
    int nb1dmats ; // Highest material numbre for 3D elements
    int nb3dcells ; // Store the total number of volume cells
    int nb2dcells ; // Store the total number of surface cells
    int nb1dcells ; // Store the total number of surface cells
    std::string filename; // Mesh file name, including .unv or .unv.gz extension
    bool fileRead; // Says if file has already been read or not
    int debuglevel ; 
    double range[6] ; // geometrical mesh range
    int nbloadsets ; // Number of load-sets of 3D face pressure
    int nbfalsv ; // Number of known 3D elements in load-sets
    int nbnolsv ; // Number of different nodes in load-sets
    int revconnect ; // Reverse connectivity flag
    int nbfaextv ; // Number of free faces (faces on the boundary)
    int nbnff ; // Number of nodes on the free faces
    int nbletsptyp[7] ; // Number of elements per element type
    int nbfreesets ; // Number of free-sets, delimitating different boundaries
    int cdim ; // Mesh topological dimension.
    std::set<UnvNode, UnvNode::compare_UnvNode> meshUnvNodes; // Full mesh nodes
    std::set<UnvElement, UnvElement::compare_UnvElement> meshUnvElements;  // Full mesh all types of elements
    std::vector<UnvFacePressure> meshUnvFacePressures; // read face pressure load sets
    std::vector<UnvFace> freeUnvFaces; // List of free faces for 3D mesh.
};


#endif
