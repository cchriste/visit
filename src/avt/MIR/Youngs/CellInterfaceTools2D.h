// Contributed by Thierry Carrard from
// Commissariat a l'Energie Atomique, (CEA)
// BP12, 91297 Arpajon, France


#ifndef __CellInterfaceTools2D_H
#define __CellInterfaceTools2D_H

#include "CellInterfaceCommon.h"

#ifndef REAL_COORD
#define REAL_COORD REAL3
#endif

FUNC_DECL
REAL makeTriangleSurfaceFunctions(
                                  const uchar3 triangle,
                                  const REAL_COORD* vertices,
                                  const REAL_COORD normal,
                                  REAL2 func[2]
                                  )
{

    // 1. load the data
    const REAL_COORD v0 = vertices[ triangle.x ];
    const REAL_COORD v1 = vertices[ triangle.y ];
    const REAL_COORD v2 = vertices[ triangle.z ];

    const REAL d0 = dot( v0 , normal );
    const REAL d1 = dot( v1 , normal );
    const REAL d2 = dot( v2 , normal );


    DBG_MESG("v0 = "<<v0.x<<','<<v0.y<<" d0="<<d0);
    DBG_MESG("v1 = "<<v1.x<<','<<v1.y<<" d1="<<d1);
    DBG_MESG("v2 = "<<v2.x<<','<<v2.y<<" d2="<<d2);


    // 2. compute

    // compute vector from point on v0-v2 that has distance d1 from Plane0
    REAL_COORD I = linearInterp( d0, v0, d2, v2 , d1 );
    DBG_MESG("I = "<<I.x<<','<<I.y);
    REAL_COORD vec = v1 - I;
    REAL length = sqrt( dot(vec,vec) );
    DBG_MESG("length = "<<length);

    // side length function = (x-d0) * length / (d1-d0) = (length/(d1-d0)) * x - length * d0 / (d1-d0)
    REAL2 linearFunc01 = make_REAL2( length/(d1-d0) , - length * d0 / (d1-d0) );
    // surface function = integral of distance function starting at d0
    func[0] = make_REAL2(0,0);
    if( d1 > d0 )
    {
        func[0]  = linearFunc01;
    }

    // side length function = (d2-x) * length / (d2-d1) = (-length/(d2-d1)) * x + d2*length / (d2-d1)
    REAL2 linearFunc12 = make_REAL2( -length/(d2-d1) , d2*length/(d2-d1) );
    // surface function = integral of distance function starting at d1
    func[1] = make_REAL2(0,0);
    if( d2 > d1 )
    {
        func[1] = linearFunc12;
    }

    return triangleSurf( v0, v1, v2 );
}

FUNC_DECL
REAL findTriangleSetCuttingPlane(
                                 const REAL_COORD normal,    // IN  , normal vector
                                 const REAL fraction,   // IN  , volume fraction
                                 const int nv,          // IN  , number of vertices
                                 const int nt,          // IN  , number of triangles
                                 const uchar3* tv,       // IN  , triangles connectivity, size=nt
                                 const REAL_COORD* vertices // IN  , vertex coordinates, size=nv
#ifdef __CUDACC__
                                 ,char* sdata           // TEMP Storage
#endif
                                 )
{
    ALLOC_LOCAL_ARRAY( derivatives, REAL2, nv-1 );
    ALLOC_LOCAL_ARRAY( index, unsigned char, nv );
    ALLOC_LOCAL_ARRAY( rindex, unsigned char, nv );

    // initialisation
    for(int i=0;i<nv;i++)
    {
        index[i] = i;
    }

    for(int i=0;i<(nv-1);i++)
    {
        derivatives[i] = make_REAL2(0,0);
    }

    // tri des sommets dans le sens de la normale
    sortVertices( nv, vertices, normal, index );

    // table d'indirection inverse
    for(int i=0;i<nv;i++)
    {
        rindex[ index[i] ] = i;
    }

    // surface totale
    REAL surface = 0;

    // construction de la fonction cubique par morceau du volume tronqu�
    for(int i=0;i<nt;i++)
    {
        // calcul de la surface de l'intersection plan/tetra aux point P1 et P2
        uchar3 triangle = sortTriangle( tv[i] , rindex );
        DBG_MESG( "\ntriangle "<<i<<" : "<<tv[i].x<<','<<tv[i].y<<','<<tv[i].z<<" -> "<<triangle.x<<','<<triangle.y<<','<<triangle.z );

        // calcul des sous fonctions cubiques du volume derriere le plan en fonction de la distance
        REAL2 triangleSurfFunc[2];
        surface += makeTriangleSurfaceFunctions( triangle, vertices, normal, triangleSurfFunc );

#ifdef DEBUG
        for(int k=0;k<2;k++)
        {
            DBG_MESG( "surf'["<<k<<"] = "<<triangleSurfFunc[k].x<<','<<triangleSurfFunc[k].y );
        }
#endif

        // surface function bounds
        unsigned int i0 = rindex[ triangle.x ];
        unsigned int i1 = rindex[ triangle.y ];
        unsigned int i2 = rindex[ triangle.z ];

        DBG_MESG( "surf(x) steps = "<<i0<<','<<i1<<','<<i2 );

        DBG_MESG( "ajout surfFunc sur ["<<i0<<';'<<i1<<"]" );
        for(unsigned int j=i0;j<i1;j++)
        {
            derivatives[j] += triangleSurfFunc[0];
        }

        DBG_MESG( "ajout surfFunc sur ["<<i1<<';'<<i2<<"]" );
        for(unsigned int j=i1;j<i2;j++)
        {
            derivatives[j] += triangleSurfFunc[1];
        }
    }

    // calcul de la fraction de surface recherchee
    REAL y = surface*fraction;
    DBG_MESG( "surface = "<<surface<<", surface*fraction = "<<y );

    // integration des fonctions et recherche de la fonction utile
    REAL sum = 0;
    REAL3 surfaceFunction = make_REAL3(0,0,0);
    REAL xmin = 0;
    REAL xmax = dot( vertices[index[0]], normal ) ;
    int s = -1;
    while( sum<y && s<(nv-2) )
    {
        xmin = xmax;
        y -= sum;
        ++ s;
        REAL3 F = integratePolynomialFunc( derivatives[s] );
        F.z = - evalPolynomialFunc( F , xmin );
        surfaceFunction = F;
        xmax = dot( vertices[index[s+1]] , normal );
        sum = evalPolynomialFunc( F, xmax );
    }
    if( s<0) s=0;

    // calcul du volume recherche

    // recherche de la portion de fonction qui contient la valeur
    DBG_MESG( "step="<<s<<", x in ["<<xmin<<';'<<xmax<<']' );

    /* chaque portion de fonction redemarre de 0,
       on calcul donc le volume recherch� dans cette portion de fonction
    */
    DBG_MESG( "surface reminder = "<< y );

    // recherche par newton
    //REAL x = quadraticFunctionSolve( funcs[s], surface, xmin, xmax );
    REAL x = newtonSearchPolynomialFunc( surfaceFunction, derivatives[s], y, xmin, xmax );

    DBG_MESG( "final x = "<< x );

    FREE_LOCAL_ARRAY( derivatives, REAL2, nv-1 );
    FREE_LOCAL_ARRAY( index, unsigned char, nv );
    FREE_LOCAL_ARRAY( rindex, unsigned char, nv );

    return x ;
}

#endif
