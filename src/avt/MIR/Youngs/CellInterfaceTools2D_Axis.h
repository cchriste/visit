// Contributed by Thierry Carrard from
// Commissariat a l'Energie Atomique, (CEA)
// BP12, 91297 Arpajon, France


#ifndef __CellInterfaceTools2D_Axis_H
#define __CellInterfaceTools2D_Axis_H

#include "CellInterfaceCommon.h"


/*
  compute the derivatives of the piecewise cubic function of the volume behind the cutting cone ( axis symetric 2D plane)
*/
FUNC_DECL
void makeConeVolumeDerivatives(
    const uchar3 triangle,
    const REAL2* vertices,
    const REAL2 normal,
    REAL3 deriv[2]
    )
{

    // 1. load the data
    const REAL2 v0 = vertices[ triangle.x ];
    const REAL2 v1 = vertices[ triangle.y ];
    const REAL2 v2 = vertices[ triangle.z ];

    // 2. compute
    const REAL d0 = dot( v0 , normal );
    const REAL d1 = dot( v1 , normal );
    const REAL d2 = dot( v2 , normal );

    DBG_MESG("v0 = "<<v0.x<<','<<v0.y<<" d0="<<d0);
    DBG_MESG("v1 = "<<v1.x<<','<<v1.y<<" d1="<<d1);
    DBG_MESG("v2 = "<<v2.x<<','<<v2.y<<" d2="<<d2);

    // compute vector from point on v0-v2 that has distance d1 from Plane0
    REAL2 I = linearInterp( d0, v0, d2, v2 , d1 );
    DBG_MESG("I = "<<I.x<<','<<I.y);
    REAL2 vec = v1 - I;
    REAL length = sqrt( dot(vec,vec) );
    DBG_MESG("length = "<<length);

    // compute truncated cone surface at d1
    REAL Isurf = REAL_CONST(M_PI) * FABS(I.y+v1.y) * length; // 2 * REAL_CONST(M_PI) * ( (I.y+v1.y) * 0.5 ) * length ;
    REAL coef;

    // build cubic volume functions derivatives
    coef = ( d1 > d0 )  ?  ( Isurf / ((d1-d0)*(d1-d0)) ) : REAL_CONST(0.0) ;
    deriv[0] = coef * make_REAL3( 1 , -2*d0 , d0*d0 ) ;

    coef = ( d2 > d1 )  ?  ( Isurf / ((d2-d1)*(d2-d1)) ) : REAL_CONST(0.0) ;
    deriv[1] = coef * make_REAL3( 1 , -2*d2 , d2*d2 ) ;
}


FUNC_DECL
REAL findTriangleSetCuttingCone(
    const REAL2 normal,    // IN  , normal vector
    const REAL fraction,   // IN  , volume fraction
    const int nv,          // IN  , number of vertices
    const int nt,          // IN  , number of triangles
    const uchar3* tv,       // IN  , triangles connectivity, size=nt
    const REAL2* vertices // IN  , vertex coordinates, size=nv
#ifdef __CUDACC__
    ,char* sdata           // TEMP Storage
#endif
    )
{
    ALLOC_LOCAL_ARRAY( derivatives, REAL3, nv-1 );
    ALLOC_LOCAL_ARRAY( index, unsigned char, nv );
    ALLOC_LOCAL_ARRAY( rindex, unsigned char, nv );

    // initialisation
    for(int i=0;i<nv;i++)
    {
        index[i] = i;
    }

    for(int i=0;i<(nv-1);i++)
    {
        derivatives[i] = make_REAL3(0,0,0);
    }

    // tri des sommets dans le sens de la normale
    sortVertices( nv, vertices, normal, index );

    // table d'indirection inverse
    for(int i=0;i<nv;i++)
    {
        rindex[ index[i] ] = i;
    }

    // construction de la fonction cubique par morceau du volume tronqu�
    for(int i=0;i<nt;i++)
    {
        // calcul de la surface de l'intersection plan/tetra aux point P1 et P2
        uchar3 triangle = sortTriangle( tv[i] , rindex );
        DBG_MESG( "\ntriangle "<<i<<" : "<<tv[i].x<<','<<tv[i].y<<','<<tv[i].z<<" -> "<<triangle.x<<','<<triangle.y<<','<<triangle.z );

        // calcul des sous fonctions cubiques du volume derriere le plan en fonction de la distance
        REAL3 coneVolDeriv[2];
        makeConeVolumeDerivatives( triangle, vertices, normal, coneVolDeriv );

        // surface function bounds
        unsigned int i0 = rindex[ triangle.x ];
        unsigned int i1 = rindex[ triangle.y ];
        unsigned int i2 = rindex[ triangle.z ];

        DBG_MESG( "surf(x) steps = "<<i0<<','<<i1<<','<<i2 );

        DBG_MESG( "ajout surfFunc sur ["<<i0<<';'<<i1<<"]" );
        for(unsigned int j=i0;j<i1;j++)
        {
            derivatives[j] += coneVolDeriv[0];
        }

        DBG_MESG( "ajout surfFunc sur ["<<i1<<';'<<i2<<"]" );
        for(unsigned int j=i1;j<i2;j++)
        {
            derivatives[j] += coneVolDeriv[1];
        }
    }

    REAL surface = 0;
    REAL xmin = 0;
    REAL xmax = dot( vertices[index[0]], normal ) ;
    for(int i=0;i<(nv-1);i++)
    {
        xmin = xmax;
        REAL4 F = integratePolynomialFunc( derivatives[i] );
        F.w = - evalPolynomialFunc( F , xmin );
        xmax = dot( vertices[index[i+1]] , normal );
        surface += evalPolynomialFunc( F, xmax );
    }

    REAL y = surface*fraction;
    DBG_MESG( "surface = "<<surface<<", surface*fraction = "<<y );

    // integration des fonctions de surface en fonctions de volume
    REAL sum = 0;
    REAL4 volumeFunction = make_REAL4(0,0,0,0);
    xmax = dot( vertices[index[0]], normal ) ;
    int s = -1;
    while( sum<y && s<(nv-2) )
    {
        xmin = xmax;
        y -= sum;
        ++ s;
        REAL4 F = integratePolynomialFunc( derivatives[s] );
        F.w = - evalPolynomialFunc( F , xmin );
        volumeFunction = F;
        xmax = dot( vertices[index[s+1]] , normal );
        sum = evalPolynomialFunc( F, xmax );
    }
    if( s<0) s=0;

    // recherche de la portion de fonction qui contient la valeur
    DBG_MESG( "step="<<s<<", x in ["<<xmin<<';'<<xmax<<']' );

    /* chaque portion de fonction redemarre de 0,
       on calcul donc le volume recherch� dans cette portion de fonction
    */
    DBG_MESG( "surface reminder = "<< y );

    // recherche par newton
    //REAL x = quadraticFunctionSolve( funcs[s], surface, xmin, xmax );
    REAL x = newtonSearchPolynomialFunc( volumeFunction, derivatives[s], y, xmin, xmax );

    DBG_MESG( "final x = "<< x );

    FREE_LOCAL_ARRAY( derivatives, REAL3        , nv-1 );
    FREE_LOCAL_ARRAY( index      , unsigned char, nv   );
    FREE_LOCAL_ARRAY( rindex     , unsigned char, nv   );

    return x ;
}

#endif
