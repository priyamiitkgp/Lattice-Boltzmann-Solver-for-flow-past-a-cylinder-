#include <iostream>
#include <fstream>
#include <stdio.h>
#include <iomanip>
#include<cmath>
#include <stdlib.h>     /* abort, NULL */
#include <math.h>
enum direction{
    ZERO_ZERO,
    PLUS_ZERO,
    ZERO_PLUS,
    MINUS_ZERO,
    ZERO_MINUS,
    PLUS_PLUS,
    MINUS_PLUS,
    MINUS_MINUS,
    PLUS_MINUS};

#define N_DV 9
#define SIZE_X 2500
#define SIZE_Y 250
#define SIZE_Z 1
const int r=5;

//define structure for the D2Q9 model
template <typename dataType>
struct ModelD2Q9
{
    dataType ci_x[N_DV];
    dataType ci_y[N_DV];
    dataType dVopposite[N_DV];
    dataType cc[N_DV];
    dataType cy2[N_DV];
    dataType cx2[N_DV];
    dataType wt[N_DV];
    dataType theta0,oneBytheta0,sqrtTheta0;
    dataType c2[N_DV];
    dataType rho;
    dataType u1,u2,theta;
    dataType inlet;
    dataType inletActual;
    dataType Ma;
    dataType fEQ[N_DV];
    dataType denominator;
};

//define grid class
class Grid
{
public:
//       double data[ (SIZE_X+2)*(SIZE_Y+2)*N_DV ];

    double *data;
    Grid();

    double  operator()(const int i,const int j, const int dv) const  { return data[  j*(SIZE_X+2)*N_DV  +  i*N_DV +  dv ];}
    double& operator()(const int i,const int j, const int dv)        { return data[  j*(SIZE_X+2)*N_DV  +  i*N_DV +  dv ];}

    void initializeGrid()
    {
        for (int i = 0; i < (SIZE_X + 2) * (SIZE_Y + 2) * (N_DV); i++)
            data[i] = i;

    }

};
//constructor
Grid::Grid()
{
    data = new double [(SIZE_X+2)*(SIZE_Y+2)*N_DV ];
}

void advection(Grid &gridLB)
{
    for(int j=1;j<=SIZE_Y;j++)
    {
        for(int i=1;i<=SIZE_X;i++)
        {//ZERO_PLUS
            gridLB(i, j, MINUS_ZERO) = gridLB(i+1, j  , MINUS_ZERO);
            //MINUS_ZERO
            gridLB(i, j, ZERO_MINUS) = gridLB(i  , j+1, ZERO_MINUS);
            //PLUS_PLUS
            gridLB(i, j, MINUS_MINUS)= gridLB(i+1, j+1, MINUS_MINUS);
            //PLUS_MINUS
            gridLB(i, j, PLUS_MINUS) = gridLB(i-1, j+1, PLUS_MINUS);
        }
    }

    for(int j=SIZE_Y;j>=1;j--)
    {
        for(int i=SIZE_X;i>=1;i--)
        {//PLUS_ZERO
            gridLB(i, j, PLUS_ZERO) = gridLB(i-1, j  , PLUS_ZERO);
            //ZERO_MINUS
            gridLB(i, j, ZERO_PLUS) = gridLB(i  , j-1, ZERO_PLUS);
            //MINUS_MINUS
            gridLB(i, j, PLUS_PLUS) = gridLB(i-1, j-1, PLUS_PLUS);
            //MINUS_PLUS
            gridLB(i, j, MINUS_PLUS)= gridLB(i+1, j-1, MINUS_PLUS);
        }
    }
}

void getfEQ(ModelD2Q9<double> &myModel)
{
    double uDotC, uSq, tempX, tempY, term1, term2, term3, term4, term5, oneByTerm3, oneByTerm4 ;

    for(int i=0;i<N_DV;i++)
    {
        uDotC = myModel.u1*myModel.ci_x[i] + myModel.u2*myModel.ci_y[i];
        uSq   = myModel.u2*myModel.u2+myModel.u1*myModel.u1 ;

        myModel.fEQ[i]= myModel.wt[i]*myModel.rho*(1.0 + uDotC*myModel.oneBytheta0 - 1.5*uSq + 4.5*uDotC*uDotC );

// 	 myModel.fEQ[i]=myModel.wt[i]*myModel.rho*(1.0 + uDotC*myModel.oneBytheta0 -1.5*(myModel.u2*myModel.u2+myModel.u1*myModel.u1) +4.5*uDotC*uDotC  );
// 	 myModel.fEQ[i]= myModel.wt[i]*myModel.rho*(1.0 + uDotC*myModel.oneBytheta0 - 1.5*uSq + 4.5*uDotC*uDotC  + 4.5*uDotC*uDotC*uDotC -4.5*uSq*uDotC );
    }

//   tempX = sqrt(1+3.0*myModel.u1*myModel.u1) ;
//   tempY = sqrt(1+3.0*myModel.u2*myModel.u2) ;
//
//   term1 = 2.0 - tempX;
//   term2 = 2.0 - tempY;
//   term5 = term1*term2;
//
//   term3 = (2.0*myModel.u1 + tempX)/(1.0-myModel.u1);
//   term4 = (2.0*myModel.u2 + tempY)/(1.0-myModel.u2);
//   oneByTerm3 = 1.0/term3;
//   oneByTerm4 = 1.0/term4;
//                                                                            ;
//   myModel.fEQ[ZERO_ZERO  ] = myModel.rho*myModel.wt[ZERO_ZERO  ]*term5             ;
//   myModel.fEQ[PLUS_ZERO  ] = myModel.rho*myModel.wt[PLUS_ZERO  ]*term5*term3       ;
//   myModel.fEQ[ZERO_PLUS  ] = myModel.rho*myModel.wt[ZERO_PLUS  ]*term5*term4       ;
//   myModel.fEQ[MINUS_ZERO ] = myModel.rho*myModel.wt[MINUS_ZERO ]*term5*oneByTerm3       ;
//   myModel.fEQ[ZERO_MINUS ] = myModel.rho*myModel.wt[ZERO_MINUS ]*term5*oneByTerm4       ;
//   myModel.fEQ[PLUS_PLUS  ] = myModel.rho*myModel.wt[PLUS_PLUS  ]*term5*term3*term4 ;
//   myModel.fEQ[MINUS_PLUS ] = myModel.rho*myModel.wt[MINUS_PLUS ]*term5*oneByTerm3*term4 ;
//   myModel.fEQ[MINUS_MINUS] = myModel.rho*myModel.wt[MINUS_MINUS]*term5*oneByTerm3*oneByTerm4 ;
//   myModel.fEQ[PLUS_MINUS ] = myModel.rho*myModel.wt[PLUS_MINUS ]*term5*term3*oneByTerm4 ;

    // Exact in temperature

//   double dot, K, ux2(myModel.u1*myModel.u1), uy2(myModel.u2*myModel.u2), uSq(ux2 + uy2), theta0(1.0/3.0), oneByTheta0(3.0);
//   double eta(myModel.theta/theta0-1.0);
//   double theta2(myModel.theta*myModel.theta);
//   double oneByTheta(1.0/myModel.theta);
//   double oneByTheta2(oneByTheta*oneByTheta) ;
//
//   myModel.fEQ[0] = myModel.rho*(1.0-myModel.theta)*(1.0-myModel.theta);
//   for (int dv = 1; dv < 5; dv++)
//     myModel.fEQ[dv] = myModel.rho*(1.0-myModel.theta)*0.5*myModel.theta;
//   for (int dv = 5; dv < 9; dv++)
//     myModel.fEQ[dv] = myModel.rho*0.25*myModel.theta*myModel.theta;
//
//   for (int dv = 0; dv < 9; dv++)
//   {
//     dot = (myModel.u1*myModel.ci_x[dv] + myModel.u2*myModel.ci_y[dv])*oneByTheta;
//     K = 2.0*myModel.theta*myModel.theta/(1.0-myModel.theta) + 0.5*myModel.cc[dv]*(1.0-3.0*myModel.theta)/(1.0-myModel.theta);
//     myModel.fEQ[dv] = myModel.wt[dv]*(1.0 + dot + 0.5*dot*dot - 0.5*uSq*oneByTheta2*K );
//   }
//


}

void periodicAll(ModelD2Q9<double> &myModel,Grid &gridLB)
{
    for(int dv=0;dv<9;dv++)
    {
        for(int j=0;j<=SIZE_Y+1;j++)
        {
            //copy from right wall to left ghost
            gridLB(0, j, dv)   = gridLB(SIZE_X, j, dv);
            //copy from left wall to right ghost
            gridLB(SIZE_X+1, j, dv)    = gridLB(1, j, dv);
        }
        for(int i=0;i<=SIZE_X+1;i++)
        {
            //copy from bottom wall to top ghost
            gridLB(i, SIZE_Y+1, dv)=gridLB(i, 1, dv);
            //copy from bottom wall to top ghost
            gridLB(i, 0, dv)=gridLB(i, SIZE_Y, dv);
        }
    }
}

void prepareWallTopBottom(ModelD2Q9<double> &myModel,Grid &gridLB)
{


    for(int dv=0;dv<N_DV;dv++)
    {
        for(int i=0;i <= SIZE_X+1;i++)
        {  //copy from bottom wall to bottom ghost
            gridLB(i, 0, dv)=gridLB(i, 1, dv);
            //copy from top wall at to top ghost
            gridLB(i, SIZE_Y+1, dv)=gridLB(i, SIZE_Y, dv);
        }

    }
}
void getMoments(ModelD2Q9<double> &myModel,Grid &gridLB,int x_coord,int y_coord)
{
    myModel.rho = 0.0;
    myModel.u1  = 0.0;
    myModel.u2  = 0.0;
    myModel.theta= 0.0;
    for(int k=0;k<N_DV;k++)
    {
        myModel.rho   += gridLB(x_coord, y_coord, k);
        myModel.u1    += gridLB(x_coord, y_coord, k) * myModel.ci_x[k];
        myModel.u2    += gridLB(x_coord, y_coord, k) * myModel.ci_y[k];
        myModel.theta += gridLB(x_coord, y_coord, k) * myModel.cc[k];
    }

    myModel.u1 /= myModel.rho;
    myModel.u2 /= myModel.rho;
//     if(x_coord>=x && x_coord<x+10&&y_coord>=y && y_coord<y+10){
//         myModel.u1 = 0.0;
//         myModel.u2 = 0.0;
//     }
// 	myModel.theta /= myModel.rho;
// 	myModel.theta -= (myModel.u1*myModel.u1 + myModel.u2*myModel.u2);
// 	myModel.theta *= 0.5;

    myModel.theta = myModel.theta0;

}
void prepareWallInlet(ModelD2Q9<double> &myModel,Grid &gridLB)
{


    for(int dv=0;dv<N_DV;dv++)
    {
        for(int j=1;j <= SIZE_Y;j++)
        {
            //copy from left wall to left ghost
           // getMoments(myModel,gridLB,2,j);
           // gridLB(2,j,0) = gridLB(2,j,0)-(myModel.rho - 1.0);
            gridLB(0, j, dv)  = gridLB(1, j,  dv);
        }

    }
//    for (int j=1;j <= SIZE_Y;j++){
//        getMoments(myModel,gridLB,0,j);
//        gridLB(0,j,0) = gridLB(0,j,0)-(myModel.rho - 1.0);
//    }
}

void prepareWallOutlet(ModelD2Q9<double> &myModel,Grid &gridLB, int i)
{
    for(int dv=0;dv<N_DV;dv++)
    {
        for(int j=1;j <= SIZE_Y;j++)
        { //copy from right wall to right ghost
            gridLB(SIZE_X+1, j, dv)   = gridLB(SIZE_X, j, dv);
        }
    }
}

void prepareWallObject(ModelD2Q9<double> &myModel,Grid &gridLB, double fLeft[][9], double fRight[][9], double fTop[][9], double fBottom[][9], int objectOriginX,int objectOriginY,int objectLengthX,int objectLengthY)
{

    int objectEndX = objectOriginX + objectLengthX - 1;
    int objectEndY = objectOriginY + objectLengthY - 1;

    for(int i2 = objectOriginY-1; i2 <= objectEndY+1; i2++)
        for(int dv=0;dv<N_DV;dv++)
        {
            fLeft[i2-objectOriginY+1][dv] = gridLB(objectOriginX-1,i2,dv);
            fRight[i2-objectOriginY+1][dv] = gridLB(objectEndX+1   ,i2,dv);
        }

    for(int i1 = objectOriginX-1; i1 <= objectEndX+1; i1++)
        for(int dv=0;dv<N_DV;dv++)
        {
            fBottom[i1-objectOriginX+1][dv] = gridLB(i1,objectOriginY-1,dv);
            fTop[i1-objectOriginX+1][dv] = gridLB(i1,objectEndY+1   ,dv);
        }

}

void prepareWallCylinder(ModelD2Q9<double> &myModel,Grid &gridLB,double fwall[][2*r+3][N_DV], int x,int y,int r)
{


    for (int i=0;i<=2*r+2;i++) {
        for (int j = 0; j <= 2 * r + 2; j++) {
            for (int dv = 0; dv < N_DV; dv++) {
                fwall[i][j][dv] = 0.0;
            }
        }
    }
    int xCentre, yCentre;
    xCentre = r+1;
    yCentre = r+1;

    for (int i=0; i<=2*r+2; i++) {
        for (int j=0; j <= 2*r + 2; j++) {
            if( (i-xCentre)*(i-xCentre) + (j-yCentre)*(j-yCentre) <= r*r)
            {
                for (int dv=0; dv < N_DV; dv++) {

                    if ( (i-xCentre-myModel.ci_x[dv])*(i-xCentre-myModel.ci_x[dv]) + (j-yCentre-myModel.ci_y[dv])*(j-yCentre-myModel.ci_y[dv]) > r*r) {
                        fwall[i-(int)myModel.ci_x[dv]][j-(int)myModel.ci_y[dv]][(int)myModel.dVopposite[dv]] = gridLB( x+i-myModel.ci_x[dv], y+j-myModel.ci_y[dv], dv);
                    }

                }


            }

        }}


}

void applyWallBC_old(ModelD2Q9<double> &myModel,Grid &gridLB)
{
    myModel.rho=1.0;
    myModel.u1=myModel.inlet;
    myModel.u2=0.0;

    getfEQ(myModel);
    myModel.denominator=1.0/(myModel.fEQ[ZERO_MINUS]+myModel.fEQ[PLUS_MINUS]+myModel.fEQ[MINUS_MINUS]);

    double rhoTemp;
    double factor,denom;

    for(int j=2;j <= SIZE_Y-1;j++)
    {
        //bounceback from right wall
        factor =  gridLB(SIZE_X+1, j, PLUS_ZERO) + gridLB(SIZE_X+1, j, PLUS_PLUS) + gridLB(SIZE_X+1, j, PLUS_MINUS)   ;
        denom  =  1.0/(myModel.wt[PLUS_ZERO] + myModel.wt[PLUS_PLUS] + myModel.wt[PLUS_MINUS] );
        gridLB(SIZE_X, j, MINUS_ZERO ) = factor*denom*myModel.wt[MINUS_ZERO ];
        gridLB(SIZE_X, j, MINUS_MINUS) = factor*denom*myModel.wt[MINUS_MINUS];
        gridLB(SIZE_X, j, MINUS_PLUS ) = factor*denom*myModel.wt[MINUS_PLUS ];

        //bounceback from left wall
        factor =  gridLB(0, j, MINUS_ZERO) + gridLB(0, j, MINUS_PLUS) + gridLB(0, j, MINUS_MINUS)   ;
        denom  =  1.0/(myModel.wt[MINUS_ZERO] + myModel.wt[MINUS_PLUS] + myModel.wt[MINUS_MINUS] );
        gridLB(1, j, PLUS_ZERO ) = factor*denom*myModel.wt[PLUS_ZERO ];
        gridLB(1, j, PLUS_MINUS) = factor*denom*myModel.wt[PLUS_MINUS];
        gridLB(1, j, PLUS_PLUS ) = factor*denom*myModel.wt[PLUS_PLUS ];
    }

    for(int i=2;i<=SIZE_X-1;i++)
    {
        //classical diffuse from bottom wall
        factor =  gridLB(i, 0, ZERO_MINUS) + gridLB(i, 0, PLUS_MINUS) + gridLB(i, 0, MINUS_MINUS)   ;
        denom  =  1.0/(myModel.wt[ZERO_MINUS] + myModel.wt[PLUS_MINUS] + myModel.wt[MINUS_MINUS] );
        gridLB(i, 1, ZERO_PLUS ) = factor*denom*myModel.wt[ZERO_PLUS ] ;
        gridLB(i, 1, MINUS_PLUS) = factor*denom*myModel.wt[MINUS_PLUS] ;
        gridLB(i, 1, PLUS_PLUS ) = factor*denom*myModel.wt[PLUS_PLUS ] ;

        //classical diffuse top wall
        factor =  gridLB(i, SIZE_Y+1, ZERO_PLUS) + gridLB(i, SIZE_Y+1, PLUS_PLUS) + gridLB(i, SIZE_Y+1, MINUS_PLUS)   ;
        gridLB(i, SIZE_Y, ZERO_MINUS) = factor*myModel.denominator*myModel.fEQ[ZERO_MINUS ];
        gridLB(i, SIZE_Y, PLUS_MINUS) = factor*myModel.denominator*myModel.fEQ[PLUS_MINUS ];
        gridLB(i, SIZE_Y, MINUS_MINUS)= factor*myModel.denominator*myModel.fEQ[MINUS_MINUS];
    }

    //BOTTOM LEFT CORNER
    gridLB(1 , 1, ZERO_PLUS )=  gridLB(1 , 0, ZERO_MINUS) ;
    gridLB(1 , 1, PLUS_PLUS )=  gridLB(1 , 0, MINUS_MINUS);
    gridLB(1 , 1, PLUS_ZERO )=  gridLB(1 , 0, MINUS_ZERO) ;
    gridLB(1 , 1, PLUS_MINUS)=  gridLB(1 , 0, MINUS_PLUS) ;
    gridLB(1 , 1, MINUS_PLUS)=  gridLB(1 , 0, PLUS_MINUS) ;
    //BOTTOM RIGHT CORNER
    gridLB(SIZE_X, 1 , ZERO_PLUS  ) = gridLB(SIZE_X, 0 , ZERO_MINUS) ;
    gridLB(SIZE_X, 1 , MINUS_PLUS ) = gridLB(SIZE_X, 0 , PLUS_MINUS) ;
    gridLB(SIZE_X, 1 , MINUS_ZERO ) = gridLB(SIZE_X, 0 , PLUS_ZERO) ;
    gridLB(SIZE_X, 1 , PLUS_PLUS  ) = gridLB(SIZE_X, 0 , MINUS_MINUS);
    gridLB(SIZE_X, 1 , MINUS_MINUS )= gridLB(SIZE_X, 0 , PLUS_PLUS) ;
    //TOP LEFT CORNER
    gridLB(1 , SIZE_Y, PLUS_PLUS  ) =  gridLB(1 , SIZE_Y+1, MINUS_MINUS) ;
    gridLB(1 , SIZE_Y, MINUS_MINUS) =  gridLB(1 , SIZE_Y+1, PLUS_PLUS);
    gridLB(1 , SIZE_Y, PLUS_ZERO  ) =  gridLB(1 , SIZE_Y+1, MINUS_ZERO )  ;
    gridLB(1 , SIZE_Y, PLUS_MINUS ) =  gridLB(1 , SIZE_Y+1, MINUS_PLUS)  ;
    gridLB(1 , SIZE_Y, ZERO_MINUS ) =  gridLB(1 , SIZE_Y+1, ZERO_PLUS)  ;

    rhoTemp = 0.0;
    for(int dv=0;dv<N_DV;dv++)
        rhoTemp += gridLB(1, SIZE_Y, dv);
    for(int dv=0;dv<N_DV;dv++)
        gridLB(1, SIZE_Y, dv) = rhoTemp*myModel.fEQ[dv];

    //TOP RIGHT CORNER
    gridLB(SIZE_X , SIZE_Y, PLUS_MINUS ) =  gridLB(SIZE_X , SIZE_Y+1, MINUS_PLUS) ;
    gridLB(SIZE_X , SIZE_Y, MINUS_PLUS ) =  gridLB(SIZE_X , SIZE_Y+1, PLUS_MINUS);
    gridLB(SIZE_X , SIZE_Y, MINUS_ZERO ) =  gridLB(SIZE_X , SIZE_Y+1, PLUS_ZERO ) ;
    gridLB(SIZE_X , SIZE_Y, MINUS_MINUS) =  gridLB(SIZE_X , SIZE_Y+1, PLUS_PLUS) ;
    gridLB(SIZE_X , SIZE_Y, ZERO_MINUS ) =  gridLB(SIZE_X , SIZE_Y+1, ZERO_PLUS) ;

    rhoTemp = 0.0;
    for(int dv=0;dv<N_DV;dv++)
        rhoTemp += gridLB(SIZE_X, SIZE_Y, dv);
    for(int dv=0;dv<N_DV;dv++)
        gridLB(SIZE_X, SIZE_Y, dv) = rhoTemp*myModel.fEQ[dv];


}


//void getMoments(ModelD2Q9<double> &myModel,Grid &gridLB,int x_coord,int y_coord)
//{
//    myModel.rho = 0.0;
//    myModel.u1  = 0.0;
//    myModel.u2  = 0.0;
//    myModel.theta= 0.0;
//    for(int k=0;k<N_DV;k++)
//    {
//        myModel.rho   += gridLB(x_coord, y_coord, k);
//        myModel.u1    += gridLB(x_coord, y_coord, k) * myModel.ci_x[k];
//        myModel.u2    += gridLB(x_coord, y_coord, k) * myModel.ci_y[k];
//        myModel.theta += gridLB(x_coord, y_coord, k) * myModel.cc[k];
//    }
//
//    myModel.u1 /= myModel.rho;
//    myModel.u2 /= myModel.rho;
////     if(x_coord>=x && x_coord<x+10&&y_coord>=y && y_coord<y+10){
////         myModel.u1 = 0.0;
////         myModel.u2 = 0.0;
////     }
//// 	myModel.theta /= myModel.rho;
//// 	myModel.theta -= (myModel.u1*myModel.u1 + myModel.u2*myModel.u2);
//// 	myModel.theta *= 0.5;
//
//    myModel.theta = myModel.theta0;
//
//}

void applyWallOutletFeq(ModelD2Q9<double> &myModel,Grid &gridLB)
{

    double rhoTemp;
    double factor,denom;

    for(int j=1;j <= SIZE_Y;j++)
    {
        //OUTLET-------RIGHT WALL
        getMoments(myModel,gridLB,SIZE_X+1,j);
        getfEQ(myModel);
//           gridLB(SIZE_X,j,MINUS_ZERO)=myModel.fEQ[MINUS_ZERO];
//           gridLB(SIZE_X,j,MINUS_MINUS)=myModel.fEQ[MINUS_MINUS];
//           gridLB(SIZE_X,j,MINUS_PLUS)=myModel.fEQ[MINUS_PLUS];
        for (int dv=0;dv<N_DV;dv++){
            gridLB(SIZE_X,j,dv)=myModel.fEQ[dv];
        }

    }

}

void applyWallObject(ModelD2Q9<double> &myModel,Grid &gridLB, double fLeft[][9], double fRight[][9], double fTop[][9], double fBottom[][9], int objectOriginX,int objectOriginY,int objectLengthX,int objectLengthY)
{


    int objectEndX = objectOriginX + objectLengthX - 1;
    int objectEndY = objectOriginY + objectLengthY - 1;

    for(int i2 = objectOriginY; i2 <= objectEndY; i2++)
    {
        gridLB(objectOriginX-1,i2,MINUS_PLUS ) =  fLeft[i2-objectOriginY+1][PLUS_MINUS ] ;
        gridLB(objectOriginX-1,i2,MINUS_ZERO ) =  fLeft[i2-objectOriginY+1][PLUS_ZERO  ] ;
        gridLB(objectOriginX-1,i2,MINUS_MINUS) =  fLeft[i2-objectOriginY+1][PLUS_PLUS  ] ;

        gridLB(objectEndX+1   ,i2,PLUS_ZERO  ) = fRight[i2-objectOriginY+1][MINUS_ZERO ] ;
        gridLB(objectEndX+1   ,i2,PLUS_MINUS ) = fRight[i2-objectOriginY+1][MINUS_PLUS ] ;
        gridLB(objectEndX+1   ,i2,PLUS_PLUS  ) = fRight[i2-objectOriginY+1][MINUS_MINUS] ;
    }

    for(int i1 = objectOriginX; i1 <= objectEndX; i1++)
    {
        gridLB(i1,objectOriginY-1,ZERO_MINUS ) = fBottom[i1-objectOriginX+1][ZERO_PLUS  ] ;
        gridLB(i1,objectOriginY-1,MINUS_MINUS) = fBottom[i1-objectOriginX+1][PLUS_PLUS  ] ;
        gridLB(i1,objectOriginY-1,PLUS_MINUS ) = fBottom[i1-objectOriginX+1][MINUS_PLUS ] ;

        gridLB(i1,objectEndY+1   ,PLUS_PLUS  ) =    fTop[i1-objectOriginX+1][MINUS_MINUS] ;
        gridLB(i1,objectEndY+1   ,ZERO_PLUS  ) =    fTop[i1-objectOriginX+1][ZERO_MINUS ] ;
        gridLB(i1,objectEndY+1   ,MINUS_PLUS ) =    fTop[i1-objectOriginX+1][PLUS_MINUS ] ;
    }


    gridLB(objectOriginX-1,objectOriginY-1,MINUS_MINUS) =  fLeft[objectOriginY-1-objectOriginY+1][PLUS_PLUS  ] ;

    gridLB(objectEndX+1   ,objectEndY+1,PLUS_PLUS  ) = fRight[objectEndY+1-objectOriginY+1][MINUS_MINUS] ;

    gridLB(objectEndX+1,objectOriginY-1,PLUS_MINUS ) = fBottom[objectEndX+1-objectOriginX+1][MINUS_PLUS ] ;

    gridLB(objectOriginX-1,objectEndY+1   ,MINUS_PLUS ) =    fTop[objectOriginX-1-objectOriginX+1][PLUS_MINUS ] ;

}

void applyWallCylinder(ModelD2Q9<double> &myModel,Grid &gridLB,double fwall[][2*r+3][N_DV],int x,int y,int r)
{


//The square enclosing the circle.
    for (int i =0; i<=2*r+2; i++){
        for (int j=0; j<=2*r+2; j++){
            for (int dv=0; dv<N_DV; dv++){

                if(fwall[i][j][dv] != 0)
                    gridLB(x+i,y+j,dv) = fwall[i][j][dv];

            }}}



}

//collision
//f=f+2*beta*(fEQ-f)
void Collide(ModelD2Q9<double> &myModel,Grid &gridLB,double twoBeta, int step)
{
    double f_i[9];
    double x_i[9];
    double ximax(0.0);

    for(int i=1;i<=SIZE_X;i++)
        for(int j=1;j<=SIZE_Y;j++)
        {
            getMoments(myModel,gridLB,i,j);
            getfEQ(myModel);
            for(int k=0;k<N_DV;k++)
                gridLB(i, j, k) = gridLB(i, j, k) + twoBeta*(myModel.fEQ[k]- gridLB(i, j, k));
        }

}

void copyToGrid(ModelD2Q9<double> &myModel,Grid &gridLB,int x_coord,int y_coord)
{
    for(int k=0;k<N_DV;k++)
        gridLB(x_coord, y_coord, k)= myModel.fEQ[k];
}

void InitialConditions(ModelD2Q9<double> &myModel,Grid &gridLB)
{
    double x;
    for(int i=0;i<=SIZE_X+1;i++)
        for(int j=0;j<=SIZE_Y+1;j++)
        {
            x = ((double)i)/((double)SIZE_X);
            myModel.rho = 1.0;
            myModel.u1  = 0.0;//*sin(2*M_PI*x);
            myModel.u2  = 0.0;
            myModel.theta= myModel.theta0;
            getfEQ(myModel);
            copyToGrid(myModel,gridLB,i,j);
            //std::cout<<i<<"\t"<<j<<std::endl;
        }
}

void massConservationCheck(Grid &gridLB, int time)
{
    double mass=0.0;
    for(int i=1;i<=SIZE_X;i++)
        for(int j=1;j<=SIZE_Y;j++)
            for(int k=0;k<N_DV;k++)
                mass+=gridLB(i, j, k);
    mass/=SIZE_X*SIZE_Y;
    std::cout<<"at time step="<<time << " global mass=" <<mass<<std::endl;
}

void setModelParameters(ModelD2Q9<double> &myModel)
{
    myModel.ci_x[ZERO_ZERO]=0.0;
    myModel.ci_x[PLUS_ZERO]=1.0;
    myModel.ci_x[ZERO_PLUS]=0.0;
    myModel.ci_x[MINUS_ZERO]=-1.0;
    myModel.ci_x[ZERO_MINUS]=0.0;
    myModel.ci_x[PLUS_PLUS]=1.0;
    myModel.ci_x[MINUS_PLUS]=-1.0;
    myModel.ci_x[MINUS_MINUS]=-1.0;
    myModel.ci_x[PLUS_MINUS]=1.0;

    myModel.ci_y[ZERO_ZERO]=0.0;
    myModel.ci_y[PLUS_ZERO]=0.0;
    myModel.ci_y[ZERO_PLUS]=1.0;
    myModel.ci_y[MINUS_ZERO]=0.0;
    myModel.ci_y[ZERO_MINUS]=-1.0;
    myModel.ci_y[PLUS_PLUS]=1.0;
    myModel.ci_y[MINUS_PLUS]=1.0;
    myModel.ci_y[MINUS_MINUS]=-1.0;
    myModel.ci_y[PLUS_MINUS]=-1.0;

    myModel.wt[ZERO_ZERO]=16.0/36.0;
    myModel.wt[PLUS_ZERO]=4.0/36.0;
    myModel.wt[ZERO_PLUS]=4.0/36.0;
    myModel.wt[MINUS_ZERO]=4.0/36.0;
    myModel.wt[ZERO_MINUS]=4.0/36.0;
    myModel.wt[PLUS_PLUS]=1.0/36.0;
    myModel.wt[MINUS_PLUS]=1.0/36.0;
    myModel.wt[MINUS_MINUS]=1.0/36.0;
    myModel.wt[PLUS_MINUS]=1.0/36.0;

    myModel.dVopposite[ZERO_ZERO  ] = ZERO_ZERO   ;
    myModel.dVopposite[PLUS_ZERO  ] = MINUS_ZERO  ;
    myModel.dVopposite[ZERO_PLUS  ] = ZERO_MINUS  ;
    myModel.dVopposite[MINUS_ZERO ] = PLUS_ZERO   ;
    myModel.dVopposite[ZERO_MINUS ] = ZERO_PLUS   ;
    myModel.dVopposite[PLUS_PLUS  ] = MINUS_MINUS ;
    myModel.dVopposite[MINUS_PLUS ] = PLUS_MINUS  ;
    myModel.dVopposite[MINUS_MINUS] = PLUS_PLUS   ;
    myModel.dVopposite[PLUS_MINUS ] = MINUS_PLUS  ;

    myModel.theta0=1.0/3.0;
    myModel.oneBytheta0=3.0;
    myModel.sqrtTheta0=sqrt(myModel.theta0);

    for(int dv=0;dv<N_DV;dv++){
        myModel.cc[dv]  = myModel.ci_x[dv]*myModel.ci_x[dv] + myModel.ci_y[dv]*myModel.ci_y[dv];
        myModel.cx2[dv] = myModel.ci_x[dv]*myModel.ci_x[dv] ;
        myModel.cy2[dv] = myModel.ci_y[dv]*myModel.ci_y[dv];
    }

}

void printVtk(ModelD2Q9<double> &myModel,Grid &gridLB,int step, double beta, double dt)
{
    double alpha(2.0), ximax(0.0), ximin(0.0), alphaMax(0.0),temp(0.0),tempPsi(0.0);
    double f_i[9];double x_i[9];
    std::ofstream file;
    char fileName[250];
    sprintf(fileName,"./results/x-velocity_%d.vtk",step);
    file.open(fileName);

    //   vtk file header
    file<<"# vtk DataFile Version 3.0"<<std::endl<<"Velocity"<<std::endl<<"ASCII"<<std::endl<<"DATASET STRUCTURED_GRID"<<std::endl;
    file<<"DIMENSIONS "<<(SIZE_X)<<" "<<(SIZE_Y)<<" "<<(SIZE_Z)<<std::endl;
    file<<"POINTS "<<(SIZE_Y)*(SIZE_X)*(SIZE_Z)<<" double"<<std::endl;
    for (int i2=1 ; i2<=SIZE_Y;  i2++)
        for (int i1=1 ; i1<=SIZE_X;  i1++)
            file<<i1<<" "<<i2<<" "<< "1"<<std::endl;

    file<<"POINT_DATA "<<(SIZE_X)*(SIZE_Y)*(SIZE_Z)<<std::endl;
    file<<"VECTORS"<<" "<<"velocity"<<" "<<"double"<<std::endl;
    for (int i3=1 ; i3<=SIZE_Z; i3++)
        for (int i2=1 ; i2<=SIZE_Y;  i2++)
            for (int i1=1 ; i1<=SIZE_X;  i1++)
            {
                if(i1==1){
                    getMoments(myModel,gridLB, i1, i2);
                    file << myModel.u1 << "     " << myModel.u2 << "     " << myModel.rho<<std::endl;//(myModel.rho*myModel.u1*myModel.u1)+(myModel.rho/3.0)<< std::endl;
                    tempPsi = 0.0;
                }
                else{
                    getMoments(myModel,gridLB, i1-1, i2);
                    temp = myModel.u2;
                    getMoments(myModel,gridLB, i1, i2);
                    file << myModel.u1 << "     " << myModel.u2 << "     " <<myModel.rho<<std::endl;//(myModel.rho*myModel.u1*myModel.u1)+(myModel.rho/3.0)<< std::endl;
                    tempPsi = tempPsi-0.5*dt*(temp+myModel.u2);
                }
            }

    file.close();
}

void printVel(ModelD2Q9<double> &myModel,Grid &gridLB,int step, double beta, double dt, int objectOriginX, int objectOriginY)
{
//   std::cout<<"Printing vertical line at x = 2"<<std::endl;
//    int x = 2;//ENTER X-COORDINATE OF THE VERTICAL LINE
//    for (int j = SIZE_Y;j>=1;j--){
//       getMoments(myModel, gridLB,x,j);
//       std::cout<<myModel.u1<<std::endl;
//    }
//    std::cout<<"Printing horizontal line at y = 50"<<std::endl;
//    int y = 50;//ENTER Y-COORDINATE OF THE HORIZONTAL LINE
//    for (int i = SIZE_X;i>=1;i--){
//        getMoments(myModel, gridLB,i,y);
//        std::cout<<myModel.u1<<std::endl;
//    }

    int i1 = 2*objectOriginX + 2*r;
    int i2 = objectOriginY + r;
    i1 = 100;
    i2 = 100;
  getMoments(myModel,gridLB, i1, i2);
  std::cout << step*dt << " "<< myModel.u1 << " " << myModel.u2 << " " << sqrt(myModel.u1*myModel.u1 + myModel.u2*myModel.u2)<<std::endl ;

    i1 = 3*objectOriginX + 2*r;
    i2 = objectOriginY + r;

   // getMoments(myModel,gridLB, i1, i2);
   // std::cout << " "<< myModel.u1 << " " << myModel.u2 << " " << sqrt(myModel.u1*myModel.u1 + myModel.u2*myModel.u2) ;

    i1 = 2*objectOriginX + 2*r;
    i2 = objectOriginY + 3*r;

   // getMoments(myModel,gridLB, i1, i2);
   // std::cout << " "<< myModel.u1 << " " << myModel.u2 << " " << sqrt(myModel.u1*myModel.u1 + myModel.u2*myModel.u2) << std::endl;
   // getMoments(myModel,gridLB, 3, 100);
   // std::cout<<myModel.rho<<std::endl;
    //i1 = 1;
    //i2 = 1;
    //int dv;
    //double sum=0;
    //getMoments(myModel,gridLB, i1, i2);
    //for (dv=0;dv<N_DV;dv++) {
    //std::cout << gridLB(i1, i2, ZERO_PLUS) << " "<<gridLB(i1,i2,ZERO_MINUS)<<std::endl;
    // std::cout << gridLB(i1, i2, MINUS_MINUS) << " "<<gridLB(i1,i2,MINUS_PLUS)<<std::endl;
    // std::cout << gridLB(i1, i2, PLUS_MINUS) << " "<<gridLB(i1,i2,PLUS_PLUS)<<std::endl;

    //   sum  = sum + gridLB(i1, i2, dv);
    // }
    // std::cout<<sum<<std::endl;
}

inline void getStress(ModelD2Q9<double> &myModel,Grid &gridLB,int x_coord,int y_coord, double  &Pxx,double  &Pyy,double  &Pxy){

    Pxx = 0.0;
    Pyy = 0.0;
    Pxy = 0.0;

    for(int k=0;k<N_DV;k++)
    {
        Pxx += gridLB(x_coord, y_coord, k) * myModel.ci_x[k] * myModel.ci_x[k];
        Pyy += gridLB(x_coord, y_coord, k) * myModel.ci_y[k] * myModel.ci_y[k];
        Pxy += gridLB(x_coord, y_coord, k) * myModel.ci_x[k] * myModel.ci_y[k];
    }


}

void applyWallTopBottom(ModelD2Q9<double> &myModel,Grid &gridLB)
{
double Pxx, Pxy, Pyy;
    for(int i=2;i<=SIZE_X;i++)
    {
        // BOTTOM WALL
        gridLB(i, 1, ZERO_PLUS ) = gridLB(i,1,ZERO_MINUS);
        gridLB(i, 1, MINUS_PLUS) = gridLB(i,1,MINUS_MINUS);
        gridLB(i, 1, PLUS_PLUS ) = gridLB(i,1,PLUS_MINUS);
       // getMoments(myModel,gridLB,i,1);
      ///  std::cout<<myModel.u1<<" .."<<std::endl;
//	      gridLB(i, 1, ZERO_PLUS ) = gridLB(i,0,ZERO_MINUS);
//	      gridLB(i, 1, MINUS_PLUS) = gridLB(i,0,MINUS_MINUS);
//	      gridLB(i, 1, PLUS_PLUS ) = gridLB(i,0,PLUS_MINUS);

        //TOP WALL

        gridLB(i, SIZE_Y, ZERO_MINUS)  = gridLB(i,SIZE_Y,ZERO_PLUS);
        gridLB(i, SIZE_Y, PLUS_MINUS)  = gridLB(i,SIZE_Y,PLUS_PLUS);
        gridLB(i, SIZE_Y, MINUS_MINUS) = gridLB(i,SIZE_Y,MINUS_PLUS);
        //getStress(myModel,gridLB,i,SIZE_Y,Pxx,Pyy,Pxy);
        //std::cout<<Pxy<<".."<<std::endl;
//
//	      gridLB(i, SIZE_Y, ZERO_MINUS)  = gridLB(i,SIZE_Y+1,ZERO_PLUS);
//	      gridLB(i, SIZE_Y, PLUS_MINUS)  = gridLB(i,SIZE_Y+1,PLUS_PLUS);
//	      gridLB(i, SIZE_Y, MINUS_MINUS) = gridLB(i,SIZE_Y+1,MINUS_PLUS);
    }

//       for(int i=2;i<=SIZE_X-1;i++)
// 	   {
// 	      BOTTOM WALL
// 	      gridLB(i, 1, ZERO_PLUS ) = gridLB(i,0,ZERO_MINUS);
// 	      gridLB(i, 1, MINUS_PLUS) = gridLB(i,0,PLUS_MINUS);
// 	      gridLB(i, 1, PLUS_PLUS ) = gridLB(i,0,MINUS_MINUS);
//
// 	      TOP WALL
// 	      gridLB(i, SIZE_Y, ZERO_MINUS)  = gridLB(i,SIZE_Y+1,ZERO_PLUS);
// 	      gridLB(i, SIZE_Y, PLUS_MINUS)  = gridLB(i,SIZE_Y+1,MINUS_PLUS);
// 	      gridLB(i, SIZE_Y, MINUS_MINUS) = gridLB(i,SIZE_Y+1,PLUS_PLUS);
// 	   }
//
//
}

inline void getGradDist(ModelD2Q9<double> &myModel, double rho, double  ux, double  uy, double P, double pXXminuspYY, double pXY){

    double tmp1 = 4.5;

    double dot(0.0), T(1.0/3.0), Tinv(3.0);

    for(int dv=0;dv<N_DV;dv++)
    {
        dot = (ux*myModel.ci_x[dv] + uy*myModel.ci_y[dv])*Tinv*rho;

        myModel.fEQ[dv] = myModel.wt[dv]*(  rho + dot + 0.5*Tinv*Tinv*(  (P-rho*T)*(myModel.cc[dv]-2.0*T) + 0.5*pXXminuspYY*(myModel.cx2[dv] - myModel.cy2[dv]) + 2.0*pXY*myModel.ci_x[dv]*myModel.ci_y[dv] )   );
    }

}


void applyWallInletFeq(ModelD2Q9<double> &myModel,Grid &lbNode,double dt, double tau) {

    double rhoTemp;
    double denom;

    double Pxx, Pyy, Pxy, rho, ux, uy, k1, k2, k3, jy, sigmaxx, sigmayy;
    double factor = 2.0 * tau / (2.0 * tau + dt);
    double dtByTwoTau = 0.5 * dt / tau;

    for (int i2 = 2; i2 <= SIZE_Y - 1; i2++) {
        getMoments(myModel, lbNode, 4, i2);
        myModel.u1 = myModel.inlet*myModel.rho;
        myModel.u2 = 0.0;
        myModel.rho = 1.0;

//        std::cout<<myModel.rho<<"  "<<myModel.u1<<std::endl;

        getfEQ(myModel);
        rho = myModel.rho;
        ux = myModel.inlet;
        uy = myModel.u2;
        //for (int dv = 0; dv < N_DV; dv++) {
        //lbNode(1, i2,PLUS_ZERO )= myModel.fEQ[PLUS_ZERO ];
        //lbNode(1, i2,PLUS_PLUS )= myModel.fEQ[PLUS_PLUS ];
        //lbNode(1, i2,PLUS_MINUS)= myModel.fEQ[PLUS_MINUS];

        //}
//
        getMoments(myModel, lbNode, 4, i2);
        getStress(myModel, lbNode, 4, i2, Pxx, Pyy, Pxy);
        sigmaxx = Pxx - myModel.rho * myModel.u1 * myModel.u1 - myModel.rho * myModel.theta0;
        sigmayy = Pyy - myModel.rho * myModel.u2 * myModel.u2 - myModel.rho * myModel.theta0;

        Pxx = 1.0/3.0 + 1.0 * sigmaxx;
        Pyy = 1.0/3.0 + 1.0 * sigmayy;

        //Pxx -=ux*ux;
        //Pyy = factor*( Pyy + dtByTwoTau*(rho*uy*uy + rho/3.0) );
        // Pxy = factor*( Pxy + dtByTwoTau*(rho*ux*uy) );
        //getGradDist(myModel, rho, ux, uy, 0.5*(Pxx+Pyy) , Pxx-Pyy, Pxy);

        for (int dv = 0; dv < N_DV; dv++) {
            k1 = myModel.wt[dv] * (Pxx - myModel.rho * myModel.theta0) / (2.0 * myModel.theta0 * myModel.theta0) *
                 (myModel.ci_x[dv] * myModel.ci_x[dv] - myModel.theta0);
            k2 = myModel.wt[dv] * (Pyy - myModel.rho * myModel.theta0) / (2.0 * myModel.theta0 * myModel.theta0) *
                 (myModel.ci_y[dv] * myModel.ci_y[dv] - myModel.theta0);
            lbNode(1, i2, dv) = myModel.fEQ[dv]+ k1 + k2;
            lbNode(2, i2, dv) = myModel.fEQ[dv] + k1 + k2;
            lbNode(3, i2, dv) = myModel.fEQ[dv] + k1 + k2;
            // lbNode(1,i2,dv) = myModel.fEQ[dv];
        }
        /*  lbNode(1,i2,PLUS_PLUS ) = myModel.fEQ[PLUS_PLUS ]+ myModel.wt[PLUS_PLUS ]*(Pxx - myModel.rho*myModel.theta0)/(2.0*myModel.theta0*myModel.theta0)*(myModel.ci_x[PLUS_PLUS ]*myModel.ci_x[PLUS_PLUS ]-myModel.theta0) ;
          lbNode(1,i2,PLUS_MINUS) = myModel.fEQ[PLUS_MINUS]+ myModel.wt[PLUS_MINUS]*(Pxx - myModel.rho*myModel.theta0)/(2.0*myModel.theta0*myModel.theta0)*(myModel.ci_x[PLUS_MINUS]*myModel.ci_x[PLUS_MINUS]-myModel.theta0) ;
          lbNode(1,i2,PLUS_ZERO ) = myModel.fEQ[PLUS_ZERO ]+ myModel.wt[PLUS_ZERO ]*(Pxx - myModel.rho*myModel.theta0)/(2.0*myModel.theta0*myModel.theta0)*(myModel.ci_x[PLUS_ZERO ]*myModel.ci_x[PLUS_ZERO ]-myModel.theta0) ;
  */
        //getMoments(myModel,lbNode,2,i2);
        //std::cout<<"U1 is "<<myModel.u1<<" U2 is "<<myModel.u2<<" rho is "<<myModel.rho<<std::endl;
        //getMoments(myModel, lbNode, 1 ,i2);
        //lbNode(1,i2,0) = lbNode(1,i2,0)-(myModel.rho - 1.0);
    }
    getMoments(myModel, lbNode, 4, 1);
    myModel.u1 = myModel.inlet*myModel.rho;
    myModel.u2 = 0.0;
    myModel.rho = 1.0;
    getfEQ(myModel);
    getMoments(myModel, lbNode, 4, 1);
    getStress(myModel, lbNode, 4, 1, Pxx, Pyy, Pxy);
    sigmaxx = Pxx - myModel.rho * myModel.u1 * myModel.u1 - myModel.rho * myModel.theta0;
    sigmayy = Pyy - myModel.rho * myModel.u2 * myModel.u2 - myModel.rho * myModel.theta0;
    Pxx = 1.0 / 3.0 + 1.0 * sigmaxx;
    Pyy = 1.0 / 3.0 + 1.0 * sigmayy;
    for (int dv = 0; dv < N_DV; dv++) {
        k1 = myModel.wt[dv] * (Pxx - myModel.rho * myModel.theta0) / (2.0 * myModel.theta0 * myModel.theta0) *
             (myModel.ci_x[dv] * myModel.ci_x[dv] - myModel.theta0);
        k2 = myModel.wt[dv] * (Pyy - myModel.rho * myModel.theta0) / (2.0 * myModel.theta0 * myModel.theta0) *
             (myModel.ci_y[dv] * myModel.ci_y[dv] - myModel.theta0);
        lbNode(1, 1, dv) = myModel.fEQ[dv]+ k1 + k2;
        lbNode(2, 1, dv) = myModel.fEQ[dv]+k1+k2;
        lbNode(3, 1, dv) = myModel.fEQ[dv]+k1+k2;

    }

    getMoments(myModel, lbNode, 1 ,1);
    lbNode(1,1,0) = lbNode(1,1,0)-(myModel.rho - 1.0);
    //getStress(myModel, lbNode, 1, 1, Pxx, Pyy, Pxy);
    ///std::cout<<Pxy<<std::endl;

   //getMoments(myModel,lbNode,2,1);
   //std::cout<<"U1 is "<<myModel.u1<<" U2 is "<<myModel.u2<<" rho is "<<myModel.rho<<std::endl;

/*
    lbNode(1,1,PLUS_ZERO) = myModel.fEQ[PLUS_ZERO];
    lbNode(1,1,PLUS_PLUS) = myModel.fEQ[PLUS_PLUS];
    lbNode(1,1,PLUS_MINUS) = myModel.fEQ[PLUS_MINUS];
    lbNode(1,1,ZERO_PLUS) = lbNode(1,1,ZERO_MINUS);

   //lbNode(1,1,ZERO_MINUS) = lbNode(1,1,ZERO_PLUS);
    lbNode(1,1,MINUS_PLUS) = -lbNode(1,1,PLUS_MINUS)+lbNode(1,1,PLUS_PLUS)+lbNode(1,1,MINUS_MINUS);
    lbNode(1,1,ZERO_ZERO) = myModel.fEQ[ZERO_ZERO];
    lbNode(1,1,ZERO_MINUS) = myModel.fEQ[ZERO_MINUS];
    lbNode(1,1,MINUS_ZERO) = myModel.fEQ[MINUS_ZERO];
    lbNode(1,1,MINUS_MINUS) = myModel.fEQ[MINUS_MINUS];


//            for (int dv=0;dv<N_DV;dv++){
//               lbNode(1,1,dv) = myModel.fEQ[dv];
//            }*/
    getMoments(myModel, lbNode, 4, SIZE_Y);
    myModel.u1 = myModel.inlet*myModel.rho;
    myModel.u2 = 0.0;
    myModel.rho = 1.0;
    getfEQ(myModel);
    getMoments(myModel, lbNode, 4 ,SIZE_Y);
    getStress(myModel, lbNode, 4, SIZE_Y, Pxx, Pyy, Pxy);
    sigmaxx = Pxx - myModel.rho * myModel.u1 * myModel.u1 - myModel.rho * myModel.theta0;
    sigmayy = Pyy - myModel.rho * myModel.u2 * myModel.u2 - myModel.rho * myModel.theta0;
    Pxx = 1.0 / 3.0 + 1.0 * sigmaxx;
    Pyy = 1.0/3.0 + 1.0 * sigmayy;

    for (int dv = 0; dv < N_DV; dv++) {
        k1 = myModel.wt[dv] * (Pxx - myModel.rho * myModel.theta0) / (2.0 * myModel.theta0 * myModel.theta0) *
             (myModel.ci_x[dv] * myModel.ci_x[dv] - myModel.theta0);
        k2 = myModel.wt[dv] * (Pyy - myModel.rho * myModel.theta0) / (2.0 * myModel.theta0 * myModel.theta0) *
             (myModel.ci_y[dv] * myModel.ci_y[dv] - myModel.theta0);
        lbNode(1, SIZE_Y, dv) = myModel.fEQ[dv]+ k1 + k2;
        lbNode(2, SIZE_Y, dv) = myModel.fEQ[dv]+ k1 + k2;
        lbNode(3, SIZE_Y, dv) = myModel.fEQ[dv]+ k1 + k2;
    }
    getMoments(myModel, lbNode, 1 ,SIZE_Y);
    lbNode(1,SIZE_Y,0) = lbNode(1,SIZE_Y,0)-(myModel.rho - 1.0);
    //getStress(myModel, lbNode, 1, SIZE_Y, Pxx, Pyy, Pxy);
   // std::cout<<Pxy<<".."<<std::endl;
    //getMoments(myModel, lbNode, 1 ,SIZE_Y);
   // std::cout<<myModel.rho<<"..."<<std::endl;
   //  getMoments(myModel,lbNode,2,SIZE_Y);
   //for (int dv=0;dv<N_DV;dv++){
   //    std::cout<<lbNode(1,SIZE_Y,dv)<<" ";
   //}
   //std::cout<<std::endl;
    // std::cout<<"U1 is "<<myModel.u1<<" U2 is "<<myModel.u2<<" rho is "<<myModel.rho<<std::endl;

/*
    //    lbNode(1,SIZE_Y,PLUS_ZERO) = myModel.fEQ[PLUS_ZERO];
//    lbNode(1,SIZE_Y,PLUS_PLUS) = myModel.fEQ[PLUS_PLUS];
//    lbNode(1,SIZE_Y,PLUS_MINUS) = myModel.fEQ[PLUS_MINUS];
//    lbNode(1,SIZE_Y,ZERO_MINUS) = lbNode(1,SIZE_Y,ZERO_PLUS);
//    lbNode(1,SIZE_Y,MINUS_MINUS) = -lbNode(1,SIZE_Y,PLUS_PLUS)+lbNode(1,SIZE_Y,PLUS_MINUS)+lbNode(1,SIZE_Y,MINUS_MINUS);
////            for (int dv=0;dv<N_DV;dv++){
////               lbNode(1,SIZE_Y,dv) = myModel.fEQ[dv];
////            }
//    lbNode(1,SIZE_Y,ZERO_ZERO) = myModel.fEQ[ZERO_ZERO];
//    lbNode(1,SIZE_Y,ZERO_PLUS) = myModel.fEQ[ZERO_PLUS];
//    lbNode(1,SIZE_Y,MINUS_PLUS) = myModel.fEQ[MINUS_PLUS];
//    lbNode(1,SIZE_Y,MINUS_ZERO) = myModel.fEQ[MINUS_ZERO];
    //getMoments(myModel,lbNode,1,SIZE_Y);
    //std::cout<<"U1 is "<<myModel.u1<<" U2 is "<<myModel.u2<<" rho is "<<myModel.rho<<std::endl;
}
*/
}
inline void  applyWallOutletFgrads( ModelD2Q9<double> &myModel, Grid &lbNode, double Pout, double dt, double tau ){

    double fTemp[N_DV];
    double rho, ux, uy, T, Pxx, Pyy, Pxy;

    double dtByTwoTau = 0.5*dt/tau;
    double factor     = 2.0*tau/(2.0*tau+dt);

    int i1 = SIZE_X+1;

    for (int i2=1 ;i2<= SIZE_Y;  i2++)
    {

        getMoments(myModel,lbNode, i1, i2);

        rho = myModel.rho;//Pout/T;
        ux = myModel.u1;
        uy = myModel.u2;

        getStress(myModel, lbNode, i1, i2, Pxx, Pyy, Pxy)   ;

        Pxx = factor*( Pxx + dtByTwoTau*(rho*ux*ux + rho/3.0) );
        Pyy = factor*( Pyy + dtByTwoTau*(rho*uy*uy + rho/3.0) );
        Pxy = factor*( Pxy + dtByTwoTau*(rho*ux*uy          ) );

        getGradDist( myModel, rho, ux, uy, 0.5*(Pxx+Pyy) , Pxx-Pyy, Pxy)	;
//         getGradDist( myModel, rho, ux, uy, 0.5*rho*(ux*ux+uy*uy) + Pout , factor*Pxx-Pyy, factor*Pxy)	;
        // 	      getGradDist( fEq, rho, ux, uy, Pout, Pxx-Pyy, Pxy)	;


//         getfEQ(myModel);

        for(int dv=0;dv<N_DV;dv++)
        {
            lbNode(SIZE_X,i2,dv)= myModel.fEQ[dv];
        }

    }
}


inline void  applyWallOutletFgradsNew( ModelD2Q9<double> &myModel, Grid &lbNode, double Pout, double dt, double tau ){

    double fTemp[N_DV];
    double rho, ux, uy, T, Pxx, Pyy, Pxy, jy,a,b,c,k,k1,k2,k3,k4;
    double sumW, sumWc, rhoPlus, jxPlus, pxxplus, jxNew, rhoNew, fZERO;
    sumW  = myModel.wt[ZERO_MINUS] + myModel.wt[ZERO_ZERO] + myModel.wt[ZERO_PLUS] + myModel.wt[PLUS_MINUS] + myModel.wt[PLUS_ZERO] + myModel.wt[PLUS_PLUS] ;//k1
    sumWc = myModel.wt[PLUS_MINUS] + myModel.wt[PLUS_ZERO] + myModel.wt[PLUS_PLUS] ;//k2
//     a11 = (2*sumW)-((myModel.wt[ZERO_MINUS]+myModel.wt[ZERO_PLUS]+myModel.wt[PLUS_ZERO]+2*(myModel.wt[PLUS_PLUS]+myModel.wt[PLUS_MINUS]))/(2*myModel.theta0));
//     a12 = sumWc/myModel.theta0;
//     a13 = (sumWc-(mymodel.theta0*sumW))/(2*myModel.theta0*myModel.theta0);
//     a21 = (2*sumWc)-(myModel.wt[PLUS_ZERO]+2*(myModel.wt[PLUS_PLUS]+myModel.wt[PLUS_MINUS])/(2*myModel.theta0));
//     a22 = a12;
//     a23 = (sumWc-(mymodel.theta0*sumWc))/(2*myModel.theta0*myModel.theta0);
    double dtByTwoTau = 0.5*dt/tau;
    double factor     = 2.0*tau/(2.0*tau+dt);

    for (int i2=2 ;i2<= SIZE_Y-1;  i2++)
    {
        //getMoments(myModel,lbNode, SIZE_X-1, i2);
        //ux = myModel.u1;

        getMoments(myModel,lbNode, SIZE_X+1, i2);
        rhoNew = lbNode(SIZE_X,i2,ZERO_ZERO)/myModel.wt[ZERO_ZERO]/( 1.0 - 0.5*(myModel.u1*myModel.u1 + myModel.u2*myModel.u2)/myModel.theta0 );
        jy = 0.5*( lbNode(SIZE_X,i2,ZERO_PLUS) - lbNode(SIZE_X,i2,ZERO_MINUS) )*myModel.theta0 / myModel.wt[ZERO_MINUS];
//         myModel.u2 = jy/rhoNew;
//         myModel.rho = rhoNew;
        rho = rhoNew;
        ux = myModel.u1;
        uy = jy/rho;
        getStress(myModel, lbNode, SIZE_X+1, i2, Pxx, Pyy, Pxy);
        Pxx = factor*( Pxx + dtByTwoTau*(rho*ux*ux + rho/3.0) );
        Pyy = factor*( Pyy + dtByTwoTau*(rho*uy*uy + rho/3.0) );
        fZERO = lbNode(SIZE_X,i2,ZERO_ZERO);
        Pxy = (myModel.theta0*myModel.theta0/2*myModel.wt[PLUS_PLUS])*((lbNode(SIZE_X,i2,PLUS_PLUS)-lbNode(SIZE_X,i2,PLUS_MINUS))-(lbNode(SIZE_X,i2,ZERO_PLUS)-lbNode(SIZE_X,i2,ZERO_MINUS)));
//         k2 = lbNode(SIZE_X,i2,PLUS_ZERO)/myModel.wt[PLUS_ZERO];
//         k = (Pxx*(1-myModel.theta0)-(Pyy*myModel.theta0))/(2*myModel.theta0*myModel.theta0);
//         k3 = 2-(1/(2*myModel.theta0));
//         k4 = k2-k;
//         a = k4/(2*myModel.theta0);
//         b = k1/myModel.theta0;
//         c = (k1*k3)-k4;
//         ux  = (-b + sqrt(b*b - (4*a*c)))/(2*a);
//         if( b*b < (4*a*c)){
//             std::cout<<"Nope"<<std::endl;
//             abort();
//         }

//         rho = k4/(ux+k3);
//         uy = jy/rho;
//         Pxx = rho*myModel.theta0 + rho*ux*ux;
//         Pyy = rho*myModel.theta0 + rho*uy*uy;

//         ux  = myModel.u1;
//         uy  = myModel.u2;
//         Pxy = rho*ux*uy          ;//factor*( Pxy + dtByTwoTau*(rho*ux*uy          ) );

//         rhoPlus = lbNode(SIZE_X,i2,ZERO_MINUS) + lbNode(SIZE_X,i2,ZERO_ZERO) + lbNode(SIZE_X,i2,ZERO_PLUS) + lbNode(SIZE_X,i2,PLUS_MINUS) + lbNode(SIZE_X,i2,PLUS_ZERO) + lbNode(SIZE_X,i2,PLUS_PLUS);
//         jxPlus = lbNode(SIZE_X,i2,PLUS_MINUS) + lbNode(SIZE_X,i2,PLUS_ZERO) + lbNode(SIZE_X,i2,PLUS_PLUS);

//         rhoNew = ( 15.0*fZERO + 18.0*jxPlus + 4.0*Pyy - 18.0*rhoPlus ) /6.0; ;
//         jxNew  = ( -357.0*fZERO - 234.0*jxPlus - 88.0*Pyy + 558.0*rhoPlus) /162.0 ;
//         Pxx    = ( 21.0*fZERO + 36.0*jxPlus + 5.0*Pyy - 36.0*rhoPlus ) /3.0;
//         std::cout << "moments are" << jxNew << "  " << rhoNew << "  ";

//         rho = fZERO/myModel.wt[ZERO_ZERO]/(1.0 - (ux*ux+uy*uy)/(2.0*myModel.theta0) );

        getGradDist(myModel, rho, ux, uy, 0.5*(Pxx+Pyy) , Pxx-Pyy, Pxy)	;
//         getGradDist( myModel, rho, ux, uy, 0.5*rho*(ux*ux+uy*uy) + Pout , factor*Pxx-Pyy, factor*Pxy)	;
        /*
           lbNode(SIZE_X,i2,MINUS_MINUS)= myModel.fEQ[MINUS_MINUS];
           lbNode(SIZE_X,i2,MINUS_ZERO )= myModel.fEQ[MINUS_ZERO ];
           lbNode(SIZE_X,i2,MINUS_PLUS )= myModel.fEQ[MINUS_PLUS ];*/
        //New idea
//      k1 = myModel.rho - (lbNode(SIZE_X,i2,PLUS_MINUS)+lbNode(SIZE_X,i2,PLUS_PLUS)+lbNode(SIZE_X,i2,PLUS_ZERO)+lbNode(SIZE_X,i2,ZERO_MINUS)+lbNode(SIZE_X,i2,ZERO_PLUS)+lbNode(SIZE_X,i2,ZERO_ZERO));
//      k2 = myModel.rho*myModel.u2 - (lbNode(SIZE_X,i2,ZERO_ZERO)+lbNode(SIZE_X,i2,ZERO_ZERO)-lbNode(SIZE_X,i2,ZERO_ZERO)-lbNode(SIZE_X,i2,ZERO_ZERO));
//      k3 = Pyy - (lbNode(SIZE_X,i2,ZERO_ZERO)+lbNode(SIZE_X,i2,ZERO_ZERO)+lbNode(SIZE_X,i2,ZERO_ZERO)+lbNode(SIZE_X,i2,ZERO_ZERO));
//      lbNode(SIZE_X,i2,MINUS_ZERO) = k1-k3;
//      lbNode(SIZE_X,i2,MINUS_PLUS) = (k2+k3)/2;
//      lbNode(SIZE_X,i2,MINUS_MINUS) = (k3-k2)/2;
//         myModel.u1 = ux;
//         myModel.u2 = uy;
//         myModel.rho = rho;
        //getfEQ(myModel);
//         std::cout<<"ux is "<<ux<<std::endl;
//         lbNode(SIZE_X,i2,MINUS_MINUS)= myModel.fEQ[MINUS_MINUS];
//         lbNode(SIZE_X,i2,MINUS_ZERO )= myModel.fEQ[MINUS_ZERO ];
//         lbNode(SIZE_X,i2,MINUS_PLUS )= myModel.fEQ[MINUS_PLUS ];
        for(int dv=0;dv<N_DV;dv++)
        {
            lbNode(SIZE_X,i2,dv)= myModel.fEQ[dv];
        }
    }

    for(int dv=0;dv<N_DV;dv++)
    {
        lbNode(SIZE_X,1,dv)      = lbNode(SIZE_X-1,1,dv);
        lbNode(SIZE_X,SIZE_Y,dv) = lbNode(SIZE_X-1,SIZE_Y,dv);
    }
}

int main()
{
    ModelD2Q9 <double> lbModel;
    Grid lbGrid;
    setModelParameters(lbModel);
    InitialConditions(lbModel,lbGrid);
    massConservationCheck(lbGrid, 0);
    //simulation parameters
    lbModel.Ma  = 0.05;
//   double n;
//   n = 1.6;
    lbModel.inlet  = lbModel.Ma*lbModel.sqrtTheta0;
    lbModel.inletActual =lbModel.Ma*lbModel.sqrtTheta0;
    double Re      = 200.0;
    double Kn      = lbModel.Ma/Re;
    double refLen  = 2.0*r;
    double kinVisc = lbModel.inletActual*refLen/Re;
    double tau     = kinVisc*lbModel.oneBytheta0;
    double dt      = refLen/(2.0*r);
    double tauNdim = tau/dt;
    double beta    = 1.0/(2.0*tauNdim+1.0);
    std::cout<<"Tau is "<<tau<<std::endl;
    int simulationTime = 50*((int)(SIZE_X/lbModel.inlet));

    const int objectOriginX = 200;
    const int objectOriginY = 120;
    const int objectLengthX = 10;
    const int objectLengthY = 10;

    double fwall[2*r+3][2*r+3][N_DV] = {0.0};

    //printVtk(lbModel,lbGrid,0,beta,dt);

    for(int time=1;time<=simulationTime;time++)
    {
        //prepareWallOutlet(lbModel,lbGrid,SIZE_X);

        Collide(lbModel,lbGrid,2.0*beta,time);
        // entropicCollide(lbModel,lbGrid,beta);

   // if(time%10==0)
       // printVel(lbModel,lbGrid,time,beta,dt,objectOriginX, objectOriginY);

        if(time%100==0)
        {
            printVtk(lbModel,lbGrid,time,beta,dt);
            massConservationCheck(lbGrid, time);
        }
//     periodicAll(lbModel,lbGrid);
        prepareWallOutlet(lbModel,lbGrid,SIZE_X);
        prepareWallTopBottom(lbModel,lbGrid);
        prepareWallInlet(lbModel,lbGrid);
        //prepareWallCylinder(lbModel, lbGrid, fwall, objectOriginX, objectOriginY, r);
        //applyWallInletFeq(lbModel,lbGrid);

        advection(lbGrid);
        //prepareWallInlet(lbModel,lbGrid);
        //applyWallTopBottom(lbModel,lbGrid);
        applyWallInletFeq(lbModel,lbGrid,dt, tau );
//        applyWallTopBottom(lbModel,lbGrid);
        //prepareWallOutlet(lbModel,lbGrid,SIZE_X);
     applyWallTopBottom(lbModel,lbGrid);
//     applyWallOutletFeq(lbModel,lbGrid);
        //applyWallOutletFgrads( lbModel, lbGrid, 1.0/3.0, dt, tau );
        applyWallOutletFgradsNew( lbModel, lbGrid, 1.0/3.0, dt, tau );
        // applyWallCylinder(lbModel, lbGrid, fwall, objectOriginX, objectOriginY, r);
       // if(time%100 == 0) {//ENTER TIME INTERVALS AFTER WHICH IT SHOULD PRINT THE VELOCITIES
            printVel(lbModel, lbGrid, time, beta, dt, objectOriginX, objectOriginY);
       // // printvert(lbModel,lbGrid,time,beta,dt,objectOriginX, objectOriginY);
        //prepareWallOutlet(lbModel,lbGrid,SIZE_X);
//        Collide(lbModel,lbGrid,2.0*beta,time);

    }
    return 0;
}
