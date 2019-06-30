#include <iostream>
#include <fstream>
#include <stdio.h>
#include <iomanip>
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
#define SIZE_X 600
#define SIZE_Y 125
#define SIZE_Z 1

//define structure for the D2Q9 model
template <typename dataType>
struct ModelD2Q9
{
  dataType ci_x[N_DV];
  dataType ci_y[N_DV];
  dataType cc[N_DV];
  dataType wt[N_DV];
  dataType theta0,oneBytheta0,sqrtTheta0;
  dataType c2[N_DV];
  dataType rho;
  dataType u1,u2,theta;
  dataType inlet;
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

void prepareWallInletOutlet(ModelD2Q9<double> &myModel,Grid &gridLB)
{	


    for(int dv=0;dv<N_DV;dv++)
	  {  
	      for(int j=1;j <= SIZE_Y;j++)
	      { //copy from right wall to right ghost
		gridLB(SIZE_X+1, j, dv)   = gridLB(SIZE_X-1, j, dv);
		//copy from left wall to left ghost
		gridLB(0, j, dv)    = gridLB(1, j,  dv);
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

void applyWallTopBottom(ModelD2Q9<double> &myModel,Grid &gridLB)
{
      for(int i=2;i<=SIZE_X-1;i++)
	   {
           //BOTTOM WALL
           //CLASSICAL BOUNDARY CONDITIONS
           getMoments(myModel, gridLB, i,0);
           getfEQ(myModel);
           double sumBottom, sumfeq,sumTop;
	       sumBottom = gridLB(i,0,ZERO_MINUS)+gridLB(i,0,PLUS_MINUS)+gridLB(i,0,MINUS_MINUS);
	       sumfeq = myModel.fEQ[ZERO_MINUS]+myModel.fEQ[PLUS_MINUS]+myModel.fEQ[MINUS_MINUS];
           gridLB(i, 1, ZERO_PLUS ) = sumBottom/sumfeq * myModel.fEQ[ZERO_PLUS];
           gridLB(i, 1, MINUS_PLUS) = sumBottom/sumfeq * myModel.fEQ[MINUS_PLUS];
           gridLB(i, 1, PLUS_PLUS ) = sumBottom/sumfeq * myModel.fEQ[PLUS_PLUS];
           //specular boundary condition
           //	      gridLB(i, 1, ZERO_PLUS ) = gridLB(i,0,ZERO_MINUS);
           //	      gridLB(i, 1, MINUS_PLUS) = gridLB(i,0,MINUS_MINUS);
           //	      gridLB(i, 1, PLUS_PLUS ) = gridLB(i,0,PLUS_MINUS);
           //TOP WALL
           //CLASSICAL BOUNDARY CONDITIONS
           getMoments(myModel, gridLB, i,0);
           getfEQ(myModel);
           sumTop = gridLB(i,0,ZERO_PLUS)+gridLB(i,0,MINUS_PLUS)+gridLB(i,0,PLUS_PLUS);
           sumfeq = myModel.fEQ[ZERO_PLUS]+myModel.fEQ[MINUS_PLUS]+myModel.fEQ[PLUS_PLUS];
           gridLB(i, 1, ZERO_MINUS ) =  sumTop/sumfeq * myModel.fEQ[ZERO_MINUS];
           gridLB(i, 1, PLUS_MINUS) =   sumTop/sumfeq * myModel.fEQ[PLUS_MINUS];
           gridLB(i, 1, MINUS_MINUS ) = sumTop/sumfeq * myModel.fEQ[MINUS_MINUS];
           //SPECULAR CONDITIONS
           //	      gridLB(i, SIZE_Y, ZERO_MINUS)  = gridLB(i,SIZE_Y+1,ZERO_PLUS);
           //	      gridLB(i, SIZE_Y, PLUS_MINUS)  = gridLB(i,SIZE_Y+1,PLUS_PLUS);
           //	      gridLB(i, SIZE_Y, MINUS_MINUS) = gridLB(i,SIZE_Y+1,MINUS_PLUS);
	   }
}

void applyWallInletOutlet(ModelD2Q9<double> &myModel,Grid &gridLB)
{	

      double rhoTemp;	
      double factor,denom;

      for(int j=1;j <= SIZE_Y;j++)
	   {
         //OUTLET-------RIGHT WALL
          getMoments(myModel,gridLB,SIZE_X+1,j);
          getfEQ(myModel);
          double x1,x2,x3,x0,rho0,rhoplus,dem;
          x1 = myModel.fEQ[MINUS_PLUS]+(gridLB(SIZE_X,j,MINUS_PLUS)-myModel.fEQ[PLUS_MINUS])+0.5*(gridLB(SIZE_X,j,ZERO_MINUS)- gridLB(SIZE_X,j,ZERO_PLUS)-myModel.fEQ[ZERO_MINUS]+myModel.fEQ[ZERO_PLUS]);
          x2 = myModel.fEQ[MINUS_ZERO]+(gridLB(SIZE_X,j,PLUS_ZERO) - myModel.fEQ[PLUS_ZERO]);
          x3 = myModel.fEQ[PLUS_MINUS]+(gridLB(SIZE_X,j,PLUS_PLUS)-myModel.fEQ[PLUS_PLUS])-0.5*(gridLB(SIZE_X,j,ZERO_MINUS)- gridLB(SIZE_X,j,ZERO_PLUS)-myModel.fEQ[ZERO_MINUS]+myModel.fEQ[ZERO_PLUS]);
          rho0 = gridLB(SIZE_X,j,ZERO_PLUS)+gridLB(SIZE_X,j,ZERO_MINUS)+gridLB(SIZE_X,j,ZERO_ZERO);
          rhoplus = gridLB(SIZE_X,j,PLUS_PLUS)+gridLB(SIZE_X,j,PLUS_ZERO)+gridLB(SIZE_X,j,PLUS_MINUS);
          dem = 1/(1-(myModel.u1));
          x0 = gridLB(SIZE_X,j,ZERO_ZERO)+myModel.rho - dem*(rho0+(2*rhoplus));
          gridLB(SIZE_X,j,MINUS_ZERO) = x2;
          gridLB(SIZE_X,j,MINUS_MINUS)= x3;
          gridLB(SIZE_X,j,MINUS_PLUS) = x1;
          gridLB(SIZE_X,j,ZERO_ZERO) = x0;

//          getMoments(myModel,gridLB,SIZE_X+1,j);
//          getfEQ(myModel);
//          gridLB(SIZE_X,j,MINUS_ZERO)=myModel.fEQ[MINUS_ZERO];
//          gridLB(SIZE_X,j,MINUS_MINUS)=myModel.fEQ[MINUS_MINUS];
//          gridLB(SIZE_X,j,MINUS_PLUS)=myModel.fEQ[MINUS_PLUS];
//          for (int dv=0;dv<N_DV;dv++){
//             gridLB(SIZE_X,j,dv)=myModel.fEQ[dv];
//          }
///////////////////////////////////////////////////////////////////////////////////
//INLET---------------------------LEFT WALL
           myModel.rho=1.0;
           myModel.u1=myModel.inlet;
           myModel.u2=0.0;

           getfEQ(myModel);
           for (int dv=0;dv<N_DV;dv++){
              gridLB(1,j,dv)=myModel.fEQ[dv];
          }
//           myModel.denominator=1.0/(myModel.fEQ[ZERO_MINUS]+myModel.fEQ[PLUS_MINUS]+myModel.fEQ[MINUS_MINUS]);
//            factor = gridLB(0,j,ZERO_MINUS)+ gridLB(0,j,PLUS_MINUS)+ gridLB(0,j,MINUS_MINUS);
//            gridLB(1,j,PLUS_ZERO)=factor*myModel.denominator*myModel.fEQ[PLUS_ZERO];
//            gridLB(1,j,PLUS_MINUS)=factor*myModel.denominator*myModel.fEQ[PLUS_MINUS];
//            gridLB(1,j,PLUS_PLUS)=factor*myModel.denominator*myModel.fEQ[PLUS_PLUS];
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
//Adding new functions

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
    for(int i=1;i<=SIZE_X;i++)
      for(int j=1;j<=SIZE_Y;j++)
      {
	    x = ((double)i)/((double)SIZE_X);
	    myModel.rho = 1.0;
	    myModel.u1  = 0.01*sin(2*M_PI*x);
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
	
	myModel.theta0=1.0/3.0;
	myModel.oneBytheta0=3.0;
	myModel.sqrtTheta0=sqrt(myModel.theta0);

	for(int dv=0;dv<N_DV;dv++)
	  myModel.cc[dv] = myModel.ci_x[dv]*myModel.ci_x[dv] + myModel.ci_y[dv]*myModel.ci_y[dv];

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
        file << myModel.u1 << "     " << myModel.u2 << "     " << "0.0" << std::endl;
	tempPsi = 0.0;
     }
     else{
        getMoments(myModel,gridLB, i1-1, i2);
        temp = myModel.u2;
        getMoments(myModel,gridLB, i1, i2);
        file << myModel.u1 << "     " << myModel.u2 << "     " << "0.0" << std::endl;
        tempPsi = tempPsi-0.5*dt*(temp+myModel.u2);	
     }
   }
   file.close();
}


void printstrouhal(ModelD2Q9<double> &myModel,Grid &gridLB,int step, double beta, double dt) {
    for (int i3 = 1; i3 <= SIZE_Z; i3++) {
        for (int i2 = 1; i2 <= SIZE_Y; i2++) {
            for (int i1 = 1; i1 <= SIZE_X; i1++) {
                if (i2 == (0.5 * SIZE_Y) && i1 == (0.5 * SIZE_X)) {
                    std::ofstream outdata;
                    outdata.open("velocity.txt", std::ios_base::app);
                    outdata <<sqrt(myModel.u1*myModel.u1 + myModel.u2*myModel.u2)<< std::endl;
                    outdata.close();
                }

            }

        }
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
  lbModel.Ma         = 0.1;
  lbModel.inlet = lbModel.Ma*lbModel.sqrtTheta0;

  const int objectOriginX = 100;
  const int objectOriginY = 60;
  const int objectLengthX =5;
  const int objectLengthY =5;
  double Re      = 250.0;
  double Kn      = lbModel.Ma/Re;
  double refLen  = objectLengthX;
  double kinVisc = lbModel.inlet*objectLengthX/Re;
  double tau     = kinVisc*lbModel.oneBytheta0;
  double dt      = 1.0;
  double tauNdim = tau/dt;
  double beta    = 1.0/(2.0*tauNdim+1.0);
  
  int simulationTime = 50*((int)(SIZE_X/lbModel.inlet));

  double fWallLeft[objectLengthY+2][9];
  double fWallRight[objectLengthY+2][9];
  double fWallTop[objectLengthX+2][9];
  double fWallBottom[objectLengthX+2][9];
  
  printVtk(lbModel,lbGrid,0,beta,dt);

  for(int time=1;time<=simulationTime;time++)
  {
    Collide(lbModel,lbGrid,2.0*beta,time);
    // entropicCollide(lbModel,lbGrid,beta);
    printstrouhal(lbModel, lbGrid, time, beta, dt);
    if(time%10000==0)
    {
        printVtk(lbModel,lbGrid,time,beta,dt);
	    massConservationCheck(lbGrid, time);
    }
//     periodicAll(lbModel,lbGrid);

    prepareWallTopBottom(lbModel,lbGrid);
    prepareWallObject(lbModel,lbGrid,fWallLeft,fWallRight,fWallTop,fWallBottom,objectOriginX,objectOriginY,objectLengthX,objectLengthX);
    prepareWallInletOutlet(lbModel,lbGrid);
    advection(lbGrid);
    applyWallTopBottom(lbModel,lbGrid);
    applyWallInletOutlet(lbModel,lbGrid);
    applyWallObject(lbModel,lbGrid,fWallLeft,fWallRight,fWallTop,fWallBottom,objectOriginX,objectOriginY,objectLengthX,objectLengthX);

}
  return 0;
}
