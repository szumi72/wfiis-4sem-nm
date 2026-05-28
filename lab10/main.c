#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double g(double x, double y){
    return (x*x-4*x + y*y - 4*y + x*y);
}

void A_grad(double(*A_T)[2],double* grad,double* vect){
    vect[0] = A_T[0][0]*grad[0]+A_T[0][1]*grad[1];
    vect[1]= A_T[1][0]*grad[0]+A_T[1][1]*grad[1];
    
}

double newt(double* r, double omega,double(*A_T)[2]){
    
    double grad[2];

    double dx=0.01;
    double dy=0.01;

    double r_n[2];
    double vect[2];
    int it;
    printf("(%g %g) ",r[0],r[1]);
    for(int i=0;i<1000;i++){
        it = i+1;
        
        
        grad[0]=(g(r[0]+dx,r[1]) - g(r[0]-dx,r[1]))/(2*dx);
        grad[1]=(g(r[0],r[1]+dy)-g(r[0],r[1]-dy))/(2*dy);
        //printf("grad %g %g\n",grad[0],grad[1]);
        A_grad(A_T,grad,vect);
        
        r_n[0]=r[0]-omega*vect[0];
        r_n[1]=r[1]-omega*vect[1];

        if(sqrt(pow((r_n[0]-r[0]),2)+pow((r_n[1]-r[1]),2))<1e-6){
            break;

        }
        r[0]=r_n[0];
        r[1]=r_n[1];


    }
    
    printf("%d %g %g\n",it,r_n[0],r_n[0]);


}

int main(){

    double A[2][2];
    A[0][0]=2;
    A[0][1]=1;
    A[1][0]=1;
    A[1][1]=2;

    double b[2]={-4,-4};

    double A_T[2][2];

    double det = (1.0)/(A[0][0]*A[1][1]-A[0][1]*A[1][0]);
    A_T[0][0]=A[1][1]*det;
    A_T[0][1]=(-1)*A[0][1]*det;
    A_T[1][0]=(-1)*A[1][0]*det;
    A_T[1][1]=A[0][0]*det;

    // for(int i=0;i<2;i++){
    //     for(int j=0;j<2;j++){
    //         printf("%10g",A_T[i][j]);
    //     }
    //     printf("\n");
    // }

    double r_min[2];

    r_min[0]=(-1)*(A_T[0][0]*b[0]+A_T[0][1]*b[0]);
    r_min[1]=(-1)*(A_T[1][0]*b[1]+A_T[1][1]*b[1]);

    printf("r_min = %g %g\n\n",r_min[0], r_min[1]);
   

    double r_0[2];

    printf("Omega = 1\n");
    //
    r_0[0]=0;
    r_0[1]=0;    
    newt(r_0,1,A_T);
    
    //
    r_0[0]=10;
    r_0[1]=-10;    
    newt(r_0,1,A_T);

    //
    r_0[0]=100;
    r_0[1]=100;    
    newt(r_0,1,A_T);

    r_0[0]=500;
    r_0[1]=500;    
    newt(r_0,1,A_T);

    printf("Różne Omegi\n");
    r_0[0]=10;
    r_0[1]=10;    
    newt(r_0,0.1,A_T);

    r_0[0]=10;
    r_0[1]=10;
    newt(r_0,0.4,A_T);

    r_0[0]=10;
    r_0[1]=10;
    newt(r_0,0.7,A_T);


}