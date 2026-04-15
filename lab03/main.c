#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "nrutil.h"
#include "nrutil.c" 

#define N 2001
double sum(double* a);

int main(){
    double *b,*d_0,*d_1,*d_2,*x_s,*x_n;
    b = dvector(1,N);
    d_0 = dvector(1,N);
    d_1 = dvector(1,N);
    d_2 = dvector(1,N);
    x_n = dvector(1,N);
    x_s = dvector(1,N);

    double omega = 1;
    double h = 0.02;
    double F0 = 0.0;
    double beta = 0.4;
    double OMEGA = 0.8;

    double a1,a2,a3;
    a1 = 1;
    a2 = omega*omega*h*h - 2 - beta*h;
    a3 = 1+beta*h;

       
    //tu był błąd wcześniej - odwrotnie miałam a3 i a1
    for(int i=1;i<=N;i++){
        d_0[i]=a3;
        d_1[i]=a2;
        d_2[i]=a1;       
        x_s[i]=0.0;
        x_n[i]=0.0;
    }

    b[1]=1;
    b[2]=0;

    for(int i=3;i<=N;i++){      
        b[i]=F0*sin(OMEGA*h*(i-1))*h*h;       
    }

    d_0[1]=1;
    d_0[2]=1;
    d_1[1]=0;
    d_1[2]=-1;
    d_2[1]=0;
    d_2[2]=0;
    

    int it=0;
   
    
    while(it<100000){
        ++it;
       
        x_n[1]=((1.0/d_0[1])*(b[1]));
        x_n[2]=((1.0/d_0[2])*(b[2]-d_1[2]*x_s[1]));        
        for(int i=3;i<=N;i++){
           x_n[i]=((1.0/d_0[i])*(b[i]-d_1[i]*x_s[i-1]-d_2[i]*x_s[i-2]));
        }
        
        double sum_n=sum(x_n);
        double sum_s=sum(x_s);
        if(fabs(sum_n-sum_s)<1.e-6){
            break;
        }
        for(int i=1;i<=N;i++){
            x_s[i]=x_n[i];
        }
       
    }
    printf("\n\n");
    
    for(int i=1;i<=N;i++){
        printf("%g  %g\n",(i-1)*h,x_n[i]);
    }
    printf("Iteracje: %d\n",it);
    


    free_dvector(b,1,N);
    free_dvector(d_0,1,N);
    free_dvector(d_1,1,N);
    free_dvector(d_2,1,N);
    free_dvector(x_n,1,N);
    free_dvector(x_s,1,N);

}




double sum(double* a){
    double sum=0;
    for(int i=1;i<=N;i++){
        sum+=(a[i]*a[i]);
    }
    return sum;
}