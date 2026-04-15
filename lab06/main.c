#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double fx(double x){
    return (x-1.2)*(x-2.3)*(x-3.3)*(x-3.3);
}

double ux(double x, double dx){
    return fx(x)/((fx(x+dx)-fx(x-dx))/(2*dx));
}

void fun(double x0, double x1){

    for(int k=1;k<=10000;k++){        
        double x2 = x1 - ((x1-x0)/(fx(x1)-fx(x0)))*fx(x1);        
        printf("%d, %g, %g, %g\n",k,x2,fabs(x2-x1),fx(x2));
        if(fabs(x2-x1)<1e-6){
            break;
        }
        x0=x1;
        x1=x2;
    }
    printf("\n\n");

}

void fun_mod(double x0, double x1,double dx){
    for(int k=1;k<=10000;k++){        
        double x2 = x1 - ((x1-x0)/(ux(x1,dx)-ux(x0,dx)))*ux(x1,dx);        
        printf("%d, %g, %g, %g\n",k,x2,fabs(x2-x1),fx(x2));
        if(fabs(x2-x1)<1e-6){
            break;
        }
        x0=x1;
        x1=x2;
    }
    printf("\n\n");

}


int main(){
    fun(0.9,1.0);
    fun(1.7,1.75);
    fun(3.7,3.65);
    fun_mod(3.7,3.65,0.1);
    fun_mod(3.7,3.65,0.001);

}