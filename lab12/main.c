#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double fun(double x){
    return log(pow(x,3)+3*pow(x,2)+x+0.1)*sin(18*x);
}

double xi(double a,int i,double h){
    return a+i*h;
}

void trapez( int n,double (*D)[n+1],double a,double b){

    for(int w=0;w<=n;w++){
        double h=(b-a)/pow(2,w);
        int N = pow(2,w);
        double sum=0.0;
        for(int i=0;i<N;i++){
            sum+=(h/2)*(fun(xi(a,i,h))+fun(xi(a,i+1,h)));
        }
        D[w][0]=sum;
    }   
}

void trzy( int n,double (*D)[n+1],double a,double b){
    for(int w=0;w<=n;w++){
        double h=(b-a)/(3*pow(2,w));
        int N = 3*pow(2,w);
        double sum=0.0;
        for(int i=0;i<N/3;i++){
            sum+=(3*h/8)*(fun(xi(a,3*i,h))+3*fun(xi(a,3*i+1,h))+3*fun(xi(a,3*i+2,h))+fun(xi(a,3*i+3,h)));
        }
        D[w][0]=sum;
    }   
}


void solve(void (*fun)(int n,double(*)[n+1],double, double)){
    
    //n
    int n=8;

    //przedział [a,b]
    double a=0.;
    double b=1.;

    //tablica
    double D[n+1][n+1];
    fun(n,D,a,b);    

    for(int k=1;k<=n;k++){
        for(int w=k;w<=n;w++){
            D[w][k]=(pow(4,k)*D[w][k-1]-D[w-1][k-1])/(pow(4,k)-1);
        }
    }

    for(int i=0;i<=n;i++){
        printf("%d %g %g\n",i,D[i][0],D[i][i]);
    }
}


int main(){
    printf("===Metoda trapezów===\nw   D[w][0]   D[w][w]\n");
    solve(trapez);

    printf("===Metoda 3/8===\nw   D[w][0]   D[w][w]\n");
    solve(trzy);

    return 0;
}