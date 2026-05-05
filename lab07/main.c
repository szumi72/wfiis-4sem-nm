#include<stdlib.h>
#include<stdio.h>
#include<math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

double W_n(double x, int n,double* xm, double* ym){
    double result=0;
    
    for(int j=0;j<=n;j++){
        double il=1;
        for(int k=0;k<=n;k++){
            if(k!=j){
                il *= (x-xm[k])/(xm[j]-xm[k]); 
            }
        }
        result += ym[j]*il;
    }

    return result;
}

double f_x(double x){
    return exp(-1*(x*x));
}

void solve(int n, FILE* fp){
    double xm[n+1];
    double ym[n+1];

    double x_min=-5.0;
    double x_max=5.0;
    double h=((x_max-x_min)/(double)n);
    for(int i = 0;i<n+1;i++){        
        xm[i]=x_min+(i*h);        
        ym[i]=f_x(xm[i]);        
    }

       
    for(double x=x_min;x<x_max;x+=0.01){
        fprintf(fp,"%g %g\n",x,W_n(x,n,xm,ym));
    }

    
}

void solve2(int n, FILE* fp){
    double xm[n+1];
    double ym[n+1];

    double x_min=-5.0;
    double x_max=5.0;
    
    for(int i = 0;i<n+1;i++){        
        xm[i]=(0.5)*((x_max-x_min)*cos(M_PI*((2*i+1)/(2.0*n+2)))+(x_min+x_max));        
        ym[i]=f_x(xm[i]);        
    }

    
    for(double x=x_min;x<x_max;x+=0.01){
        fprintf(fp,"%g %g\n",x,W_n(x,n,xm,ym));
    }
}



int main(){

    FILE *fp = fopen("zad_1.dat","w");
    solve(5,fp);
    fprintf(fp,"\n\n");
    solve(10,fp);
    fprintf(fp,"\n\n");
    solve(15,fp);
    fprintf(fp,"\n\n");
    solve(20,fp);
    fprintf(fp,"\n\n");
    fclose(fp);


    FILE* fp2 = fopen("zad_2.dat","w");
    solve2(5,fp2);
    fprintf(fp2,"\n\n");
    solve2(10,fp2);
    fprintf(fp2,"\n\n");
    solve2(15,fp2);
    fprintf(fp2,"\n\n");
    solve2(20,fp2);
    fprintf(fp2,"\n\n");
    fclose(fp2);
}