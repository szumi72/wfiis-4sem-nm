#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define N 201

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define frand() ((double)rand())/(RAND_MAX+1.0)




double fun_clean(double x){

    double x_min = -4.0;
    double x_max = 4.0;
    double sigma = (x_max - x_min)/16.0;
    double x_0 = 2.0;
    return sin((14*M_PI*x)/(x_max-x_min))*(exp((-1)*((x-x_0)*(x-x_0))/(2*sigma*sigma))+exp((-1)*((x+x_0)*(x+x_0))/(2*sigma*sigma)));
}

double fun(double x){   
    return fun_clean(x) + ((frand()-0.5)/5);
}

double F_approx(double(*phi)[201],int k,double* yw,int m){
    double f=0;
    for(int j=0;j<=m;j++){
        double Cj=0;
        double Sj=0;
        for(int i=0;i<N;i++){
            Cj+=yw[i]*phi[j][i];
            Sj+=phi[j][i]*phi[j][i];
            
        }
        f+=(Cj/Sj)*phi[j][k];       

    }
    return f;
}

int main(){
    double xw[N];
    double yw[N];
    double x_min = -4.0;
    double x_max = 4.0;
    double h = (x_max-x_min)/(double)(N-1.0);
    double x_0 = 2.0;
    double sigma = (x_max - x_min)/16.0;

    for(int i=0;i<N;i++){
        xw[i]=(x_min+(h*i));
        //printf("%d %g\n",i,xw[i]);
    }
    for(int i=0;i<N;i++){
        yw[i]=fun(xw[i]);
    }


    double phi[51][N];

    for(int i=0;i<N;i++){
        phi[0][i] = 1;       
    }

    double A1=0;
    double mian=0;
    double licz=0;
    for(int i=0;i<N;i++){
        mian += 1;
        licz += xw[i];
    }
    A1=licz/mian;

    for(int k=0;k<N;k++){
        phi[1][k] = (xw[k]-A1)*phi[0][k];
    }


    double Aj1=0;
    double licz_Aj1=0;
    double mian_Aj1=0;   
    double Bj=0;
    double licz_Bj=0;
    double mian_Bj=0;

    for(int j=1;j<50;j++){
        licz_Aj1=0;
        mian_Aj1=0;
        licz_Bj=0;
        mian_Bj=0;
        for(int i=0;i<N;i++){
            licz_Aj1+=xw[i]*(phi[j][i]*phi[j][i]);
            mian_Aj1+=(phi[j][i]*phi[j][i]);
            licz_Bj+=xw[i]*phi[j-1][i]*phi[j][i];
            mian_Bj+=phi[j-1][i]*phi[j-1][i];
        }
        Aj1 = licz_Aj1/mian_Aj1;
        Bj = licz_Bj/mian_Bj;
        //printf("%d %g %g\n",j,Aj1,Bj);

        for(int k=0;k<N;k++){
            phi[j+1][k] = (xw[k]-Aj1)*phi[j][k]-(Bj*phi[j-1][k]);
        }
    }
    
    FILE *fg = fopen("Gram.dat","w");
   for(int k=0;k<N;k++){
        fprintf(fg,"%g ", xw[k]);  
        for(int j=0;j<7;j++){
            fprintf(fg,"%g ", phi[j][k]/phi[j][0]); 
        }

        fprintf(fg,"\n");
    }

    fclose(fg);

    FILE *fpkt = fopen("pkt.dat","w");

    for(int i=0;i<N;i++){
        fprintf(fpkt,"%g %g\n",xw[i],yw[i]);
    }
    fclose(fpkt);

    FILE *fp = fopen("approx.dat","w");

    fprintf(fp,"\n\n");
    
    for(int k=0;k<N;++k){
        fprintf(fp,"%g %g\n",xw[k],F_approx(phi,k,yw,10));
    }
    fprintf(fp,"\n\n");
    
    for(int k=0;k<N;++k){
        fprintf(fp,"%g %g\n",xw[k],F_approx(phi,k,yw,30));
    }
    fprintf(fp,"\n\n");
    
    for(int k=0;k<N;++k){
        fprintf(fp,"%g %g\n",xw[k],F_approx(phi,k,yw,50));
    }

    fclose(fp);


}