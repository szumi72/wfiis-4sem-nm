#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_fft_complex.h>
#define frand() ((double)rand())/(RAND_MAX+1.0)

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

double yi(int i, double omega){
    return sin(omega*i)+sin(2*omega*i)+sin(3*omega*i);
}

void solve(int k){
    int N=pow(2,k);
    double y[N*2];
    double omega = 2*(2*M_PI)/N;
    
    for(int i=0;i<N;i++){
        int j=2*i;    
        double szum = (2*frand())-1;
        y[j]=yi(i,omega)+szum;
        y[j+1]=0.;
    }

    FILE *f = fopen("y12.dat","w");
    for(int i=0;i<N;i++){
        fprintf(f,"%d %g\n",i,y[2*i]);
    }

    fprintf(f,"\n\n");

   


    gsl_fft_complex_radix2_forward(y,1,N);

    FILE *fp=fopen("fft8.dat","w");
    for(int i=0;i<N;i++){
        fprintf(fp,"%d %g %g\n",i,y[2*i],y[2*i+1]);
    }
    fclose(fp);

    double max=sqrt((pow(y[0],2)+pow(y[1],2)));

    for(int i=0;i<N;i++){
        double temp = sqrt((pow(y[i*2],2)+pow(y[i*2+1],2)));
        if(max<temp){
            max = temp;
        }        
    }

    printf("max = %g\n",max);

    for(int i=0;i<N;i++){

        double modul = sqrt((pow(y[i*2],2)+pow(y[i*2+1],2)));
        if(modul<(max/2.0)){
            y[i*2]=0;
            y[2*i+1]=0;
        }
    }

    

    gsl_fft_complex_radix2_backward(y,1,N);
    
    for(int i=0;i<N;i++){
        fprintf(f,"%d %g\n",i,y[2*i]/N);
    }
    fclose(f);

}


int main(){

    //solve(8);    
    //solve(10);
    solve(12);
    

    return 0;
}