#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include <gsl/gsl_linalg.h>




double f_1(double x){
    return (1.0/(1+x*x));
}

double f_2(double x){
    return cos(2*x);
}

void print_matrix(gsl_matrix* A, int n){
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            printf("%5g",gsl_matrix_get(A,i,j));
        }
        printf("\n");
    }
}

void print_vector(gsl_vector* v, int n){
    for(int j=0;j<n;j++){
        printf("%g ",gsl_vector_get(v,j));
    }
    printf("\n");
}

void wyzM(double* xw, double* yw,double* m, int n,double ALPHA,double BETA){
    gsl_matrix *A = gsl_matrix_calloc(n,n);
    gsl_vector *d = gsl_vector_calloc(n);
    m[0] = ALPHA;
    m[n-1] = BETA;

    double x_min=-5;
    double x_max=5;
    double h=(x_max-x_min)/((double)(n-1.0));   

    for(int i=1;i<n-1;i++){       
        double di = (6.0/(h+h))*(((yw[i+1]-yw[i])/h)-((yw[i]-yw[i-1])/h));
        gsl_vector_set(d,i,di);
        
    }
    //print_vector(d,n);

    for(int i=1;i<n-1;i++){
        gsl_matrix_set(A,i,i,2);
        gsl_matrix_set(A,i,i-1,0.5);
        gsl_matrix_set(A,i-1,i,0.5);       
    }
    gsl_matrix_set(A,0,0,1);
    gsl_matrix_set(A,n-1,n-1,1);
    gsl_matrix_set(A,0,1,0);  
    gsl_matrix_set(A,n-2,n-1,0.5);  

    //print_matrix(A,n);

    gsl_linalg_HH_svx(A,d);

    print_vector(d,n);

    for(int i=0;i<n;i++){
        m[i] = gsl_vector_get(d,i);
    }


    gsl_matrix_free(A);
    gsl_vector_free(d);
    return;
}

double wyzSx(double* xw, double* yw,double* m, int n,double x){
    int i=1;
   double x_min=-5;
    double x_max=5;
    double h=(x_max-x_min)/((double)(n-1.0));

    for(i=1;i<n;i++){
        if(x >= xw[i-1] && x <= xw[i]){
            break;
        }
    }

    double Ai = (yw[i] - yw[i-1])/h - (h/6.0)*(m[i] - m[i-1]);
    double Bi = yw[i-1] - m[i-1]*h*h/6.0;


    double s=m[i-1]*pow((xw[i]-x),3)/(6.0*h)+m[i]*pow((x-xw[i-1]),3)/(6.0*h) + Ai*(x-xw[i-1]) + Bi;

    return s;
}

void pochodne_fun(int n, double (*fun)(double), FILE *plik){

    double xw[n], yw[n], m[n];

    double x_min = -5, x_max = 5;
    double h = (x_max - x_min)/(n-1);

    for(int i=0;i<n;i++){
        xw[i] = x_min + i*h;
        yw[i] = fun(xw[i]);
    }

    wyzM(xw,yw,m,n,0.0,0.0);

    for(int i=0;i<n;i++){
        double d2 = (fun(xw[i]-0.01) - 2*fun(xw[i]) + fun(xw[i]+0.01))/(0.01*0.01);

        fprintf(plik,"%g %g %g\n", xw[i], m[i], d2);
    }
}

void solve(int n, double (*fun)(double), FILE *plik){

    double ALPHA =0.0;
    double BETA = 0.0;

    double xw[n];
    double yw[n];    
    double x_min=-5;
    double x_max=5;
    double x_delta=(x_max-x_min)/((double)(n-1));
    for(int i=0;i<n;i++){
        xw[i]=x_min+x_delta*i;
    }
    for(int i=0;i<n;i++){
        yw[i]=fun(xw[i]);
    }
    
    printf("\n\n");
    
    double m[n];

    wyzM(xw,yw,m,n,ALPHA,BETA);

    for(double x = x_min;x<=x_max; x+=0.01){
        fprintf(plik,"%g %g\n",x,wyzSx(xw,yw,m,n,x));
    }
    fprintf(plik,"\n\n");

 
}

int main(){

    FILE *f1_file = fopen("f1.dat","w");
    FILE *f2_file = fopen("f2.dat","w");
    FILE *pochodne = fopen("pochodne.dat","w");

    solve(5, f_1, f1_file);
    solve(8, f_1, f1_file);
    solve(21, f_1, f1_file);

    solve(5, f_2, f2_file);
    solve(8, f_2, f2_file);
    solve(21, f_2, f2_file);

    pochodne_fun(10, f_1, pochodne);

    fclose(f1_file);
    fclose(f2_file);
    fclose(pochodne);

}