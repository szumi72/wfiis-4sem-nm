#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "tqli.c"
#include "pythag.c" 
#include "tred2.c"
#include "nrutil.h"
#include "nrutil.c" 


#define N 5

void printMatrix(float **A);

int main(){

    float** A,**A_og;
    A=matrix(1,N,1,N);
    A_og=matrix(1,N,1,N);
    for(int i=1;i<=N;++i){
        for(int j=1;j<=N;++j){
            A[i][j]=sqrt(i+j);
            A_og[i][j]=sqrt(i+j);
        }
    
    }
    printf("Mtrix A:\n");
    printMatrix(A);

    float *d, *e;
    d=vector(1,N);
    e=vector(1,N);

    tred2(A,N,d,e);
    printf("Mtrix P:\n");
    printMatrix(A);

    
    
    float**T;
    T=matrix(1,N,1,N);
    for(int i=1;i<=N;++i){
        for(int j=1;j<=N;++j){
            T[i][j]=0;            
        }    
    }

    for(int i=1;i<=N;i++){
        for(int j=1;j<=N;j++){            
            if(i==j){
                T[i][j]=d[i];
                if(i<N && j<N){                    
                    T[i+1][j]=e[i+1];
                    T[i][j+1]=e[i+1];                    
                }             
            }
        }
    
    }
    printf("Mtrix T:\n");
    printMatrix(T);

    tqli(d,e,N,A);
    printf("wekrory wlasne macierzy A:\n");
    printMatrix(A);

    printf("wartości własne A:\n");
    for(int k=1;k<=N;k++){
        printf("%g\n",d[k]);
    }
    printf("\n");


    float**C;
    C=matrix(1,N,1,N);

    for(int i=1;i<=N;i++){
        for(int j=1;j<=N;j++){
            C[i][j]=0;
            for(int k=1;k<=N;k++){
                C[i][j]+=(A_og[i][k]*A[k][j]);
            }
        }
    }

    printf("Mtrix C:\n");
    printMatrix(C);
    

    float* beta;
    beta=vector(1,N); 
    float* beta2;
    beta2=vector(1,N); 
    for(int i=1;i<=N;i++){
        beta[i]=0;
        beta2[i]=0;
    }
    for(int i=1;i<=N;i++){        
        for(int k=1;k<=N;k++){
            beta[k]+=C[i][k]*A[i][k];
        }
    }   
    for(int i=1;i<=N;i++){        
        for(int k=1;k<=N;k++){
            beta2[k]+=(A[i][k]*A[i][k]);
        }
    }   
    for(int i=1;i<=N;i++){
        beta[i]/=beta2[i];
    }

    printf("beta:\n");
    for(int k=1;k<=N;k++){
        printf("%g\n",beta[k]);
    }
    

    free_vector(d,1,N);
    free_vector(e,1,N);
    free_matrix(A,1,N,1,N);
    free_matrix(T,1,N,1,N);
    free_matrix(A_og,1,N,1,N);
    free_matrix(C,1,N,1,N);
    free_vector(beta,1,N);
    free_vector(beta2,1,N);
}

void printMatrix(float **A){
    for(int i=1;i<=N;++i){
        for(int j=1;j<=N;++j){
            printf(" %10g ",A[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}
