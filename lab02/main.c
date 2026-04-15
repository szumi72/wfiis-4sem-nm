#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "nrutil.h"
#include "nrutil.c" 
#include "ludcmp.c"
#include "lubksb.c"

#define N 4

void printMatrix(float **A);
float findMatrixMax(float **A);

int main(){
    float **A, **A_og;
    A=matrix(1,N,1,N);
    A_og=matrix(1,N,1,N);


    int *indx = ivector(1,N);
    float d;
    int delta = 0;

    float **L,**U;
    L=matrix(1,N,1,N);
    U=matrix(1,N,1,N);


    for(int i=1;i<=N;++i){
        for(int j=1;j<=N;++j){
            A[i][j]=1.0/(i+j+delta);
            A_og[i][j]=1.0/(i+j+delta);
        }
    }
    printf("Matrix A:\n");
    printMatrix(A);

    ludcmp(A,N,indx,&d);

    printf("Matrix LU:\n");
    printMatrix(A);   

    for(int i=1;i<=N;i++){
        for(int j=i;j<=N;j++){
            U[i][j]=A[i][j];
            L[i][j]=0; 
        }
        for(int j=1;j<i;j++){
            U[i][j]=0;
             
        }
    }

    for(int i=1;i<=N;i++){
        for(int j=1;j<=i;j++){
            L[i][j]=A[i][j];  
                 
        }
        L[i][i]=1;
    }


    printf("Matrix U:\n");
    printMatrix(U);

    printf("Matrix L:\n");
    printMatrix(L);

    float detA=1;
    for(int i=1;i<=N;++i){
        printf(" %g ",A[i][i]);
        detA *=A[i][i];
    }
    detA *=d;
    printf("\ndet(A): %g \n",detA);

    float *b1;
    b1 = vector(1,N);
    
    float **A_inv;
    A_inv = matrix(1,N,1,N);

    for(int i=1;i<=N;++i){
        for(int i=1;i<=N;++i){
            b1[i]=0;        
        }
        b1[i]=1;
        lubksb(A,N,indx,b1);        
        A_inv[1][i] = b1[1];
        A_inv[2][i] = b1[2];
        A_inv[3][i] = b1[3];
        A_inv[4][i] = b1[4];
    }
    printf("Matrix A(-1):\n");
    printMatrix(A_inv);

    float **Id = matrix(1,N,1,N);
    
    for(int i=1;i<=N;i++){
        for(int j=1;j<=N;j++){
            Id[i][j]=0;
            for(int k=1;k<=N;k++){
                Id[i][j]+=A_og[i][k]*A_inv[k][j];
            }
        }
    }
    printf("Matrix A*A(-1)\n");
    printMatrix(Id);

    float m1;
    m1= findMatrixMax(A_og);
    float m2 = findMatrixMax(A_inv);

    
    printf("cond A = %g, cond A(-1) = %g\n",m1,m2);
    printf("cond A*A(-1) = %g\n",m1*m2);


    free_matrix(A,1,N,1,N);
    free_matrix(A_inv,1,N,1,N);
    free_matrix(A_og,1,N,1,N);  
    free_matrix(Id,1,N,1,N);
    free_matrix(L,1,N,1,N);
    free_matrix(U,1,N,1,N);
    free_vector(b1,1,N);    
    free_ivector(indx,1,N);
}

void printMatrix(float **A){
    for(int i=1;i<=N;++i){
        for(int j=1;j<=N;++j){
            printf(" %8g ",A[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}


float findMatrixMax(float **A){
    float max=fabs(A[1][1]);
    for(int i=1;i<=N;i++){
        for(int j=1;j<=N;j++){
            if(fabs(A[i][j])>max){
                max = fabs(A[i][j]);
            }
        }
    }

    return max;
}