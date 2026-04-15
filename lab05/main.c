#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 7
#define IT_MAX 12

void printMatrix(double (*A)[N]);
double skalar(double *v1,double *v2);
double norm(double *v1);
void vectxMatrix(double *v1, double (*A)[N], double* res);
void multiplyMatrix(double (*A)[N],double (*B)[N],double (*C)[N]);
void transposiotionMatrix(double (*A)[N]);

int main(){
    double XT[N][N]={0};
    double A[N][N]={0};   
    double xp[N]={0};
    double xk[N]={0};
    double lambda[N]={0};     

    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            A[i][j]=(1.0)/(sqrt(2+fabs(i-j)));
        }
    }
    printMatrix(A);

    for(int k=0;k<N;k++){
        // zerowanie xp
        for(int t=0;t<N;t++){
            xp[t]=1.0;
        }        
        for(int i=0;i<IT_MAX;i++){ 
            
            vectxMatrix(xp,A,xk);
            for(int j=0;j<k;j++){
                double scalar = skalar(xk,XT[j]);
                for(int w=0;w<N;w++){
                    xk[w]=xk[w]-(scalar*XT[j][w]);
                }
            }
            lambda[k]=(skalar(xk,xp)/skalar(xp,xp));
            
            printf("%d %g\n",i+1,lambda[k]);            
            double norm_xk = norm(xk);
            for(int w=0;w<N;w++){
                xp[w]=(xk[w]/norm_xk);
            }
        }
        for(int i=0;i<N;i++){
            XT[k][i]=xp[i];
        }
        printf("\n\n");
        
    }
    printf("Transpownowana macierz X\n");
    printMatrix(XT);

    double D[N][N]={0};
    multiplyMatrix(XT,A,D);

    //teraz macierz XT nie jest juz transponowana
    transposiotionMatrix(XT);
    printf("macierz X\n");
    printMatrix(XT);

    double D1[N][N]={0};
    multiplyMatrix(D,XT,D1);
    printf("Macierz D\n");
    printMatrix(D1);

}

void printMatrix(double (*A)[N]){
    for(int i=0;i<N;++i){
        for(int j=0;j<N;++j){
            printf(" %6g ",A[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

double skalar(double *v1,double *v2){
    double sum=0;
    for(int i=0;i<N;i++){
        sum+=v1[i]*v2[i];
    }
    return sum;
}

double norm(double *v1){
    return(sqrt(skalar(v1,v1)));
}

void vectxMatrix(double *v1, double (*A)[N], double* res){
    for(int i=0;i<N;i++){
        res[i]=0;
    }    
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            res[i]+=A[i][j]*v1[j];
        }
    }    
    
}

void multiplyMatrix(double (*A)[N],double (*B)[N],double (*C)[N]){
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            C[i][j]=0;
            for(int k=0;k<N;k++){
                C[i][j]+=(A[i][k]*B[k][j]);
            }
        }
    }
}

void transposiotionMatrix(double (*A)[N]){
    double temp;
    for(int i=0;i<N;i++){
        for(int j=i;j<N;j++){
            temp = A[i][j];
            A[i][j]=A[j][i];
            A[j][i]=temp;
        }
    }
}
