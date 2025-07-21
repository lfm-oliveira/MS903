#include <stdio.h>
#include <stdlib.h>

double funcao(int n,double* x);
void deriva(int n, double* x);

int main(){
    int n;
    double *x;
    scanf("%d",&n);
    x= malloc(n*sizeof(double));
    if(x == NULL)
        exit(1);
    for(int k = 0; k<n;k++)
        scanf("%lf", &x[k]);
    deriva(n,x);
    return 0;
}

double funcao(int n,double* x){
    double sum = 0;
    for(int i = 0; i<n;i++){
        sum+=x[i]*x[i];
    }
    return sum;
}

void deriva(int n, double* x){
    double *grad, **hess, h, xminus[n], xplus[n];
    int i,j;
    h = 1e-3;
    grad = malloc(n*sizeof(double));
    hess = malloc(n*sizeof(double));
    if(grad == NULL || hess == NULL)
        exit(1);
    printf("Gradiente de f em x = [");

    for(i=0;i<n;i++){
        hess[i] = malloc(sizeof(double));
        if(hess[i] == NULL)
            exit(1);
    }
    for(i = 0; i<n;i++){
        xplus[i] = x[i];
        xminus[i] = x[i];
        printf(" %lf",x[i]);
    }
    printf("]: \n");

    for(j = 0; j<n;j++){
        xplus[j] +=h;
        xminus[j]-=h;
        grad[j] = (funcao(n,&xplus[0]) - funcao(n,&xminus[0]))/(2*h);
        hess[j][j] = (funcao(n,&xplus[0]) - 2*funcao(n,x)+f(n,&xminus))/(h*h);
        
        for(i = 0; i<j-1;i++){
            hess[i][j] = (funcao(n,&xplus[0]) - 2*funcao(n,x)+f(n,&xminus))/(4*h*h);
        }
        
        
        printf("%lf\n",grad[j]);



        xplus[j] = xminus[j] = x[j];
    }


    free(grad);
    for(i = 0; i<n;i++)
        free(hess[i]);
    free(hess);
}

