#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

clock_t start_time;

void tic() {
    start_time = clock();
}

void tac() {
    clock_t end_time = clock();
    double elapsed_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    printf("Elapsed time: %f seconds\n", elapsed_time);
}


#define M 2000
/*Dado x in R^n, epsilon > 0, M > 0
k = 0
Enquanto ||grad(f(x))|| >= epsilon e k < M, faça:
    Escolha d in R^n tal que grad(f(x))^t d < 0 (podemos tomar d = - grad(f(x))
    t =1
    Enquanto f(x+td) >= f(x), faça t <- t/2
    x<-x+td
    k<-k+1*/
double* quadratica(int n, double *x, int flag);

double* rosenbrook(int n, double *x, int flag);

double max(double *x, int n);

int main(){
    tic();
    int n,k,i,j,l;
    double *x,*d,*a,*b, aux, y_ant;
    double t, epsilon = 1e-4, alpha = 1e-4, maxgrad;
    double quad = 0;
    scanf("%d", &n);
    double y[n];
    x = malloc(n*sizeof(double));
    if(x == NULL)
        exit(1);
    for(i = 0; i<n;i++)
        scanf("%lf", &x[i]);
    k=0;
    maxgrad = max(rosenbrook(n,x,1),n);
    b = rosenbrook(n,x,0); /*b aponta para f(x)*/
    while(maxgrad >= epsilon && k < M){
        /*Calculo da direcao de descida e y para t = 1*/
        d = rosenbrook(n,x,1); /*d aponta para + grad f (x)*/
        t = 1;
        quad = 0;
        for(l = 0; l<n;l++){
            quad -=d[l]*d[l];
            y[l] = x[l]-d[l];
        }
        a = rosenbrook(n,y,0); /*a aponta para f(y)*/
        aux = *b + alpha*quad;
        if(*a <= aux){
            while(*a <= aux){
                t*=2;
                for(j = 0; j<n;j++)
                    y[j] = x[j]-t*d[j];
                y_ant = *a;
                free(a);
                /*free(b);*/
                a = rosenbrook(n,y,0);
                aux = *b + alpha*t*quad;
            }
            t/=2;
            for(j = 0; j<n;j++)
                    y[j] = x[j]-t*d[j];
            *a = y_ant;
        }
        else {
            while(*a > aux){
            t/=2;
            for(j = 0; j<n;j++)
                y[j] = x[j]-t*d[j];
            free(a);
            /*free(b);*/
            a = rosenbrook(n,y,0);
            /*b = function(n,x,0);*/
            aux = *b + alpha*t*quad;
            }
        }
        free(x);
        free(b);
        b = a;
        x= malloc(n*sizeof(double));
        for(j = 0; j<n;j++)
            x[j] = y[j];

        maxgrad = max(rosenbrook(n,x,1),n);
        k++;
        printf("Valor da Funcao: %lf; Norma do Gradiente: %lf; t: %lf; k: %d\n",*b,maxgrad,t, k);
        free(d);
    }
    printf("x: [ ");
    for(j=0;j<n;j++)
        printf("%lf ",x[j]);
    printf("]\n");
    free(x);
    free(b);
    tac();
    return 0;
}

double* quadratica(int n, double *x, int flag){
    int i;
    double *f, *g;
    if(flag == 0){
        f = malloc(sizeof(double));
        *f = 0;
        for(i = 0;i<n;i++){
            *f+=0.5*(i+1)*x[i]*x[i];
        }
        return f; 
    }

    if(flag == 1){
        g = malloc(n*sizeof(double));
        for(i = 0;i<n;i++){
            g[i]=(i+1)*x[i];
        }
        return g;
    }
     
}

double* rosenbrook(int n, double *x, int flag){
    int i;
    double *f, *g;
    /*if(flag == 0){
        f = malloc(sizeof(double));
        *f = 0;
        for(i=0;2*i<n;i++){
            *f+=10*(x[2*i+1]-x[2*i]*x[2*i])*(x[2*i+1]-x[2*i]*x[2*i])+(x[2*i]-1)*(x[2*i]-1);
        }
        return f; 
    }*/
    if(flag == 0){
        f = malloc(sizeof(double));
        *f = 0;
        *f = (x[1]-x[0]*x[0])*(x[1]-x[0]*x[0])+(x[0]-1)*(x[0]-1);
        return f; 
    }

    if(flag == 1){
        g = malloc(n*sizeof(double));
        /*printf("Gradiente de f:\n");*/
        g[0] = 2*(x[1]-x[0]*x[0])*(-2*x[0]) + 2*(x[0]-1);
        g[1] = 2*(x[1]-x[0]*x[0]);
        return g;
    }
     
}

double max(double *x, int n){
    double m = 0;
    double mod;
    for(int i = 0; i<n;i++){
        mod = fabs(x[i]);
        if(mod > m)
            m=mod;
    }
    return m;
}