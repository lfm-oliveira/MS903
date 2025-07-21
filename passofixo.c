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


#define M 15000
#define N 50
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

double* fun(int n, double *x, int flag);

double u[N],v[N],z[N];

int main(){
    tic();
    int n,k,i,j,l;
    double *x,*d,*a,*b, aux, y_ant;
    double t, epsilon = 1e-4, alpha = 1e-4, maxgrad;
    double quad = 0;
    scanf("%lf",&t);
    n = 3;
    double y[n];
    /*u = malloc(N*sizeof(double));
    v = malloc(N*sizeof(double));
    z = malloc(N*sizeof(double));
    if(u == NULL || v == NULL || z == NULL)
        exit(1);*/
    for(i = 0; i<N;i++){
        scanf("%lf", &u[i]);
        scanf("%lf", &v[i]);
        scanf("%lf", &z[i]);
    }
    x = malloc(n*sizeof(double));
    if(x == NULL)
        exit(1);
    for(i=0;i<3;i++)
        x[i] = 0;
    k=0;
    maxgrad = max(fun(n,x,1),n);
    b = fun(n,x,0); /*b aponta para f(x)*/
    while(maxgrad >= epsilon && k < M){
        /*Calculo da direcao de descida e y para t = 1*/
        d = fun(n,x,1); /*d aponta para + grad f (x)*/
        quad = 0;
        for(l = 0; l<n;l++){
            quad -=d[l]*d[l];
            y[l] = x[l]-t*d[l];
        }
        a = fun(n,y,0); /*a aponta para f(y)*/
        //aux = *b + alpha*quad;
        if(*a < *b){
            for(j = 0; j<n;j++)
                y[j] = x[j]-t*d[j];
            free(a);
            /*free(b);*/
            a = fun(n,y,0);
            /*b = function(n,x,0);*/
            //aux = *b + alpha*t*quad;
        }
        free(x);
        free(b);
        b = a;
        x= malloc(n*sizeof(double));
        for(j = 0; j<n;j++)
            x[j] = y[j];

        maxgrad = max(fun(n,x,1),n);
        k++;
        printf("Valor da Funcao: %lf; Norma do Gradiente: %lf; t: %lf; k: %d\n",*b,maxgrad,t, k);
        //}
        //else{
        //    free(a);
        //}
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


double* fun(int n, double *x, int flag){
    int i,j,p=5;
    double *f, *g;
    double aux1,aux2,aux3;
    if(flag == 0){
        f = malloc(sizeof(double));
        *f = 0;
        for(j = 0; j<N;j++){
            aux1 = exp(x[3]*u[j]) ;
            aux2 = x[2]*aux1;
            aux3 = x[1]+aux2 - v[j];
            *f+=aux3*aux3;
        }
        *f = *f/(2*N);
        return f; 
    }
    if(flag == 1){
        g = malloc(3*sizeof(double));
        /*printf("Gradiente de f:\n");*/
        for(j = 0; j<n; j++)
            g[j] = 0;
        //for(j = 0; j<N; j++){
            for(i = 0;i < p; i++){
                j = rand() % N;
                
                //printf("j: %d\n", j);
                aux1 = exp(x[3]*u[j]) ;
                aux2 = x[2]*aux1;
                aux3 = x[1]+aux2 - v[j];
                g[0]+=aux3;
                g[1] += 2*aux3*aux1;
                g[2] += aux3*u[j]*aux2;
            }
        //}
        g[0]/= p;
        g[1]/= p;
        g[2]/= p;
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