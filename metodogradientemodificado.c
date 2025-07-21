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
    Escolha d in R^n tal que grad(f(x))^T d < 0 (podemos tomar d = - grad(f(x))
    t = 1
    Enquanto f(x+td) >= f(x), faça t <-- t/2
    x <-- x+td
    k <-- k+1*/
double* quadratica(int n, double *x, int flag);

double* rosenbrook(int n, double *x, int flag);

double max(double *x, int n);

double raiz(double x);

int main(){
    tic();
    int n,k,i,j,l;
    double *x,*x_ant,*d,*d_ant,*a,*b, aux, fy_ant;
    double t,taux, epsilon = 1e-4, alpha = 1e-4,sigma=0.5, maxgrad;
    double nx,ng,quad = 0;
    scanf("%d", &n);
    double y[n];
    x = malloc(n*sizeof(double));
    if(x == NULL)
        exit(1);
    for(i = 0; i<n;i++)
        scanf("%lf", &x[i]);
    k=0;
    d = quadratica(n,x,1); /*d aponta para + grad(f(x))*/
    maxgrad = max(d,n);
    /*free(d);*/
    b = quadratica(n,x,0); /*b aponta para f(x)*/
    d_ant = NULL;
    x_ant = NULL;
    while(maxgrad >= epsilon && k < M){
        if(k == 0){
            t = 1;
        }
        else{ 
        /*Quando o valor de k eh n, o algoritmo esta calculando a iteracao k = n+1 a partir dos dados da anterior)*/ 
        /*Programa recebe x_ant = x^{n-1}, x = x^{n}, d_ant = grad(f(x^{n-1})), d = grad(f(x^{n})), b = f(x^{n})*/
        /*Calculo da direcao de descida e y para t = 1*/
            for(i=0;i<n;i++){
            x_ant[i] -=x[i];
            x_ant[i] *=-1;
            d_ant[i] -= d[i];
            d_ant[i]*=-1;
            }
            nx=0;
            ng=0;
            for(i=0;i<n;i++){
                nx +=x_ant[i]*x_ant[i];
                ng +=d_ant[i]*d_ant[i];
            }
            aux = nx/ng;
            t = sqrt(aux); /* Gradiente espectral na iteracao k+1: t_k = sqrt( ||x^{k}-x^{k-1}|| / ||grad_f(x^{k}) - grad_f(x^{k-1})|| )*/
        }
        quad = 0;
        for(l = 0; l<n;l++){
            quad -=d[l]*d[l]; /* quad = grad(f(x)) \cdot d */
            y[l] = x[l]-t*d[l];
        }
        a = quadratica(n,y,0); /*a aponta para f(y)*/
        aux = *b + alpha*t*quad;
        if(*a <= aux){
            while(*a <= aux){ /*Busca maior tamanho de passo possível no formato 2^n \times t_espectral (ou t_passo0)*/
                t*=2;
                for(j = 0; j<n;j++)
                    y[j] = x[j]-t*d[j];
                fy_ant = *a;
                free(a);
                /*free(b);*/
                a = quadratica(n,y,0); /*Aponta para o novo f(y)*/
                aux = *b + alpha*t*quad;
            }
            t/=2;
            for(j = 0; j<n;j++)
                    y[j] = x[j]-t*d[j];
            *a = fy_ant;
        }
        else {
            while(*a > aux){
                t*=sigma;
                for(j = 0; j<n;j++)
                    y[j] = x[j]-t*d[j];
                free(a);
                /*free(b);*/
                a = quadratica(n,y,0);
                /*b = function(n,x,0);*/
                aux = *b + alpha*t*quad;
            }
        } /*Após as condicionais e seus respectivos loops, obtivemos novos t, y e *a = f(y) que respeitam a condicao de Armijo*/
        free(x_ant);
        x_ant = x;
        free(b);
        b = a;
        x= malloc(n*sizeof(double));
        for(j = 0; j<n;j++)
            x[j] = y[j];
        free(d_ant);
        d_ant=d;
        d = quadratica(n,x,1);
        maxgrad = max(d,n);
        k++;
        printf("Valor da Funcao: %lf; Norma do Gradiente: %lf; t: %lf; k: %d\n",*b,maxgrad,t, k);
        
    } /*Ao final, teremos as seguintes variaveis alocadas ativas: x_ant, x, d_ant, d, b (= a)*/
    free(x_ant);
    free(d_ant);
    printf("x: [ ");
    for(j=0;j<n;j++)
        printf("%lf ",x[j]);
    printf("]\n");
    free(x);
    free(b);
    free(d);
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

