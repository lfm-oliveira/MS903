#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

typedef struct dupla {
    double x;
    double y;
} dupla;
clock_t start_time;

void tic() {
    start_time = clock();
}

void tac() {
    clock_t end_time = clock();
    double elapsed_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    printf("Elapsed time: %f seconds\n", elapsed_time);
}


#define M 500
#define pi 3.141592653589793
/*Dado x in R^n, epsilon > 0, M > 0, 0 < sigma < 1
k = 0
Enquanto ||grad(f(x))|| >= epsilon e k < M, faça:
    Escolha d in R^n tal que grad(f(x))^t d < 0 (podemos tomar d = - grad(f(x))
    t =1
    Enquanto f(x+td) >= f(x), faça t <- sigma*t
    x<-x+td
    k<-k+1*/

double* fun(double *x, double u1, double u2, int flag);

double fineura(double *x, int flag,int ind);

double* erroneura(int n, double *x, int flag);

double max(double *x, int n);

double *errolambda(int n, double*x, double lambda, int flag);

double sqrtmine(double x);

double atanmine(double x);

double expmine(double x);

double tanhmine(double x);

dupla *u;
double *v, N = 2000 ;
/*Variaveis globais de entrada*/
double *a,*b,*v;
int main(){
    tic();
    int n,k,i,j,l;
    int acertoum, acertomenosum, umtotal, menosumtotal;
    double perc;
    double *x,*x_ant,*d,*d_ant,*c,*f_x,*e_teste,*grad_teste,*ima, aux, fy_ant;
    double t,epsilon = 1e-5, alpha = 1e-4,sigma=0.5,tol = 1e-4, maxgrad,lambda=0;//1e-6;
    double nx,ng,quad = 0;
    n = 9;
    u = malloc(N*sizeof(dupla));
    v = malloc(N*sizeof(double));
    if(u == NULL || v == NULL)
        exit(1);
    /*Le arquivos de entrada (dados de teste e de treino)*/
    FILE *input = fopen("curva_treina", "r");
    FILE *teste = fopen("curva_testa", "r");

    /*Cria gnuplot*/
    FILE *gnuplot = popen("gnuplot", "w");

    /*Cria dois arquivos para cada amostra de dados*/
    /*um e umteste possuirao as coordenadas dos pontos com indicador positivo (pontos azuis)*/
    /*menosum e menosumteste possuirao as coordenadas dos pontos com indicador negativo (pontos vermelhos)*/
    FILE *um = fopen("um.txt", "w");
    FILE *menosum = fopen("menosum.txt", "w");
    FILE *umteste = fopen("umteste.txt", "w");
    FILE *menosumteste = fopen("menosumteste.txt", "w");
    FILE *umerro = fopen("umerro.txt", "w");
    FILE *menosumerro = fopen("menosumerro.txt", "w");

    for(i=0;i<N;i++){
        fscanf(input,"%lf %lf %lf",&u[i].x, &u[i].y, &v[i]);
    }

    double y[n];
    x = malloc(n*sizeof(double));
    if(x == NULL)
        exit(1);
    srand( (unsigned)time(NULL) );
    printf("x = [ ");
    for(i = 0; i<n;i++){
        x[i] = (rand() % 30000);
        x[i]/= 30000;
        x[i] = 10*(2*(x[i])-1);
        printf("%lf ",x[i]);
    }
    printf("]\n");
    k=0;
    d = errolambda(n,x,lambda,1); /*d aponta para + grad(f(x))*/
    printf("[ ");
    for(i = 0; i<n;i++)
        printf("%e ",d[i]);
    printf("]\n");
    maxgrad = max(d,n);
    f_x = errolambda(n,x,lambda,0); /*f_x aponta para f(x)*/
    d_ant = NULL;
    x_ant = NULL;
    while(maxgrad >= epsilon && k < M){
        if(k == 0){
            t = 1;
        }
        else { 
        /*Programa recebe x_ant = x^{n-1}, x = x^{n}, d_ant = grad(f(x^{n-1})), d = grad(f(x^{n})), f_x = f(x^{n})*/
        /*Passo espectral*/
            for(i=0;i<n;i++){
            x_ant[i] -=x[i];
            d_ant[i] -= d[i];
            }
            nx=0;
            ng=0;
            for(j=0;j<n;j++){
                nx +=x_ant[j]*(x_ant[j]);
                ng +=d_ant[j]*(d_ant[j]);
            }
            t = sqrt(nx)/sqrt(ng);
        }
        quad = 0;
        for(l = 0; l<n;l++){
            quad -=d[l]*(d[l]); /* quad = grad(f(x)) \cdot d */
            y[l] = x[l]-t*(d[l]);
        }
        c = errolambda(n,y,lambda,0); /*c aponta para f(y)*/
        aux = *f_x + alpha*t*quad;
        while(*c > aux){
            t*=sigma;
            for(j = 0; j<n;j++)
                y[j] = x[j]-t*(d[j]);
            free(c);
            c = errolambda(n,y,lambda,0);
            aux = *f_x + alpha*t*quad;
        }
         
/*Após o loop, obtivemos novos t, y e *c = f(y) que respeitam a condicao de Armijo*/
        free(x_ant);
        x_ant = x;
        free(f_x);
        f_x = c;
        x= malloc(n*sizeof(double));
        for(j = 0; j<n;j++)
            x[j] = y[j];
        free(d_ant);
        d_ant=d;
        d = errolambda(n,x,lambda,1);
        maxgrad = max(d,n);
        k++;
        printf("Valor da Funcao: %e; Norma do Gradiente: %e; t: %e; k: %d\n",*f_x,maxgrad,t, k);
    } /*Ao final, teremos as seguintes variaveis alocadas ativas: x_ant, x, d_ant, d, f_x (= c)*/
    free(x_ant);
    free(d_ant);
    printf("x: [ ");
    for(j=0;j<n;j++)
        printf("%lf ",x[j]);
    printf("]\n");
    tac();

    nx=0;
    for(i=0;i<n;i++){
        nx+=x[i]*(x[i]);
    }
/*Escrita do codigo do gnuplot*/

    fprintf(gnuplot, "set terminal x11\n");
    //fprintf(gnuplot, "set xrange [-13:13]\n");
    //fprintf(gnuplot,"set yrange [-13:13]\n");
    fprintf(gnuplot, "set style line 1 pointtype 7 pointsize 1 linecolor rgb '#dd181f'\n");
    fprintf(gnuplot, "set style line 2 pointtype 7 pointsize 1 linecolor rgb '#0060ad'\n");
    fprintf(gnuplot, "set style line 3 \
    linecolor rgb '#ffa500' \
    linetype 1 linewidth 2 \
    pointtype 7 pointsize 1\n");
    fprintf(gnuplot, "set style line 4 pointtype 7 pointsize 1 linecolor rgb '#800080'\n");
    fprintf(gnuplot, "set key box\n");
    /*fprintf(gnuplot,
    "fione(x) = 2*atan(x)/pi");
    fprintf(gnuplot,
    "fitwo(x) = tanh(x)");*/
    menosumtotal = umtotal = 0;
    acertomenosum = acertoum = 0;
    for(i=0;i<N;i++){
        ima = fun(x, u[i].x, u[i].y,0);
        if(v[i] == -1.0){
            menosumtotal++; 
            if(*ima < 0){
                fprintf(menosum, "%lf %lf %lf\n", u[i].x, u[i].y, *ima);
                acertomenosum++;  
            }
            else {
                fprintf(menosumerro, "%lf %lf %lf\n", u[i].x, u[i].y,*ima);
            }
        }
        else if(v[i] == 1.0) { 
            umtotal++;
            if(*ima > 0){
                fprintf(um, "%lf %lf %lf\n", u[i].x, u[i].y, *ima);
                acertoum++;
            }
            else
                fprintf(umerro, "%lf %lf %lf\n", u[i].x, u[i].y,*ima);
        }
        free(ima);
    }

    perc = (acertoum + acertomenosum)/(double) N;
    fprintf(gnuplot, "set title 'Treino com k1=3 e k2=2 (Acertos (1) = %d/%d | Acertos(-1) = %d/%d) | Total: %lf' font ',26'\n",acertoum,umtotal,acertomenosum,menosumtotal,perc);
    fprintf(gnuplot, "plot 'um.txt'  with points title 'v = 1' linestyle 2\n");
    fprintf(gnuplot, "replot 'menosum.txt' with points title 'v = -1' linestyle 1\n");
    fprintf(gnuplot, "replot 'menosumerro.txt'  with points title 'v = -1 (Erro)' linestyle 3\n");
    fprintf(gnuplot, "replot 'umerro.txt'  with points title 'v = 1 (Erro)' linestyle 4\n");
    fprintf(gnuplot,"set terminal png giant size 900,600 enhanced\n");
    fprintf(gnuplot,"set termoption enhanced\n");
    fprintf(gnuplot,"set output 'treino.png'\n");
    fprintf(gnuplot,"replot\n");
    fflush(gnuplot);
    free(x);
    free(f_x);
    free(d);
    free(u);
    free(v);
    return 0;
}

double fineura(double *x, int flag,int ind){
    double f, g;
    double a;
    if(ind==0){
        f = *x;
        g =1.0;
    } 
    else if(ind == 1) {
        a = 2/pi;
        f = a*atan(*x);
        g = a/ (double) (1+(*x)*(*x));
    }
    else if(ind==2) {
        f = tanh(*x);
        g = 1-f*f;
    }
    else if(ind==3) {
        a = sqrt(1+(*x)*(*x));
        f = 0.5*(*x+a);
        g = f/(double) a;
    }
    if(flag == 0){
        return f; 
    }
    else if(flag == 1) {
        return g;
    }
}

double* fun(double *x,double u1, double u2, int flag){ /*(u1,u2) = (x,y)*/
    double *f, *g;
    double *a,*a1,*a2,*aux,*aux1,*aux2;
    int k1,k2;
    //k1 = 1; k2=1;
    //k1 = 2; k2 = 2; /*(Modificacao)*/
    //k1 = 1; k2=2; /*(Modificacao)*/
    k1 = 1; k2=3;
    a1 = malloc(sizeof(double));
    a2 = malloc(sizeof(double));
    a = malloc(sizeof(double));
    if(a1 == NULL || a2 == NULL || a == NULL)
        exit(1);
    *a1 = x[0]*u1+x[1]*u2+x[4];
    *a2 = x[2]*u1+x[3]*u2+x[5];
    *a = x[6]*fineura(a1,0,k1)+x[7]*fineura(a2,0,k1)+x[8];
    if(flag == 0){
        f = malloc(sizeof(double));
        *f = fineura(a,0,k2);
        if(f == NULL)
            exit(1);
        free(a1);
        free(a2);
        free(a);
        return f; 
    }

    if(flag == 1){
        g = malloc(9*sizeof(double));
        aux = malloc(sizeof(double));
        aux1 = malloc(sizeof(double));
        aux2 = malloc(sizeof(double));
        if(g == NULL || aux == NULL || aux1 == NULL || aux2 == NULL)
            exit(1);
        *aux = fineura(a,1,k2);
        *aux1 = fineura(a1,1,k1);
        *aux2 = fineura(a2,1,k1);
        g[8] = *aux;  
        g[7] = (*aux)*fineura(a2,0,k1);
        g[6] = (*aux)*fineura(a1,0,k1);
        g[5] = (*aux)*(x[7])*(*aux2);
        g[4] = (*aux)*(x[6])*(*aux1);
        g[3] = (*aux)*(x[7])*(*aux2)*u2;
        g[2] = (*aux)*(x[7])*(*aux2)*u1;
        g[1] = (*aux)*(x[6])*(*aux1)*u2;
        g[0] = (*aux)*(x[6])*(*aux1)*u1;
        free(a1);
        free(a2);
        free(a);
        free(aux);
        free(aux1);
        free(aux2);
        return g;
    }
}


double* erroneura(int n, double *x, int flag){
    int i,j;
    double *f, *g,*ima,*grad;
    double aux,aux2,aux3,a;
    if(flag == 0){
        f = malloc(sizeof(double));
        if(f == NULL)
            exit(1);
        aux = 0;
        for(j=0;j<N;j++){
            ima = fun(x,u[j].x,u[j].y,0);
            aux+=(*ima-v[j])*(*ima-v[j]);
            free(ima);
        }
        *f=aux/(2.0*N);
        return f; 
    }

    if(flag == 1){
        g = malloc(n*sizeof(double));
        if(g == NULL)
            exit(1);
        for(i = 0; i< n; i++){
            aux = 0;
            for(j = 0;j<n;j++){
                ima = fun(x,u[j].x,u[j].y,0);
                grad = fun(x,u[j].x,u[j].y,1);
                aux+=(*ima-v[j])*(grad[i]);
                free(ima);
                free(grad);
            }
            g[i] = aux/(double) N;
        }
        return g;
    }
}

double *errolambda(int n, double*x, double lambda, int flag){
    double *f,*g,*aux;
    double nx = 0, aux1;
    int i;
    if(flag == 0) { 
        f = erroneura(n,x,0);
        aux1 = *f;
        for(i = 0; i < n; i++)
            nx+=x[i]*(x[i]);
        aux1+=0.5*lambda*nx;
        *f = aux1;
        return f;
    } else if(flag == 1) {
        g = erroneura(n,x,1);
        for(i = 0; i < n; i++)
            g[i]=g[i]+lambda*(x[i]);
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

double atanmine(double x){
    double soma;
    double aux, sign = 1;
    int i;
    if(fabs(x) > 1){
        if(x < 0)
            sign = -1;
        soma = sign*pi/2.0;
        aux = 1/x;
        for(i =0; i<100;i++){
            soma-=aux/(2*i+1.0);
            aux*=-1/(x*x);
        }
        soma-=aux/(2*i+1.0);
    }
    else {
        soma = 0;aux = x;
        for(i =0; i<100;i++){
            soma+=aux/(2*i+1.0);
            aux*=-1*x*x;
        }
        soma+=aux/(2*i+1.0);
    }
    return soma;
}

double expmine(double x){
    double soma = 1;
    double aux = 1.0,aux2 = x,prod= 1.0;
    int total = 7000,i;
    for(i = 1; fabs(aux2)/prod >= 1e-6; i++){
        aux2*=x;
        prod*=(i+1.0);
    }
    if(i < x)
        i = abs(x)+2; /*a fim de que proximos termos da serie sejam estritamente decrescentes*/
    if(i > 7000)
        total = i;
    for(i = 0; i<total;i++){
        aux*=x/(i+1.0);
        soma+=aux;
    }
    return soma;
}

double tanhmine(double x){
    double soma = 1,aux = -1;
    int sign = 1;
    if(x < 0)
        sign = -1;
    for(int i =0; i<100;i++){
        aux *=-1;
        soma-=2*aux*expmine(-2*sign*(i+1)*x);
    }
    soma*=sign;
    return soma;
}

double sqrtmine(double x){
    double an,bn,sq;
    int i;
    if(x == 0)
        return 0;
    for(i = 0; i*i<x;i++);
    sq = i-1;
    if(sq == 0)
        sq = 1e-2;
    while(fabs(sq*sq - x) > 1e-8){
        an = (x - sq*sq)/(2.0*sq);
        bn = sq + an;
        sq = bn - an*an/(2.0*bn);
    }
    return sq;
}

