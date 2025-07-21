#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

typedef struct dupla {
    double x;
    double y;
} dupla;

typedef struct tripla {
    double erro;
    double norma;
    double* x;
} tripla;

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

double fineura(double *x, int flag,int k);

double* erroneura(int n, double *x, int flag, int N_t);

double max(double *x, int n);

double *errolambda(int n, double*x, double lambda, int flag, int N_t);

tripla* mainneural(int n,dupla* u, double* v,double lambda, double* x,FILE* input, FILE* teste,FILE *gnuplot, int N_t);

dupla *u;
double *v, N = 500, N_teste = 300;
/*Variaveis globais de entrada*/
double *a,*b,*v;
int main(){
    tic();
    int i,n = 9;
    double *xini, lambda = 100,*xteste,*ima;
    tripla* s;
    u = malloc(N*sizeof(dupla));
    v = malloc(N*sizeof(double));
    if(u == NULL || v == NULL)
        exit(1);
    xteste = malloc(n*sizeof(double));
    if(xteste == NULL)
        exit(1);
    /*Le arquivos de entrada (dados de teste e de treino)*/

    /*Cria gnuplot*/
    FILE *gnuplot = popen("gnuplot", "w");

    /*Cria dois arquivos para cada amostra de dados*/
    /*um e umteste possuirao as coordenadas dos pontos com indicador positivo (pontos azuis)*/
    /*menosum e menosumteste possuirao as coordenadas dos pontos com indicador negativo (pontos vermelhos)*/
    /*FILE *um = fopen("um.txt", "w");
    FILE *menosum = fopen("menosum.txt", "w");
    FILE *umteste = fopen("umteste.txt", "w");
    FILE *menosumteste = fopen("menosumteste.txt", "w");
    FILE *umerro = fopen("umerro.txt", "w");
    FILE *menosumerro = fopen("menosumerro.txt", "w");*/

    /*Le arquivos de entrada (dados de teste e de treino)*/
    FILE *input = fopen("dado_treino", "r");
    FILE *teste = fopen("dado_teste", "r");

    for(i=0;i<N;i++){
        fscanf(input,"%lf %lf %lf",&u[i].x, &u[i].y, &v[i]);
    }

    double y[n];
    double m_aux[7][2];
    double l[7];
    xini = malloc(n*sizeof(double));
    if(xini == NULL)
        exit(1);
    printf("x = [");
    for(i = 0; i<n;i++){
        xini[i] = rand() % 30000;
        xini[i]/= 30000;
        xini[i] = 10*(2*(xini[i])-1);
        printf("%lf ",xini[i]);
    }
    printf("]\n");
    s = mainneural(n,u,v,0,xini,input,teste,gnuplot,N);
    xini = s->x;
    free(s);
    for(i=0;i<6;i++){
        s = mainneural(n,u,v,lambda,xini,input,teste,gnuplot,N);
        m_aux[i][0] = s->norma;
        m_aux[i][1] = s->erro;
        l[i] = lambda;
        if(lambda == 0.01){
            for(i=0;i<n;i++)
                xteste[i] = s->x[i];
        }
        lambda = (double) lambda/10.0;
        free(s);
    }
    s = mainneural(n,u,v,0,xini,input,teste,gnuplot,N);
    l[6] = 0;
    m_aux[6][0] = s->norma;
    m_aux[6][1] = s->erro;
    free(s);
    free(u);
    free(v);
    u = malloc(N_teste*sizeof(dupla));
    v = malloc(N_teste*sizeof(double));
    if(u == NULL || v == NULL)
        exit(1);
    for(i=0;i<N_teste;i++){
        fscanf(teste,"%lf %lf %lf",&u[i].x, &u[i].y, &v[i]);
    }
    FILE*  output = fopen("teste_comparado.txt","w");
    
    FILE* gnuplot2 = popen("gnuplot","w");
    fprintf(gnuplot, "set terminal x11\n");
    //fprintf(gnuplot, "set xrange [-13:13]\n");
    //fprintf(gnuplot,"set yrange [-13:13]\n");
    fprintf(gnuplot, "set style line 1 pointtype 7 pointsize 1 linecolor rgb '#dd181f'\n");
    fprintf(gnuplot, "set style line 2 pointtype 7 pointsize 1 linecolor rgb '#0060ad'\n");
    fprintf(gnuplot, "set style line 3 \
    linecolor rgb '#ffa500' \
    linetype 1 linewidth 2 \
    pointtype 7 pointsize 1\n");
    for(i=0;i<7;i++){
        printf("F(x): %lf; Norma de x: %lf, lambda: %lf\n",m_aux[i][1],m_aux[i][0],l[i]);
    }
    for(i=0;i<N_teste;i++){
        ima = fun(xteste,u[i].x,u[i].y,0); 
        fprintf(output,"%lf %lf\n",v[i],*ima);
        free(ima);
    }
    double* erro;
    erro = erroneura(n,xteste,0,N_teste);
    printf("Erro: %lf \n", *erro);
    free(erro);
    fprintf(gnuplot,"n=300 #number of intervals\n");
fprintf(gnuplot,"max=50. #max value\n");
fprintf(gnuplot,"min=0. #min value\n");
fprintf(gnuplot,"width=(max-min)/n #interval width\n");
fprintf(gnuplot,"#function used to map a value to the intervals\n");
fprintf(gnuplot,"hist(x,width)=width*floor(x/width)+width/2.0\n");
fprintf(gnuplot,"set boxwidth width*0.9\n");
fprintf(gnuplot,"set style fill solid 0.5 # fill style\n");
fprintf(gnuplot,"#count and plot\n");
fprintf(gnuplot,"plot 'teste_comparado.txt' u (hist($1,width)):(1.0) smooth freq w boxes lc rgb'green' notitle\n");
fprintf(gnuplot,"set output teste.png\n");
    free(u);
    free(v);
    return 0;
}

tripla* mainneural(int n,dupla* u, double* v,double lambda, double* xini,FILE* input, FILE *teste,FILE *gnuplot,int N_t){
    int i,k,j,l;
    int acertoum, acertomenosum, umtotal, menosumtotal;
    double *x_ant,*d,*d_ant,*c,*f_x,*e_teste,*grad_teste,*ima,*erroini, aux,perc;
    tripla* s;
    double t,epsilon = 1e-5, alpha = 1e-4,sigma=0.5,tol = 1e-4, maxgrad;//1e-6;
    double nx,ng,quad = 0, normaini = 0;
    n = 9;
    /*u = malloc(N*sizeof(dupla));
    v = malloc(N*sizeof(double));
    if(u == NULL || v == NULL)
        exit(1);*/


    /*Cria gnuplot*/
    //FILE *gnuplot = popen("gnuplot", "w");

    /*Cria dois arquivos para cada amostra de dados*/
    /*um e umteste possuirao as coordenadas dos pontos com indicador positivo (pontos azuis)*/
    /*menosum e menosumteste possuirao as coordenadas dos pontos com indicador negativo (pontos vermelhos)*/
    FILE *um = fopen("um.txt", "w");
    FILE *menosum = fopen("menosum.txt", "w");
    FILE *umteste = fopen("umteste.txt", "w");
    FILE *menosumteste = fopen("menosumteste.txt", "w");
    FILE *umerro = fopen("umerro.txt", "w");
    FILE *menosumerro = fopen("menosumerro.txt", "w");
    double *x = malloc(n*sizeof(double));
    double y[n];
    if(x == NULL)
        exit(1);
    for(i=0;i<n;i++)
        x[i] = xini[i];
    k=0;
    d = errolambda(n,x,lambda,1,N_t); /*d aponta para + grad(f(x))*/
    printf("[ ");
    for(i = 0; i<n;i++)
        printf("%e ",d[i]);
    printf("]\n");
    maxgrad = max(d,n);
    f_x = errolambda(n,x,lambda,0,N_t); /*f_x aponta para f(x)*/
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
        c = errolambda(n,y,lambda,0,N_t); /*c aponta para f(y)*/
        aux = *f_x + alpha*t*quad;
        while(*c > aux){
            t*=sigma;
            for(j = 0; j<n;j++)
                y[j] = x[j]-t*(d[j]);
            free(c);
            c = errolambda(n,y,lambda,0,N_t);
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
        d = errolambda(n,x,lambda,1,N_t);
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
    s = malloc(sizeof(tripla));
    if(s == NULL)
        exit(1);
    s->erro = *f_x - 0.5*lambda*nx;
    s->norma = nx;
    s->x = x;
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
    fprintf(gnuplot,"set terminal png giant size 900,600 enhanced\n");
    fprintf(gnuplot,"set termoption enhanced\n");
    /*fprintf(gnuplot,"set output 'treino%.2e.png'\n",lambda);*/
    /*fprintf(gnuplot,"replot\n");*/
    fflush(gnuplot);
    fflush(um);
    fflush(umerro);
    free(f_x);
    free(d);
    free(x);
    return(s);
}


double fineura(double *x, int flag,int k){
    double f, g;
    double a;
    if(k==0){
        f = *x;
        g =1.0;
    } 
    else if(k == 1){
        a = 2/pi;
        f = a*atan(*x);
        g = a/(1+(*x)*(*x));
    }
    else if(k==2){
        f = tanh(*x);
        g = 1-f*f;
    }
    else if(k==3){
        a = sqrt(1+(*x)*(*x));
        f = 0.5*(*x+a);
        g = f/a;
    }
    if(flag == 0){
        return f; 
    }
    else if(flag == 1){
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
    k1 = 3; k2=3;
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
        g[5] = (*aux)*x[7]*(*aux2);
        g[4] = (*aux)*x[6]*(*aux1);
        g[3] = (*aux)*x[7]*(*aux2)*u2;
        g[2] = (*aux)*x[7]*(*aux2)*u1;
        g[1] = (*aux)*x[6]*(*aux1)*u2;
        g[0] = (*aux)*x[6]*(*aux1)*u1;
        free(a1);
        free(a2);
        free(a);
        free(aux);
        free(aux1);
        free(aux2);
        return g;
    }
}


double* erroneura(int n, double *x, int flag, int N_t){
    int i,j,k;
    double *f, *g,*ima,*grad;
    double aux,aux2,aux3,a;
    if(flag == 0){
        f = malloc(sizeof(double));
        if(f == NULL)
            exit(1);
        *f = 0;
        for(j=0;j<N_t;j++){
            ima = fun(x,u[j].x,u[j].y,0);
            *f+=(*ima-v[j])*(*ima-v[j]);
            free(ima);
        }
        *f/=(2*N_t);
        return f; 
    }

    if(flag == 1){
        g = malloc(n*sizeof(double));
        if(g == NULL)
            exit(1);
        for(i=0;i<n;i++){
            g[i] = 0;
        }
        for(j=0;j<N_t;j++){
            ima = fun(x,u[j].x,u[j].y,0);
            grad = fun(x,u[j].x,u[j].y,1);
            for(i=0;i<n;i++){
                g[i] += (*ima - v[j])*grad[i];
            }
            free(ima);
            free(grad);
        }
        for(i=0;i<n;i++)
            g[i]/=N_t;
        return g;
    }
}

double *errolambda(int n, double*x, double lambda, int flag,int N_t){
    double *f,*g,*aux;
    double nx = 0;
    int i;
    if(flag == 0){ 
        f = erroneura(n,x,0,N_t);
        for(i = 0; i < n; i++)
            nx+=x[i]*(x[i]);
        *f+=0.5*lambda*nx;
        return f;
    } else if(flag == 1){
        g = erroneura(n,x,1,N_t);
        for(i = 0; i < n; i++)
            g[i]+=lambda*(x[i]);
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

