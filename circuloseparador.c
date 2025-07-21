#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

/*Funcoes para calcular o tempo de execucao do codigo*/

clock_t start_time;
void tic() {
    start_time = clock();
}

void tac() {
    clock_t end_time = clock();
    double elapsed_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    printf("Elapsed time: %f seconds\n", elapsed_time);
}

/*Parametros globais*/
#define M 1000
#define N 500
#define pi 3.141592653589793

/*Funcoes que serao utilizadas*/
/*Variavel "n": tamanho do vetor x*/
/*Todas as funcoes possuem uma flag de entrada*/
/*flag == 0 -> retorna funcao; flag == 1 --> retorna gradiente*/

double* fun(int n, double *x,double ai,double bi,int flag);

double* erro(int n,double *x, int flag);

double max(double *x, int n);

/*Variaveis globais de entrada*/
double *a,*b,*v;

int main(){
    tic();
    int n,k,i,j,l;
    double *x,*x_ant,*d,*d_ant,*c,*f_x,*e_teste,*grad_teste, aux, fy_ant;
    double t,epsilon = 1e-5, alpha = 1e-4,sigma=0.5, maxgrad;
    double nx,ng,quad = 0;
    n = 3;
    a = malloc(N*sizeof(double));
    b = malloc(N*sizeof(double));
    v = malloc(N*sizeof(double));
    if(a == NULL || b == NULL || v == NULL)
        exit(1);
    /*Le arquivos de entrada (dados de teste e de treino)*/
    FILE *input = fopen("circulo_treina", "r");
    FILE *teste = fopen("circulo_testa", "r");

    /*Cria gnuplot*/
    FILE *gnuplot = popen("gnuplot", "w");

    /*Cria dois arquivos para cada amostra de dados*/
    /*fora e forateste possuirao as coordenadas dos pontos com indicador positivo (pontos azuis)*/
    /*dentro e dentroteste possuirao as coordenadas dos pontos com indicador negativo (pontos vermelhos)*/
    FILE *fora = fopen("um.txt", "w");
    FILE *dentro = fopen("menosum.txt", "w");
    FILE *forateste = fopen("umteste.txt", "w");
    FILE *dentroteste = fopen("menosumteste.txt", "w");

    for(i=0;i<N;i++){
        fscanf(input,"%lf %lf %lf",&a[i], &b[i], &v[i]);
        if(v[i] == -1.0){
            fprintf(dentro, "%lf %lf\n", a[i], b[i]);
        }
        else
            fprintf(fora, "%lf %lf\n", a[i], b[i]);
    }

    double y[n];
    x = malloc(n*sizeof(double));
    if(x == NULL)
        exit(1);
    for(i = 0; i<n;i++){
        x[i] = 1; /*Condicao inicial*/
    }
    k=0;
    d = erro(n,x,1); /*d aponta para + grad(f(x))*/
    maxgrad = max(d,n);
    f_x = erro(n,x,0); /*f_x aponta para f(x)*/
    d_ant = NULL;
    x_ant = NULL;
    while(maxgrad >= epsilon && k < M){
        if(k == 0){
            t = 1;
        }
        else{ 
        /*Programa recebe x_ant = x^{n-1}, x = x^{n}, d_ant = grad(f(x^{n-1})), d = grad(f(x^{n})), f_x = f(x^{n})*/
        /*Passo espectral*/
            for(i=0;i<n;i++){
            x_ant[i] -=x[i];
            d_ant[i] -= d[i];
            }
            nx=0;
            ng=0;
            for(i=0;i<n;i++){
                nx +=x_ant[i]*x_ant[i];
                ng +=d_ant[i]*d_ant[i];
            }
            t = sqrt(nx)/sqrt(ng);
        }
        quad = 0;
        for(l = 0; l<n;l++){
            quad -=d[l]*d[l]; /* quad = grad(f(x)) \cdot d */
            y[l] = x[l]-t*d[l];
        }
        c = erro(n,y,0); /*c aponta para f(y)*/
        aux = *f_x + alpha*t*quad;
        if(*c <= aux){
            while(*c <= aux){ /*Busca maior tamanho de passo possível no formato 2^n vezes t_inicial*/
                t*=2;
                for(j = 0; j<n;j++)
                    y[j] = x[j]-t*d[j];
                fy_ant = *c;
                free(c);
                c = erro(n,y,0); /*Aponta para o novo f(y)*/
                 aux = *f_x + alpha*t*quad;
            }
            t/=2; /*Resgata o t que satisfaz Armijo*/
            for(j = 0; j<n;j++)
                    y[j] = x[j]-t*d[j];
            *c = fy_ant;
        }
        else {
            while(*c > aux){
                t*=sigma;
                for(j = 0; j<n;j++)
                    y[j] = x[j]-t*d[j];
                free(c);
                c = erro(n,y,0);
                aux = *f_x + alpha*t*quad;
            }
        } 
/*Após as condicionais e seus respectivos loops, obtivemos novos t, y e *c = f(y) que respeitam a condicao de Armijo*/
        free(x_ant);
        x_ant = x;
        free(f_x);
        f_x = c;
        x= malloc(n*sizeof(double));
        for(j = 0; j<n;j++)
            x[j] = y[j];
        free(d_ant);
        d_ant=d;
        d = erro(n,x,1);
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

/*Escrita do codigo do gnuplot*/

    fprintf(gnuplot, "set terminal x11\n");
    fprintf(gnuplot, "set xrange [-13:13]\n");
    fprintf(gnuplot,"set yrange [-13:13]\n");
    fprintf(gnuplot, "set style line 1 pointtype 7 pointsize 1 linecolor rgb '#dd181f'\n");
    fprintf(gnuplot, "set style line 2 pointtype 7 pointsize 1 linecolor rgb '#0060ad'\n");
    fprintf(gnuplot, "set style line 3 \
    linecolor rgb '#ffa500' \
    linetype 1 linewidth 2 \
    pointtype 5 pointsize 1.5\n");
    fprintf(gnuplot, "set key box\n");

/*Para plotar o circulo, duas funcoes foram criadas: uma para o hemisferio superior e outra para o inferior*/

    fprintf(gnuplot,
    "hsuperior(x) = x>=%lf-%lf? x<=%lf+%lf? %lf + sqrt(%lf**2 - (x-%lf)**2) : 1/0 : 1/0\n",x[0],
    x[2],x[0],x[2],x[1],x[2],x[0]);
    fprintf(gnuplot,
    "hinferior(x) = x>=%lf-%lf? x<=%lf+%lf? %lf - sqrt(%lf**2 - (x-%lf)**2) : 1/0 : 1/0\n",x[0],
    x[2],x[0],x[2],x[1],x[2],x[0]);
    
/*Plot do primeiro grafico*/

    fprintf(gnuplot, "plot hsuperior(x) with lines title 'Circulo' linestyle 3\n");
    fprintf(gnuplot, "replot hinferior(x) with lines notitle linestyle 3\n");
    fprintf(gnuplot, "replot 'menosum.txt' with points title 'v = -1' linestyle 1\n");
    fprintf(gnuplot, "replot 'um.txt'  with points title 'v = 1' linestyle 2\n");
    fprintf(gnuplot,"set terminal png giant size 900,600 enhanced\n");
    fprintf(gnuplot,"set termoption enhanced\n");
    fprintf(gnuplot,"set output 'treino.png'\n");
    fprintf(gnuplot,"replot\n"); 

/*Construcao do grafico com os dados de teste*/

    for(i=0;i<N;i++){
        fscanf(teste,"%lf %lf %lf",&a[i], &b[i], &v[i]);
        if(v[i] == -1.0){
            fprintf(dentroteste, "%lf %lf\n", a[i], b[i]);
        }
        else
            fprintf(forateste, "%lf %lf\n", a[i], b[i]);
    }
    e_teste = erro(n,x,0);
    grad_teste = erro(n,x,1);
    printf("Erro no teste: %e; Gradiente no teste: %e\n",*e_teste, *grad_teste);
    free(e_teste);
    free(grad_teste);

/*Plot do novo grafico*/

    fprintf(gnuplot, "plot hsuperior(x) with lines title 'Circulo' linestyle 3\n");
    fprintf(gnuplot, "replot hinferior(x) with lines notitle linestyle 3\n");
    fprintf(gnuplot, "replot 'menosumteste.txt' with points title 'v = -1' linestyle 1\n");
    fprintf(gnuplot, "replot 'umteste.txt'  with points title 'v = 1' linestyle 2\n");
    fprintf(gnuplot, "set output 'teste.png'\n");
    fprintf(gnuplot,"replot\n");
    fflush(gnuplot);
    free(x);
    free(f_x);
    free(d);
    free(a);
    free(b);
    free(v);
    return 0;
}

/*Algoritmos das funcoes auxiliares*/

double* fun(int n, double *x,double ai,double bi,int flag){
    int i,j;
    double *f,r;
    r = (x[0]-ai)*(x[0]-ai)+(x[1]-bi)*(x[1]-bi) - (x[2])*(x[2]);
    if(flag == 0){
        f = malloc(sizeof(double));
        if(f == NULL)
            exit(1);
        *f=2/pi*atan(r);
    }
    else if(flag == 1){
        f = malloc(n*sizeof(double));
        if(f == NULL)
            exit(1);
        r = 1+r*r; 
        f[0] = 4*(x[0]-ai)/(pi*r);
        f[1] = 4*(x[1]-bi)/(pi*r);
        f[2] = -4*(x[2])/(pi*r);
    }
    return f;
}

double* erro(int n,double *x, int flag){
    int i,j,k;
    double *f, *g,*ima,*grad;
    double aux,aux2,aux3;
    if(flag == 0){
        f = malloc(sizeof(double));
        if(f == NULL)
            exit(1);
        *f = 0;
        for(j=0;j<N;j++){
            ima = fun(n,x,a[j],b[j],0);
            *f+=(*ima-v[j])*(*ima-v[j]);
            free(ima);
        }
        *f/=(2*N);
        return f; 
    }

    if(flag == 1){
        g = malloc(n*sizeof(double));
        if(g == NULL)
            exit(1);
        for(i=0;i<n;i++){
            g[i] = 0;
        }
        for(j=0;j<N;j++){
            ima = fun(n,x,a[j],b[j],0);
            grad = fun(n,x,a[j],b[j],1);
            for(i=0;i<n;i++){
                g[i] += (*ima - v[j])*grad[i];
            }
            free(ima);
            free(grad);
        }
        for(i=0;i<n;i++)
            g[i]/=N;
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