#include <stdio.h>
#include <stdlib.h>

double *procuraautovalor(double ** matriz, int n);

int main(){

    return 0;
}

double *procuraautovalor(double ** matriz,int n){
    /*procurar autovalores a partir dos discos de Gershgorin (centro + raio de maior modulo e centro - raio de menor modulo)*/
    int i,j;
    double centro,raio = 0;
    for(i=0;i<n;i++){
        centro = matriz[i][i];
        for(j=0;j<n;j++){
            raio+=dabs(matriz[i][j]);
            if(j==i)
                j++;
        }
    }
}