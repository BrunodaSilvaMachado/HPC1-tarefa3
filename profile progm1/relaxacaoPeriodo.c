#include <stdio.h>
#include <string.h>
#include <math.h>

#define N 300

double matriz[N][N];
double mascara[N][N];


void contorno()
{
    int i;

    for(i = 0; i < N;i++)
    {
        matriz[N-1][i] = 1;
    }
    
    for(i = 0; i < N/2; i++)
    {
        mascara[i][N/2] = 1;
    }
}

void relaxacao()
{
    int i,j;

    for(i = 1;i < N - 1; i++)
    {
        for(j = 1; j < N - 1;j++)
        {
            if(mascara[i][j] == 0)
                matriz[i][j] = (matriz[i - 1][j] + matriz[i + 1][j] + matriz[i][j - 1] + matriz[i][j + 1])/4.0;
        }
    }
}

void periodo()
{
    int i;
    
    for(i = 0; i < N; i++)
    {
        matriz[i][0] = matriz[i][N -2];
        matriz[i][N - 1] = matriz[i][1];
    }
}

void imprime()
{
    FILE *arq = fopen("potencialp.txt","w+");
    int i,j;

    for(i = 0;i < N; i++)
    {
        for(j = 0; j < N;j++)
        {
            fprintf(arq,"%lf\t",matriz[i][j]);
        }
        fprintf(arq,"\n");
    }
    fclose(arq);
}

double trs()
{
    int i;
    double sum = 0;

    for(i = 0;i < N; i++)
    {
       sum += matriz[i][i];
    }

    return sum;
}

int main()
{
    int i = 0;
    double s;
    contorno();

    do
    {
        s = trs();
        relaxacao();
        periodo();
        i++;
    }while(fabs(s - trs()) > 1e-5);

    printf("ite: %d\n",i);


    imprime();

    return 0;
}

