#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

double *substituicaoRegressiva(double **m, int dim)
{
    double *root = (double*)malloc(dim * sizeof(double));
    double sum;
    int i,j,n;

    n = dim - 1;
    root[n] = m[n][dim]/(double)m[n][n];

    for(i = n - 1; i >= 0; i--)
    {
        sum = 0;

        for(j = i + 1; j <= n; j++ )
        {
            sum += m[i][j] * root[j];
        }

        root[i] = (m[i][dim] - sum)/(double)m[i][i];
    }
    return root;
}

void triangularSuperior_p(double **m, int dim)
{
    int i,troca;
    double n;

    double *aux;

	aux = (double*)malloc((dim + 1) * sizeof(double));

    for(i = 0; i < dim; i++)
    {

		troca = -1;
		for(int l = i; l < dim; l++)
		{
			if(m[i][i]*m[i][i] < m[l][i]*m[l][i])
			{
				troca = l;
			}
		}

		if(troca != -1)
		{
			memcpy(aux, m[troca], (dim + 1)*sizeof(double));
			memcpy(m[troca], m[i], (dim + 1)*sizeof(double));
			memcpy(m[i], aux, (dim + 1)*sizeof(double));
		}

        for(int j = i + 1; j < dim; j++)
        {
            n = m[j][i]/(double)m[i][i];

            for(int k = 0; k < dim + 1; k++)
            {
                m[j][k] = m[j][k] - n * m[i][k];
            }
        }
    }
}

double **lerMatrizCompleta(const char *arg, int *dim)
{
    double **m;
    int i,j;

    FILE *arq = fopen(arg,"r");

    if(arq == NULL)
    {
        printf("arquivo nao encontrado\n");
        exit(1);
    }

    fscanf(arq,"%d",dim);

    m = (double**)malloc((*dim) * sizeof(double*));

    for(i = 0; i< *dim; i++)
    {
        m[i] = (double*)malloc((*dim + 1) * sizeof(double));
    }

    for(i = 0; i < *dim; i++)
    {
        for(j = 0; j < *dim + 1; j++)
        {
            fscanf(arq,"%lf", &m[i][j]);
        }
    }
    fclose(arq);
    return m;

}

void imprimeMatrizCompleta(double **m, int dim)
{
    int i,j;

    for(i = 0; i < dim; i++)
    {
        for(j = 0; j < dim + 1; j++)
        {
            printf("%5.2lf\t",m[i][j]);
        }
        puts("");
    }
    puts("________");
}

void imprimeRaiz(double *r, int dim)
{
    int i;

    puts("\n__roots__\n");

    for(i = 0; i < dim; i++)
    {
        printf("x[%d] = %6.5lf\n",i + 1,r[i]);
    }

    puts("________");
}

int main(int argc, char **argv)
{
    double **matriz;
    double *root;
    int dim;

    matriz = lerMatrizCompleta(argv[1],&dim);
    printf("Matriz de entrada\n");
    imprimeMatrizCompleta(matriz,dim);
    printf("pivotamento\n");
    triangularSuperior_p(matriz,dim);
    printf("Matriz de saida\n");
    imprimeMatrizCompleta(matriz,dim);

    root = substituicaoRegressiva(matriz,dim);
    imprimeRaiz(root,dim);

    return 0;
}
