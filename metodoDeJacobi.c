#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

void imprimeMatrizCompleta(double **m, int dim);

double norm(double *v,int dim)
{
	int i;
	double sum = 0;
	for(i = 0; i < dim; i++)
	{
		sum += v[i] * v[i];
	}
	
	sum = sqrt(sum);
	
	return sum;
}

double *metodoJacobi(double **A, int linha, int coluna)
{
	double *y;
	double *x;
	
	y = (double*)malloc(coluna*sizeof(double));
	x = (double*)malloc(coluna*sizeof(double));
	
	memset(y,0,coluna*sizeof(double));

	double E = 0.00001;//Precisão
	int m = 1000;//Numero máximo de iterações
	int ni = 0;//Contador de iterações
	int continuar = 1;

	while (continuar && ni < m) {
	    for (int i=0; i < linha; i++) {
	        double soma = 0;
	        for (int j = 0; j < coluna; j++) {
	            if (j != i) {
	                soma = soma + A[i][j]*y[j]/A[i][i];
	            }
	            x[i] = (A[i][coluna]/A[i][i]) - soma;
	        }
		}
	    if (fabs(norm(x,coluna) - norm(y,coluna)) < E) {
	        continuar = 0;
		} else {
	        memcpy(y,x,linha*coluna*sizeof(double));
		}
	    ni = ni + 1;
	}
	return y;
}

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

int triangularSuperior_p(double **m, int dim)
{
    int i,j,k,l,troca,cont = 0;
    double n;

    double *aux;

	aux = (double*)malloc((dim + 1) * sizeof(double));

    for(i = 0; i < dim; i++)
    {

		troca = -1;
		for(l = i; l < dim; l++)
		{
			if(m[i][i] < fabs(m[l][i]))
			{
				troca = l;
			}
		}

		if(troca != -1)
		{
			memcpy(aux, m[troca], (dim + 1) * sizeof(double));
			memcpy(m[troca], m[i], (dim + 1) * sizeof(double));
			memcpy(m[i], aux, (dim + 1) * sizeof(double));
			cont++;
		}
		imprimeMatrizCompleta(m,dim);

        for(j = i + 1; j < dim; j++)
        {
            n = m[j][i]/(double)m[i][i];

            for(k = 0; k < dim + 1; k++)
            {
                m[j][k] = m[j][k] - n * m[i][k];
            }
        }
    }
    return cont;
}


void triangularSuperior(double **m, int dim)
{
    int i,j,k;
    double n;

    for(i = 0; i < dim; i++)
    {
        for(j = i + 1; j < dim; j++)
        {
            n = m[j][i]/(double)m[i][i];

            for(k = 0; k < dim + 1; k++)
            {
                m[j][k] = m[j][k] - n * m[i][k];
            }
        }
    }
}

double determinante(double **m, int dim,int t)
{
	double det = 1;
	int i;
	
	for(i = 0; i < dim; i++)
	{
		det *= m[i][i];
	}
	
	return pow(-1,t) * det;
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
            (j == dim)? printf("| %6.3lf\t",m[i][j]): printf("%6.3lf\t",m[i][j]);
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
    int dim,trocas;

    matriz = lerMatrizCompleta(argv[1],&dim);
    printf("Matriz de entrada\n");
    imprimeMatrizCompleta(matriz,dim);
    printf("pivotamento\n");
    trocas = triangularSuperior_p(matriz,dim);
    printf("Matriz de saida\n");
    imprimeMatrizCompleta(matriz,dim);

    root = substituicaoRegressiva(matriz,dim);
    imprimeRaiz(root,dim);
	imprimeRaiz(metodoJacobi(matriz,dim,dim),dim);
    printf("determinate:\t%lf\n",determinante(matriz,dim,trocas));
	
	for(;dim > 0;dim--)
		free(matriz[dim - 1]);
	free(root);

    return 0;
}
