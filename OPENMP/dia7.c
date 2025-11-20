#include <omp.h>
void ejercicio209()
{
#pragma omp parallel
    {
#pragma omp sections(s1)
        {
#pragma omp section
//    minx = minimo(x,n); /* T1 */
#pragma omp section
//    maxx = maximo(x,n); /* T2 */
#pragma omp section
            //    calcula_y(y,x,n); /* T4 */
        }

#pragma omp sections(s2)
        {
#pragma omp section
//    calcula_z(z,minx,maxx,n); /* T3 */
#pragma omp section
            //    calcula_x(x,y,n); /* T5 */
        }
    }
    //    calcula_v(v,z,x); /* T6 */
}

void ejercicio210()
{
    double fun1(double a[], int n, double v0)
    {
        int i;
        a[0] = v0;
        // no se puede paralelizar
        for (i = 1; i < n; i++)
            a[i] = genera(a[i - 1], i);
    }

    double compara(double x[], double y[], int n)
    {
        int i;
        double s = 0;
#pragma omp parallel for reduction(+ : s)
        for (i = 0; i < n; i++)
            s += fabs(x[i] - y[i]);
        return s;
    }

    int i, n = 10;
    double a[10], b[10], c[10], x = 5, y = 7, z = 11, w;
#pragma omp parallel
    {
#pragma omp sections
        {
#pragma omp section
            // fun1(a, n, x);        /* T1 */
#pragma omp section
            // fun1(b, n, y);        /* T2 */
#pragma omp section
            // fun1(c, n, z);        /* T3 */
        }

#pragma omp sections
        {
#pragma omp section
            // x = compara(a, b, n); /* T4 */
#pragma omp section
            // y = compara(a, c, n); /* T5 */
#pragma omp section
            // z = compara(c, b, n); /* T6 */
        }
    }
    // w = x + y + z;        /* T7 */
    // printf("w:%f\n", w);
}

int M, N;
double funcion(double A[M][N])
{
    int i, j;
    double suma;
    #pragma omp parallel
    {
        #pragma omp for private(i)
        for (j = 0; j < N; j++)
        {
            for (i = 0; i < M - 1; i++)
            {
                A[i][j] = 2.0 * A[i + 1][j];
            }
        }

        suma = 0.0;
        #pragma omp for private(j) reduction(+:suma)
        {
            for (i = 0; i < M; i++)
            {
                for (j = 0; j < N; j++)
                {
                    suma = suma + A[i][j];
                }
            }
        }
    }
    return suma;
}