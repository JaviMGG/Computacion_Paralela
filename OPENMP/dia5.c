// a)
#include <stdio.h>
void ejercicio()
{
    int M, N;
    double ej(double x[M], double y[N], double A[M][N])
    {
        int i, j;
        double aux, s = 0.0;
#pragma omp parallel {
#pragma omp for nowait
        for (i = 0; i < M; i++)
            x[i] = x[i] * x[i];
#pragma omp for
        for (i = 0; i < N; i++)
            y[i] = 1.0 + y[i];
#pragma omp for private(j, aux) reduction(+ : s)
        for (i = 0; i < M; i++)
            for (j = 0; j < N; j++)
            {
                aux = x[i] - y[j];
                A[i][j] = aux;
                s += aux;
            }
    }
    return s;
}

/**
(a) Paralelízala eficientemente mediante OpenMP, usando para ello una sola región paralela.
(b) Calcula el número de flops de la función inicial y de la función paralelizada.
(c) Determina el speedup y la eficiencia

    // a)

    #pragma omp parallel for private(j, s)
    for (i = 0; i < n; i++)
    {
        s = 0;
        for (j = 0; j < i; j++)
        {
            s += A[i][j] * b[j];
        }
        c[i] = s;
    #pragma omp atomic
        x[ind[i]] += s;
    }

    // b)

    for (i = 0; i < n; i++)
    {
        s = 0;
    #pragma omp parallel for reduction(+ : s)
        for (j = 0; j < i; j++)
        {
            s += A[i][j] * b[j];
        }
        c[i] = s;
        x[ind[i]] += s;
    }
    // c)
    }
*/
#pragma omp parallel for private(j, s) schedule(static, 1)
for (i = 0; i < n; i++)
{
    s = 0;
    for (j = 0; j < i; j++)
    {
        s += A[i][j] * b[j];
    }
    c[i] = s;
#pragma omp atomic
    x[ind[i]] += s;
}

/*
(a) Realiza una implementación paralela mediante OpenMP, en la que se reparten las iteraciones del
bucle externo.
(b) Realiza una implementación paralela mediante OpenMP, en la que se reparten las iteraciones del
bucle interno.
(c) Para la implementación del apartado (a), indica si cabe esperar que haya diferencias de prestaciones
dependiendo de la planificación empleada. Si es así, ¿qué planificaciones serían mejores y por qué?

*/

// a)
int funcion(int n, double v[])
{
    int i, pos_max = -1;
    double suma, norma, aux, max = -1;
    suma = 0;
#pragma omp parallel
    {
#pragma omp for reduction(+ : suma)
        for (i = 0; i < n; i++)
            suma = suma + v[i] * v[i];
        norma = sqrt(suma);
#pragma omp for
        for (i = 0; i < n; i++)
            v[i] = v[i] / norma;
#pragma omp for private(aux)
        for (i = 0; i < n; i++)
        {
            aux = v[i];
            if (aux < 0)
                aux = -aux;
            if (aux > max)
            {
#pragma omp critical
                if (aux > max)
                {
                    pos_max = i;
                    max = aux;
                }
            }
        }
    }

    return pos_max;
}

// c)

int funcion(int n, double v[])
{
    int i, pos_max = -1;
    double suma, norma, aux, max = -1;
    suma = 0;
#pragma omp parallel
    {
#pragma omp for reduction(+ : suma) schedule(static, 2)
        for (i = 0; i < n; i++)
            suma = suma + v[i] * v[i];
        norma = sqrt(suma);
#pragma omp for schedule(static, 2)
        for (i = 0; i < n; i++)
            v[i] = v[i] / norma;
#pragma omp for private(aux) schedule(static, 2)
        for (i = 0; i < n; i++)
        {
            aux = v[i];
            if (aux < 0)
                aux = -aux;
            if (aux > max)
            {
#pragma omp critical
                if (aux > max)
                {
                    pos_max = i;
                    max = aux;
                }
            }
        }
    }

    return pos_max;
}
/*
(a) Paralelízala con OpenMP, usando una única región paralela.
(b) ¿Tendría sentido poner una cláusula nowait a alguno de los bucles? ¿Por qué? Justifica cada bucle
separadamente.
No, no tendria sentido

(c) ¿Qué añadirías para garantizar que en todos los bucles las iteraciones se reparten de 2 en 2 entre
los hilos?
*/

// a)

double transferencias(double saldos[], int origenes[],
                      int destinos[], double cantidades[], int n)
{
    int i, i1, i2;
    double dinero, maxtransf = 0;
#pragma omp parallel for private(i1, i2, dinero) reduction(max : maxtransf)
    for (i = 0; i < n; i++)
    {
        /* Procesar transferencia i: La cantidad transferida es
         * cantidades[i], que se mueve de la cuenta origenes[i]
         * a la cuenta destinos[i]. Se actualizan los saldos de
         * ambas cuentas y la cantidad maxima */
        i1 = origenes[i];
        i2 = destinos[i];
        dinero = cantidades[i];
#pragma omp atomic
        saldos[i1] -= dinero;
#pragma omp atomic
        saldos[i2] += dinero;
        if (dinero > maxtransf)
        {                // se puede hacer esto o empleando la reduccion de arriba del todo, ambas son validas, ambas no, lo dejo para que se vean
#pragma omp for critical // es mas eficiente hacer la reduccion
            if (dinero > maxtransf)
                maxtransf = dinero;
        }
    }
    return maxtransf;
}

// b)

double transferencias(double saldos[], int origenes[],
                      int destinos[], double cantidades[], int n)
{
    int i, i1, i2;
    double dinero, maxtransf = 0;
#pragma omp parallel for private(i1, i2, dinero)
    for (i = 0; i < n; i++)
    {
        /* Procesar transferencia i: La cantidad transferida es
         * cantidades[i], que se mueve de la cuenta origenes[i]
         * a la cuenta destinos[i]. Se actualizan los saldos de
         * ambas cuentas y la cantidad maxima */
        i1 = origenes[i];
        i2 = destinos[i];
        dinero = cantidades[i];
#pragma omp atomic
        saldos[i1] -= dinero;
#pragma omp atomic
        saldos[i2] += dinero;
        if (dinero > maxtransf)
        {
#pragma omp for critical
            if (dinero > maxtransf)
            {
                maxtransf = dinero;
                ind_transf = i;
            }
        }
    }
    printf("Maxima transferencia: %d\n", ind_transf);
    return maxtransf;
}

/*
(a) Paraleliza la función de forma eficiente mediante OpenMP.
(b) Modifica la solución del apartado anterior para que se imprima el índice de la transferencia con más
dinero.
*/