#include <omp.h>
int N = 20;
int funcion(double A[N][N], double posiciones[][2])
{
    int k = 0;
    double maximo;
    int num_max, myID;
#pragma omp parallel private(num_max, myID) // ambas son privadas
    {
        myID = omp_get_thread_num();
        num_max = 0;
#pragma omp for private(j) reduction(max : maximo)
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; i < N; j++) // j posible variable compartida
            {
                if (A[i][j] > maximo)
                    maximo = A[i][j];
            }
        }

        #pragma omp for private(j)
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; i < N; j++)
            {
                if (A[i][j] == maximo)
                {
                    #pragma omp critical(pos) // no se puede poner reduction() porque cada hilo hace una fila (las 3 instrucciones son criticas)
                    {
                        posiciones[k][0] = i;
                        posiciones[k][1] = j;
                        k = k + 1;
                    }
                    num_max += 1; // cada hilo interactua con el num_max de manera independiente
                }
            }
        }
        printf("soy el hilo %d y he encontrado %d maximos", myID, num_max);
    }
}

// ejercicio 13
double producto_escalar(double x[], double y[], int n)
{
    // no se puede usar for con reduccion
    int i;
    double suma = 0, suma_local, myID, num_hilos;
#pragma omp parallel private(suma_local, myID, i)
    {
        suma_local = 0;
        // #pragma omp for => APARTADO A
        // for (int i = 0; i < n; i++) => apartado A

        // no se puede usar for
        myID = omp_get_thread_num();
        num_hilos = omp_get_num_threads();

        for (int i = myID; i < n; i += num_hilos)
        {
            suma_local += x[i] * y[i];
        }
#pragma omp atomic // si es atomic, hay que poner +=, susituir suma = suma +...
        suma += suma_local;
    }
    return suma;
}

// ejercicio 2.6
int n;
double a, b[3];
    #pragma omp parallel sections
    {
    #pragma omp section
    {
        //codigo 
    }
    #pragma omp section
    {
        //codigo 2
    }
    #pragma omp section
    {
        //codigo 3
    }


    }

// ejercicio 2.7b
void fun(){
    #pragma omp parallel sections
    {
        #pragma omp section
        {
            //f1(a,b);
            //f2(b,b);
        }
        #pragma omp section
        {
            //f3(c,d);
            //f4(d,d);
        }
        #pragma omp single
        //operacion larga

    }
}