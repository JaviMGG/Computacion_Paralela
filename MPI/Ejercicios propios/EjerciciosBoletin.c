#include <stdio.h>
#include <mpi.h>
#define DIM 1000

/**
 * Ejercicio 2-09
 * Haz una versión paralela MPI del programa anterior, utilizando operaciones de comunicación colectiva cuando sea posible.
 * La función leer deberá ser invocada solo por el proceso 0.
 * Se puede asumir que DIM es divisible entre el número de procesos.
 * Nota: hay que escribir el programa completo, incluyendo la declaración de las variables y las llamadas necesarias para iniciar y cerrar MPI.
 */

void leer(double A[DIM][DIM], double *x)
{
    /*...*/
}

int ejercicio2_09(int argc, char *argv[])
{
    MPI_Init(&arcg, &argv);
    double A[DIM][DIM], Alocal[DIM][DIM], x;
    int i, j, cont, myid, nprocs, contTotal = 0.0;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    if (myid == 0) // solamente el proceso 0 lee
    {
        leer(A, &x);
    }

    MPI_Bcast(&x, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD); // distribuyo los datos leido a todos los procesos

    int k = DIM / nprocs;
    MPI_Scatter(A, k * DIM, MPI_DOUBLE, Alocal, k * DIM, MPI_DOUBLE, 0, MPI_COMM_WORLD); // distribuir por filas

    cont = 0;
    for (i = 0; i < k; i++)
        for (j = 0; j < DIM; j++)
            if (Alocal[i][j] == x)
                cont++;

    MPI_Reduce(&cont, &contTotal, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); // sumo los cont y los guardo en contTotal

    if (myid == 0)
    {
        printf("%d ocurrencias\n", contTotal);
    }

    MPI_Finalize();
    return 0;
}

/**
 * Cuestión 3–1
 * Supongamos definida una matriz de enteros A[M][N].
 * Escribe el fragmento de código necesario para el envío desde P0 y la recepción en P1 de los datos que se especifican en cada caso,
 * usando para ello un solo mensaje.
 * En caso necesario, define un tipo de datos derivado MPI.
 * (a) Envío de la tercera fila de la matriz A.
 * (b) Envío de la tercera columna de la matriz A.
 */
#define A [M][N] // suponer una matriz de enteros definida

void ejercicio3_01()
{
    int myid;
    MPI_Comm_rank(MPI_COMM_WORLD, myid);
    //apartado b
    MPI_Datatype columnas
    MPI_Type_vector(M, 1, N, MPI_INT, &columnas);
    MPI_Commit(&columnas);

    if (myid == 0)
    {
        //apartado a
        MPI_Send(&A[2], N, MPI_INT, 1, 0, MPI_COMM_WORLD);
        //apartado b
        MPI_Send(&A[0][2], 1, columnas, 1, 0, MPI_COMM_WORLD);

    }else if (myid == 1)
    {
        //apartado a
        MPI_Recv(&A[2], N, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_COMM_WORLD);
        //apartado b
        MPI_Recv(&A[0][2], 1, columnas, 0, 0, MPI_COMM_WORLD, MPI_COMM_WORLD);
        
    }
    MPI_Type_free(&columnas);
    
}


/**
 * Sea A un array bidimensional de números reales de doble precisión, de dimensión N×N. Define un tipo de
datos derivado MPI que permita enviar una submatriz de tamaño 3 × 3. Por ejemplo, la submatriz que
empieza en A[0][0] serían los elementos marcados con ⋆:
A =
[⋆ ⋆ ⋆ · · ·]
[⋆ ⋆ ⋆ · · ·]
[⋆ ⋆ ⋆ · · ·]
[· · · · · ·]
[· · · · · ·]
[· · · · · ·]

(a) Realiza las correspondientes llamadas para el envío desde P0 y la recepción en P1 del bloque de la figura.
(b) Indica qué habría que modificar en el código anterior para que el bloque enviado por P0 sea el que
empieza en la posición (0,3), y que se reciba en P1 sobre el bloque que empieza en la posición (3,0)
*/

void ejercicio3_03()
{
    int myid;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Datatype bloque;
    MPI_Type_vector(3, 3, N, MPI_INT, &bloque);
    MPI_Type_commit(&bloque);

    if (myid == 0)
    {
        //apartado a
        MPI_Send(&A[0][0], 1, bloque, 1, 0, MPI_COMM_WORLD);
        //apartado b
        MPI_Send(&A[0][3], 1, bloque, 1, 0, MPI_COMM_WORLD);
    }
    else if(myid == 1)
    {
        //apartado b
        MPI_Recv(&A[0][0], 1, bloque, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        //apartado b
        MPI_Recv(&A[3][0], 1, bloque, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    MPI_Type_free(&bloque);
}