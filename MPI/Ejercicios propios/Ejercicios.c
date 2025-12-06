#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>

/**
 * 1)
 * Paralelizar el siguiente programa mediante MPI de modo que sea el proceso 0 quien
 * pida un número al usuario y se lo envíe al proceso 1,
 * el cual lo mostrará por pantalla
 * tras recibirlo.
 */
void ejercicio_1(int argc, char **argv)
{
    int myid;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    double numero;
    if (myid == 0)
    {
        printf("Dame un número: ");
        scanf("%lf", &numero);
        MPI_Send(&numero, 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
    }
    else if (myid == 1)
    {
        MPI_Recv(&numero, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        printf("El número proporcionado es el %f\n", numero);
    }
    MPI_Finalize();
    // return 0;
}

/**
 * 2)
 * Dado el siguiente programa, paralelizarlo mediante MPI de manera que el proceso 0
 * inicialice el vector a y se lo envíe al proceso 1, quien lo mostrará por pantalla tras
 * recibirlo.
 */

#define N 5
void ejercicio_2(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    int a[N], i, myid;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    if (myid == 0)
    {
        for (i = 0; i < N; i++)
            a[i] = i;
        MPI_Send(a, N, MPI_INT, 1, 0, MPI_COMM_WORLD);
    }
    else if (myid == 1)
    {
        MPI_Recv(a, N, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        printf("Los elementos del vector son: \n");
        for (i = 0; i < N; i++)
            printf("%d\n", a[i]);
    }
    MPI_Finalize();
    // return 0;
}

/**
 * 3. El siguiente programa suma los n primeros números naturales, siendo dicho valor
 * de n proporcionado como un dato de entrada al programa. Paralelizarlo mediante MPI
 * de manera que el proceso 0 muestre el resultado de la suma por pantalla.
 */

void ejercicio_3(int argc, char *argv[])
{
    int n, i, suma, myid;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    if (argc == 2)
        n = atoi(argv[1]);
    else
        n = 10;
    suma = 0;
    for (i = 1; i <= n; i++)
        suma += i;

    if (myid == 0)
    {
        printf("Valor de la suma = %d\n", suma);
    }
    MPI_Finalize();
    // return 0;
}

/**
 * Desarrollar un programa en el que todos los procesos disponen de un vector de N
 * componentes de tipo entero que envían al proceso 0, el cual almacena
 * consecutivamente todos los elementos que recibe en otro vector. Se supone que el
 * proceso 0 dispone de espacio necesario para almacenar todos los elementos que
 * recibe siempre y cuando el número de procesos sea menor o igual a 32.
*/
#define N = 5;
void ejercicio_5(int v[N])
{
    int myid, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    if(myid == 0){
        if(size <= 32){
            int vNuevo [N];
            for (int i = 1; i < size; i++)
            {
                MPI_Recv( vNuevo, N, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            
        }
    }
    else{
        MPI_Send(v, N, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }
}