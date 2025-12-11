#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

/**
 * Asignacion de tareas
 * ------------------------
 * Determinar en que unidades de procesamiento y en que orden se ejecutará la tarea.
 * ------------------------
 * 
 * Proceso:
        Unidad logica de computo capaz de ejecutar tareas computacionales.
 * Procesador:
        Unidad hardware capaz de realizar operaciones.
 * ------------------------ 
 * Depende de la asignacion, el tiempo de computo puede ser mayor o menor.
 * TIempo de computacion = maximizar la concurrencia.
 * Tiempo de comunicacion = asignar las tareas que tengan relacion entre ellas en el mismo proceso (evita comunicaciones innecesarias).
 * ------------------------ 
 * Asignacion estatica = siempre es identica en cada ejecucion (decisiones se toman antes de la ejecución).
 *      Ventajas:
                No añade ninguna sobrecarga.
                Diseño e implementacion simple.
 *      Desventajas:
                No puede adaptarse a las situaciones cambiantes.
                No puede beneficiarse de la concurrencia.
 * ------------------------
 * Asignacion dinamica = Reparto de trabbajo computacional se realiza en tiempo de ejecucion.
        Ventajas:
                No es necesario conocer el comportamiento a priori del problema.
                Puede beneficiarse de la concurrencia.
        Desventajas:
                Añade sobrecarga.
                Diseño e implementacion mas compleja.
 * ------------------------
 * Agrupamiento:
        Agrupar tareas en unidades logicas.

 * Replicacion:
        Parte de los calculos o datos se repiten en todos los procesos .
        Consiste en copiar datos de acceso comun en los distintos procesos con objetivo reducir las comunicaciones.
 * ------------------------
*/

//ejercicio 2-14
/**
 * coste T1, T2, T3 = 7n (CADA UNA)
 * coste T5, t6 = n
 * coste T4 = 2
 * coste total = 23n + 2 ~ 23n.
 * 
 * Coste paralelo:
    - t(n, 3) = ta(n,3) + tc(n, 3)
    - ta(n,3) = 7n + 2 + n + n ~ 9n + 2
    - tc(n, 3) = 2 * (ts + tw) => 9n + 2(ts + tw)
    - Speedup => t(n, 1) / t(n, 3) = (23n + 2) / (9n + 2(ts + tw)) ==> 2,56
*/
double ejercicio2_14(int val[3])
{
    int myid;
    double a,b,c,d,e,f;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    if(myid == 0){
        a = T1(val[0]);
        aux = a;
        //opcion b
        //MPI_Reduce(&a, &d, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    }
    else if(myid == 1)
    {
        b = T2(val[1]);
        aux = b;
        //opcion b
        //MPI_Reduce(&b, &d, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    }
    else if(myid == 2)
    {
        c = T3(val[2]);
        aux = c;
        //opcion b
        //MPI_Reduce(&c, &d, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    }
    double aux;
    //opcion a
    MPI_Reduce(&aux, &d, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if(myid == 0)
    {
        e = T5(val[2],d);
        f = T6(val[0], val[1], e);
    }
    return f;
}

//ejercicio 1-10
/**
 * Coste T1, T3, T4 = n
 * Coste T2, T5 = 2n
 * Coste T6 = 1
 * Coste total = 3*n + 2*2n +1 +1+1+1 (por operaciones aritmeticas) = 7n + 4 ~ 7n.
 * 
 * COSTE PARALELO asumiendo 2 procesos:
  t(n, 2) = ta(n,2) + tc(n, 2)
  ta(n,2) = 2n + n + 2n + 1 + 1 + 1 ~ 5n + 3
  tc(n,2) = 3(ts + tw) => 5n + 3(ts + tw)
  Speedup => t(n, 1) / t(n, 2) = 7n / (5n + 3(ts + tw))
*/
double ejercicio1_10(int i, int j)
{
    int myid;
    double a,b,c,d,e;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    if (myid == 0)
    {
        a = T1(i);
        MPI_Recv(&b, 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        c = T3(a+b,i);
        MPI_Send(&c, 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
        d = T4(a/c);
        MPI_Recv(&e, 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    else if (myid == 1) //si pidiera mas de 2 procesos, se pondria solo else
    {
        b = T2(j);
        MPI_Send(&b, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        MPI_Recv(&c, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        e = T5(b/c);
        MPI_Send(&e, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }
    return d+e;
}

//ejercicio 3-12
void ejercicio3_12(double A[N][N], int root, MPI_Comm comm)
{
    int myid;
    MPI_Datatype tridiagonal;
    MPI_Comm_rank(comm, &myid);
    MPI_Type_vector(N-2, 3, N+1, MPI_DOUBLE, &tridiagonal);
    MPI_Type_commit(&tridiagonal);
    MPI_Bcast(&A[1][0], 1, tridiagonal, root, comm);
    MPI_Type_free(&tridiagonal);
}

//ejercicio 3-10
void ejercicio3_10(double A[NF][NC], double Aloc[NF][NC/2])
{
    int myid;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Datatype mitad;
    MPI_Type_vector(NF, NC/2, NC, MPI_DOUBLE, &mitad);
    MPI_Type_commit(&mitad);
    if(myid == 0)
    {
        MPI_Sendrecv(&A[0][0], 1, mitad, 0, 0, MPI_COMM_WORLD, &Aloc[0][0], NF*NC/2, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Send(&A[0][NC/2], 1, mitad, 1, 0, MPI_COMM_WORLD);
    }
    else if (myid == 1)
    {
        MPI_Recv(&Aloc[0][0], NF*NC/2, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); //el tipo de dato NO PUEDE EMPLEARSE EN LA RECEPCION
    }
    MPI_Type_free(&mitad);
}

void ejercicioFinal()
{
    int i, n, j;
    double *v,*w,*z, sv, sw, x, res = 0.0;
    leer(&n, &v, &w, &z);
    
    calcula_v(n, v);
    calcula_w(n, w);
    calcula_z(n, z);
    
    for(j = 0; j < n; j++)
    {
        sv = 0;
        for(i = 0; i < n; i++)
        {
            sv += v[i]*w[i];
        }
        for(i = 0; i < n; i++)
        {
            v[i] = sv*v[i];
        }
    }

    //tarea 5
    for(j = 0; j < n; j++)
    {
        sw = 0;
        for(i = 0; i < n; i++)
        {
            sw += v[i]*z[i];
        }
        for(i = 0; i < n; i++)
        {
            v[i] = sw*v[i];
        }
    }

    for(j = 0; j < n; j++)
    {
        x = v[j];
        for(i = 0; i < n; i++)
        {
            v[i] = x*v[i];
        }
    }

    for(j = 0; j < n; j++)
    {
        res += v[j];
    }

    printf("%f\n", res);
}