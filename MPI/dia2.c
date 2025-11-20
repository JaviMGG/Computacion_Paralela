#include <mpi.h>

/**
 * Modo de Envio estándar:
 * MPI_Send(buf, count, datatype, dest, tag, comm)
 *      buf         = buffer (direccion de memoria de la variable donde queremos enviar los datos, indicamos solo posicion inicial)
 *      count       = nº datos que enviamos / longitud del vector que enviamos
 *      datatype    = tipo de dato que enviamos (entero, float, char...)
 *      dest        = proceso destino
 *      tag         = etiqueta para los mensajes para poder diferenciarlos (tipo entero)
 *      comm        = comunicador donde se situa el envio
 *
 * MPI_Ssend(buf, count, datatype, dest, tag, comm) //similar a lo anterios pero es sincronico
 * MPI_Bsend(buf, count, datatype, dest, tag, comm) //similar a lo anterior pero usa un buffer
 *
 * MPI_ISend(buf, count, datatype, dest, tag, comm)
 *      buf         = buffer (direccion de memoria de la variable donde queremos enviar los datos, indicamos solo posicion inicial)
 *      count       = nº datos que enviamos / longitud del vector que enviamos
 *      datatype    = tipo de dato que enviamos (entero, float, char...)
 *      dest        = proceso destino
 *      tag         = etiqueta para los mensajes para poder diferenciarlos (tipo entero)
 *      comm        = comunicador donde se situa el envio
 *
 * caundo se ejecuta, el emisor no se bloquea, continua en ejecucion
 */

/**
 * ------------------------------------
 * TOPOLOGIA EN MALLA
 * ------------------------------------
 * P0 ---> P1 ---> P2 ---> P3 ---> ... ---> Pn
 *
 * if (rank == 0){                                                  //SI SOMOS EL PRIMER PROCESO
 *      MPI_Send(&val, 1, MPI_DOUBLE, rank+1, 0, comm);
 * }else if (rank == p-1){                                          //SI SOMOS EL ULTIMO PROCESO
 *      MPI_Recv(&val, 1, MPI_DOUBLE, rank-1, 0, comm, &status);
 * }else{                                                           //SI SOMOS EL PROCESO DEL MEDIO
 *      MPI_Send(&val, 1, MPI_DOUBLE, rank+1, 0, comm);             //SI EL PROCESO ULTIMO PIDE RECIBIR, EL PROCESO PENULTIMO ENVIA
 *      MPI_Recv(&val, 1, MPI_DOUBLE, rank-1, 0, comm, &status);    //UNA VEZ HA ENVIADO, SE PONE A RECIBIR LO QUE HACE QUE EL PROCESO ANTEPENULTIMO LE ENVIE Y ASI HASTA EL INICIO
 * }
 *
 * ------------------------------------
 * Posible modificacion del codigo
 * ------------------------------------
 * if (rank == 0) //SI SOMOS EL PRIMER PROCESO
 *      prev = MPI_PROC_NULL;
 * else prev = rank-1:
 *
 * if (rank == p-1)
 *      next = MPI_PROC_NULL;
 *
 * else next = rank+1;
 *
 * MPI_Send(&val, 1, MPI_DOUBLE, rank+1, 0, comm);
 * MPI_Recv(&val, 1, MPI_DOUBLE, rank-1, 0, comm, &status);
 * A NIVEL DE PRESTACIONES, NO SUPONE NINGUNA MODIFICACION (SI FUERAN CON BUFFER Y CON POCOS DATOS, EL PROCESO NO SE QUEDARIA ESPERANDO EN EL SEND)
 *
 * ------------------------------------
 * Modificar tipo de envio (ahora no es secuencial, se hace por bloques de pares e impares)
 * ------------------------------------
 * if (rank == 0) //SI SOMOS EL PRIMER PROCESO
 *      prev = MPI_PROC_NULL;
 * else prev = rank-1:
 *
 * if (rank == p-1)
 *      next = MPI_PROC_NULL;
 *
 * else next = rank+1;
 *
 * if(rank%2 == 0){//los procesos pares se ponen a enviar y luego a recibir
 *      MPI_Send(&val, 1, MPI_DOUBLE, rank+1, 0, comm);
 *      MPI_Recv(&val, 1, MPI_DOUBLE, rank-1, 0, comm, &status);
 * }else{//los impares se ponen a recibir y luego a enviar
 *      MPI_Recv(&tmp, 1, MPI_DOUBLE, rank-1, 0, comm, &status);
 *      MPI_Send(&val, 1, MPI_DOUBLE, rank+1, 0, comm);
 *      val = tmp;
 * }
 *
 * ------------------------------------
 * Hacerlo mas "eficiente" (usando la funcion MPI_Sendrecv_replace)
 * ------------------------------------
 *
 * ------------------------------------
 * Desplazamiento en anillo
 * ------------------------------------
 * if (rank == 0) //SI SOMOS EL PRIMER PROCESO
 *      prev = p-1;
 * else prev = rank-1:
 *
 * if (rank == p-1)
 *      next = 0;
 * else next = rank+1;
 *
 * MPI_Send(&val, 1, MPI_DOUBLE, rank+1, 0, comm);
 * MPI_Recv(&val, 1, MPI_DOUBLE, rank-1, 0, comm, &status);
*/

// ejercicio 1-11
int N = 0;
void intercambiar(double x[N], int proc1, int proc2)
{
    int myid, i;
    double y[N];
    MPI_Comm_rank(MPI_COMM_WORLD, &myid); // devuelve el id de mi porceso
    if (myid == proc1)
    {
        MPI_Send(x, N, MPI_DOUBLE, proc2, 0, MPI_COMM_WORLD);
        MPI_Recv(x, N, MPI_DOUBLE, proc2, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    else if (myid == proc2)
    {
        MPI_Recv(y, N, MPI_DOUBLE, proc1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Send(x, N, MPI_DOUBLE, proc1, 0, MPI_COMM_WORLD);
        for (i = 0; i < N; i++)
        {
            x[i] = y[i];
        }
        memcpy(x, y, N * sizeof(double));
        // copia desde la direccion de memoria de Y hasta N*sizeof(double) al vector X desde direccion inicial hasta el N*sizeof(double)
    }
}

// ejercicio 1-7
void ejercicio1_7()
{
    MPI_status stat;
    MPI_Request request;
    int sbuf[N], rbuf[N], rank, size, src, dst;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    src = (rank == 0) ? size - 1 : rank - 1;
    dst = (rank == size - 1) ? 0 : rank + 1;

    /*solucion 1*/
    MPI_Bsend(sbuf, N, MPI_INT, dst, 111, MPI_COMM_WORLD);
    MPI_Recv(rbuf, N, MPI_INT, src, 111, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    /*solucion 2*/
    if (rank % 2 == 0)
    {
        MPI_Ssend(sbuf, N, MPI_INT, dst, 111, MPI_COMM_WORLD);
        MPI_Recv(rbuf, N, MPI_INT, src, 111, MPI_COMM_WORLD, &stat);
    }
    else
    {
        MPI_Recv(rbuf, N, MPI_INT, src, 111, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Ssend(sbuf, N, MPI_INT, dst, 111, MPI_COMM_WORLD);
    }

    /*solucion 3*/
    MPI_Sendrecv(sbuf, N, MPI_INT, dst, 111, rbuf, N, MPI_INT, src, 111, MPI_COMM_WORLD, &stat);

    /*solucion 4*/
    MPI_Isend(sbuf, N, MPI_INT, dst, 111, MPI_COMM_WORLD, &request);
    MPI_Recv(rbuf, N, MPI_INT, src, 111, MPI_COMM_WORLD, &stat);
    MPI_Wait(&request, MPI_STATUS_IGNORE);
}

// ejercicios 3

#include <stdio.h>
#include <stdlib.h>
void ejercicio_3(int argc, char *argv[])
{
    int n, i, suma, myid, nprocesos, suma_local;
    MPI_Init(&argc, &argv);

    if (argc == 2)
        n = atoi(argv[1]);
    else
        n = 10;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocesos);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    suma_local = 0;
    for (i = myid + 1; i <= n; i += nprocesos)
        suma_local += i;

    // DESDE AQUI
    if (myid == 0)
    {
        suma = suma_local;
        for (int i = 1; i < nprocesos; i++)
        {
            MPI_Recv(&suma_local, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // MPI_ANY_SOURCE SIEMPRE ES EN RECEPCION, NUNCA SE PUEDE PONER EN ENVIO
            suma += suma_local;
        }
    }
    else
    {
        MPI_Send(&suma_local, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }
    // HASTA AQUI ES EQUIVALENTE A UN REDUCTION DE OPENMP

    if (myid == 0)
    {
        printf("Valor de la suma = %d\n", suma);
    }

    MPI_Finalize();
    // return 0;
}

/**
 * Operaciones de comunicacion colectiva:
 *
 * Sincronizacion   (Barrier)
 * Difusion         (Bcast)
 * Reparto          (Scatter)
 * Recogida         (Gather)
 * Multi-recogida   (Allgather)
 * Todos a todos    (Alltoall)
 * Reduccion        (Reduce)
 * Prefijacion      (Scan)
 */

/**
 * Barrier
 * MPI_Barrier(comm): todos los procesos de comm se detienen hasta que todos han invocado esta operacion.
 *
 * Reduce
 * MPI_reduce(sendbuf, recvbuf, count, datatype, op, root, comm) //LA TIENEN QUE INVOCAR TODOS LOS PROCESOS
 *      sendbuf     = direccion de donde sale la informacion
 *      recvbuf     = direccion donde se queda la informacion
 *      count       = cantidad de datos que cada proceso tiene que enviar al raiz
 *      datatype    = tipo de datos
 *      op          = tipo de operacion (MPI_SUM, MPI_MAX,...)
 *      root        = quien recibe todos los datos
 *      comm        = comunicador
 *
 * AllReduce
 * MPI_Allreduce(sendbuf, recvbuf, count, datatype, op, comm) //LA TIENEN QUE INVOCAR TODOS LOS PROCESOS, asi reciben TODOS
 *      sendbuf     = direccion de donde sale la informacion
 *      recvbuf     = direccion donde se queda la informacion
 *      count       = cantidad de datos que cada proceso tiene que enviar al raiz
 *      datatype    = tipo de datos
 *      op          = tipo de operacion (MPI_SUM, MPI_MAX,...)
 *      comm        = comunicador
*/

void ejercicio_4(int argc, char *argv[])
{
    int n, i, suma, myid, nprocesos, suma_local;
    MPI_Init(&argc, &argv);

    if (argc == 2)
        n = atoi(argv[1]);
    else
        n = 10;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocesos);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    suma_local = 0;
    for (i = myid + 1; i <= n; i += nprocesos)
        suma_local += i;

    // HACIENDO UNA REDUCCION TODOS A UNO
    MPI_Reduce(&suma_local, &suma, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    // HACIENDO UNA REDUCCION TODOS A TODOS
    MPI_Allreduce(&suma_local, &suma, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    if (myid == 0)
    {
        printf("Valor de la suma = %d\n", suma);
    }

    MPI_Finalize();
    // return 0;
}

// ejercicio 1-15
void desplazar(double vloc[], int mb)
{
    int myid, nprocesos, next, prev;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocesos);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    // si soy primero
    if (myid == 0)
    {
        prev = nprocesos - 1;
    }
    else
    {
        prev = myid - 1;
    }
    // si soy ultimo
    if (myid == nprocesos - 1)
    {
        next = 0;
    }
    else
    {
        next = myid + 1;
    }

    double elemento_a_enviar = vloc[mb - 1];
    for (int i = mb - 1; i > 0; i--)
    {
        vloc[i] = vloc[i - 1];
    }
    MPI_Sendrecv(&elemento_a_enviar, 1, MPI_DOUBLE, next, 0, vloc, 1, MPI_DOUBLE, prev, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); //o poner &vloc[0]
}