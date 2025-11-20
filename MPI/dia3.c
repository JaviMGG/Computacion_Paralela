#include <mpi.h>
#include <stdio.h>
/**
 * Esquemas algoritmicos
 * ---------------------
 * Paralelismo de datos:
 *      - memoria compartida   => se paralelizan los bucles ( a veces no es necesario paralelizar debido a que los datos ya estan distribuidos)
 *      - paso de mensajes     => se realiza particionado de datos explicito
 * 
 * Paralelismo de tarea:
 *      - Asignacion estatica no es viable/ tiene problemas de carga desequilibrada.
 *      - Asignacion dinámica: se van asignando tareas a los procesos a medida que quedan ociosos.
 *      - Esquema asimétrico (maestro-esclavo) => un proceso maestro lleva cuenta de las tareas hechas/por hacer. 
 *          Trabajadores reciben las tareas y notifican cuando han terminado.
 *          Depende de la velocidad de los trabajadores, no se puede saber quien trabaja más
 *      - Esquema simétrico (todos iguales) => todos los procesos son iguales y se asignan tareas de forma dinámica. 
*/

/**
 * Operaciones de comunicacion colectiva:
 *      - Sincronizacion   (Barrier)
 *      - Difusion         (Bcast)
 *      - Reparto          (Scatter)
 *      - Recogida         (Gather)
 *      - Multi-recogida   (Allgather)
 *      - Todos a todos    (Alltoall)
 *      - Reduccion        (Reduce)
 *      - Prefijacion      (Scan)
*/

/**
 * Difusion 
 *  MPI_Bcast(buffer, count, datatype, root, comm)
 *      buffer: buffer de datos
 *      count: número de datos
 *      datatype: tipo de datos
 *      root: proceso raíz
 *      comm: comunicador
 *
 * ejemplo sencillos: reparto de examen <=> Uno reparte a TODOS
*/

void ejemploDifusion(){
    double val;
    MPI_Status status;
    int rank, p, i;
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    /*if(rank == 0){
        read_value(&val); //valor que difundiremos
        for(i = 1; i < p; i++){
            MPI_Send(&val, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
        }
    }else{
        MPI_Recv(&val, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
    }*/

    //todo esto sustituye a lo de arriba
    if(rank == 0){
        read_value(&val);
    }
    MPI_Bcast(&val, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);    
    
}

/**
 * Reparto:
 * ----------------------------
 *  Distribuye una serie de fragmentos consecutivos del buffer al resto de procesos (incluyendo a él mismo).
 *  Ejemplo sencillo: reparto de cartas de una baraja: al tener posiciones consecutivas, el proceso 0 tiene la carta 0, el proceso 1 la carta 1, etc.
 * 
 *  MPI_Scatter(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm)
 *      de envio
 *      --------------
 *      sendbuf: buffer de datos que tiene los datos/direccion inicial del vector
 *      sendcount: número de datos que se envian a cada uno de los procesos (todos tendran la misma cantidad de información). NO ES LA TOTALIDAD DEL VECTOR COMPLETO
 *      sendtype: tipo de datos 
 *      
 *      de recepcion
 *      --------------
 *      recvbuf: buffer de datos 
 *      recvcount: número de datos que recibe
 *      recvtype: tipo de datos que recibe
 *      root: proceso raíz (el que envia)
 *      comm: comunicador 
 * 
 * Version asimétrica:
 * ----------------------------
 * MPI_Scatterv(sendbuf, sendcounts, displs, sendtype, recvbuf, recvcount, recvtype, root, comm)
 *      sendbuf: buffer de datos
 *      sendcounts: número de datos
 *      displs: desplazamiento
 *      sendtype: tipo de datos
 *      recvbuf: buffer de datos
 *      recvcount: número de datos
 *      recvtype: tipo de datos
 *      root: proceso raíz
 *      comm: comunicador
 *
*/

//ejemplo de reparto:
void ejemplo_Envio (int argc, char *argv[]){
   int i, myproc;
   int a[15], b[5];
   MPI_Init(&argc, &argv);
   MPI_Comm_rank(MPI_COMM_WORLD, &myproc);
   if(myproc == 0){
      for(i = 0; i < 15; i++){
         a[i] = i;
      }
   }
   MPI_Scatter(a, 5, MPI_INT, b, 5, MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Finalize();

}

/**
 * Recogida:
 * ----------------------------
 * Cada proceso envia un mensaje a root el cual almacena de forma ordenada de acuerdo al indice del proceso en el buffer de recepcion.
 * 
 * MPI_Gather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm)
 *      sendbuf: buffer de datos por parte de cada uno de los procesos (direccion del proceso A, B,...)
 *      sendcount: número de datos que envia cada proceso 
 *      sendtype: tipo de datos (MPI_INT,...)
 *      recvbuf: buffer de datos por parte del proceso raíz (direccion del proceso 0)
 *      recvcount: número de datos que recibimos de cada proceso (debe coincidir con sendcount)
 *      recvtype: tipo de datos 
 *      root: proceso raíz
 *      comm: comunicador
 *
 * si quisieramos que todos recogieran/multi recogida, se usaria MPI_Allgather
 * MPI_Allgather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm)
*/

//ejercicio 2-11
double scalarprod(double X[], double Y[], int n)
{
    double Xlcl[MAX], Ylcl[MAX];
    int nproc, k;
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    k = n/nproc;
    
    //esto ha sido reparto 
    MPI_Scatter(X, k, MPI_DOUBLE, Xlcl, k, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatter(Y, k, MPI_DOUBLE, Ylcl, k, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    double sum_De_CADA_VECTOR = 0.0;
    for(int i = 0; i < k; i++){
        sum_De_CADA_VECTOR += Xlcl[i]*Ylcl[i];
    }

    //esto es la recogida, al juntar cada vector, obtenemos uno distinto a los de cada proceso 
    double global_sum;
    MPI_Reduce(&sum_De_CADA_VECTOR, &global_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    
    //apartado c: modifica el codigo para que todos reciban la suma
    MPI_AllReduce(&sum_De_CADA_VECTOR, &global_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    return global_sum;
}

//ejercicio 2-1
void ejercicio2_1(){
    int i, j, myid, nproc,k;
    int A[N][N], v[N], x[N], Alocal[N][N], xlocal[N];
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    if(myid == 0)
        leer(A,v);

    MPI_Bcast(v,N,MPI_INT, 0, MPI_COMM_WORLD);
    k = N/nproc;
    MPI_Scatter(A, k*N, MPI_INT, Alocal, k*N, MPI_INT, 0, MPI_COMM_WORLD);
    
    for(i = 0; i <k; i++){
        xlocal[i] = 0;
        for(j = 0; j < N j++){ //la N se mantiene igual porque las columnas NO SE MODIFICAN
            x[i] += Alocal[i][j]*v[j];
        }
    }

    MPI_Gather(xlocal, k, MPI_INT, x, k, MPI_INT, 0, MPI_COMM_WORLD);
    //si pidiera que todos recibieran los vectores, se usaria AllGather
    MPI_ALLgather(xlocal, k, MPI_INT, x, k, MPI_INT, MPI_COMM_WORLD);

    if(myid == 0)
        escribir(x);
}


//ejercicio 2-2
void ejercicio2_3(){
    int i, j;
    double s, norm, A[N][N] Alocal[N][N];
    int myid, nproc;
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    if(myid == 0){
        leermat(A);
    }
    k = N/nproc;
    for(int i = 0; i < k; i++){
        MPI_Scatter(A[i*nproc], N, MPI_DOUBLE, Alocal[i], N, MPI_DOUBLE, 0, MPI_COMM_WORLD) //si pongo solo 1 elemento, indica el primer elemento de la fila de la matriz, se puede poner tambien &A[i*nproc][0]
    }

    s = 0.0;
    for(i = 0; i < k; i++){
        for(j = 0; j < N j++){ //la N se mantiene igual porque las columnas NO SE MODIFICAN
            s += Alocal[i][j]*Alocal[i][j];
        }
    }
    double s_local = 0.0;
    MPI_Reduce(&slocal, &s, 1, MPI_DOUBLE,MPI_SUM, 0, MPI_COMM_WORLD);
    if(myid == 0){
        norm = sqrt(s);
        prinf("norm=%f\n", norm);
    }
}