//Autor: JaviMGG
#include <mpi.h>
/**
 * Evaluacion de prestaciones
 * ----------------------------
 * Tiempo de ejecucion Paralelo:
 *      tiempo desde que comienza a ejecutar el PRIMERO hasta que acaba el ULTIMO
 * EJEMPLO: si tenemos 3 procesos que tardan 1,2 y 3 segundos respectivamente, el tiempo de ejecucion paralelo es 3 segundos (el tiempo es el máximo con respecto a los otros dos)
 * 
 * Hay una etapa de solapamiento (overlap) que ocurre cuando un proceso termina y otro empieza
 * Tproceso = ta + tb - Toverlap
 * 
 * Tiempo de establecimiento de la comunicacion (ts):
 *      tiempo que tarda en establecer la comunicacion
 * 
 * Ancho de banda (W):
 *      ancho de banda de la comunicacion
 * 
 * tiempo de envio de 1 byte => tw = 1/W
 * 
 * tiempo necesario para envia un mensaje de n bytes => ts + tw => ts + n/W
 * 
 * speedup      = Tsecuencial / Tparalelo
 * Tsecuencial  = tiempo que tarda en ejecutar el proceso en secuencial
 * Tparalelo    = tiempo que tarda en ejecutar el proceso en paralelo
 * 
 * casos posibles del speedup:
 *      "caso speed-down"   => speedup < 1     => algoritmo paralelo es mas lento que el secuencial
 *      "caso sublineal"    => 1 < speedup < p => algoritmo paralelo es mas rapido que el secuencial PERO no aprovecha la capacidad de los procesos
 *      "caso lineal"       => speedup = p     => algoritmo paralelo es mas rapido que el secuencial Y aprovecha la capacidad de los procesos al 100%
 *      "caso superlineal"  => speedup > p     => situacion anomala, el algortmo paralelo tiene menor coste que el secuencial
 * 
 * eficiencia = speedup / p
 * 
 *------------------------------
 * Ley de Amdahl
 *------------------------------
 * dado un algoritmo secuencial, descomponemos en dos partes:
 *      - parte paralela (p)
 *      - parte secuencial (1-p)
 * 
 * tiempo minimo paralelo alccanzable será: ts + tp/p
 * 
 * speedup minimo = 1 / (ts + tp/p) 
 * 
 * Ejemplo: tn = 100 segundos, 25 es tiempo secuencial, 75 es tiempo paralelo
 * t(n,p) = 25 + 75/p
 * speedup maximo = lim p->infinito 100 / (25 + 75/p)
 * speedup maximo = 100 / (25 + 0)
 * speedup maximo = 4 
 * (siempre valdra 4 en este ejemplo porque la parte secuencial no se puede paralelizar)
 * 
 * -------------------------
 * Coste de comunicacion
 * -------------------------
 * enviar y recibir datos
 * tc = ts + n*tw (siendo n el numero de datos y tw el tiempo de envio de 1 dato de 1 tipo concreto)
 * 
 * --------------------------
 * Broadcast (difusion)
 * --------------------------
 * tc = (p-1)* (ts + n*tw)
 
 * ##########################
 * # IMPORTANTE: Si fuera una red de tipo bus, donde uno envia y el resto escucha, con un envio seria suficiente:
 * #             tc = ts + n*tw
 * ##########################
 * 
 * --------------------------
 * Reduce (reduccion)
 * --------------------------
 * n es NUMERO DE COMPONENTES DEL VECTOR, SI FUERA UNA MATRIZ, N SERIA EL NUMERO DE ELEMENTOS DE LA MATRIZ (FILAS * COLUMNAS)
 * tc = (p-1)* (ts + n*tw) + (p-1)*n
 * (p-1) = numero de procesos que envian al raiz
 * (ts + n*tw) = tiempo de envio de datos
 * (p-1)*n = tiempo de recepcion de datos (p-1) porque el raiz no recibe nada. Es p-1 porque si tengo por ejemplo 3 + 5 + 7, se produce 2 sumas
 * 
 * --------------------------
 * Allreduce (reduccion)
 * --------------------------
 * tc = 2* (p-1)* (ts + n*tw) + (p-1)*n
 * 2* (p-1)* (ts + n*tw) = tiempo de envio de datos (se multiplica por 2 porque envia a todos)
 * (p-1)*n = tiempo de recepcion de datos. (p-1) porque el raiz no recibe nada
 * 
 * ##########################
 * # IMPORTANTE: si se tratara de una red de tipo bus: 
 *              tc = p * (ts + n*tw) + (p-1)*n
 * ##########################
 * 
 * --------------------------
 * Scatter (distribucion)
 * --------------------------
 * Proceso raiz dispone de n datos que reparte entre los p procesos. Se realizan p-1 envios de n/p datos a cada uno de ellos
 * tc = (p-1)* (ts + n/p*tw)
 * 
 * --------------------------
 * Gather (recopilar)
 * --------------------------
 * Cada proceso envia n/p datos al raiz. Se realizan p-1 envios de n/p datos
 * tc = (p-1)* (ts + n/p*tw)

 * --------------------------
 * Allgather (recopilar)
 * --------------------------
 * Cada proceso envia n/p datos a todos los procesos. Se realizan p-1 envios de n/p datos
 * tc = (p-1)* (ts + n/p*tw) + (p-1)*(ts+n*tw)
 * 
 * ##########################
 * # IMPORTANTE: si se tratara de una red de tipo bus: 
 *              tc = (p-1)* (ts + n/p*tw) + (ts+n*tw)
 * ##########################
*/



//ejercicio 2-1:    
/**
    coste secuencial: sumatorio desde 0 hasta n-1 de (sumatorio de 0 a n-1 de 2) => 2n^2
    coste paralelo: coste de scatter + coste de Bcast + coste de Gather + coste de la operacion local
    coste de scatter: (p-1)* (ts + n/p*n*tw)    => p*ts + n^2*tw
    coste de Bcast: (p-1)* (ts + n*tw)          => p*ts + p*n*tw
    coste de Gather: (p-1)* (ts + n/p*tw)       => p*ts + n*tw
    coste de la operacion local: 2n^2/p         => 2n^2/p
*/    
void ejercicio2_1(){
    int i, j, myid, nproc,k;
    int A[N][N], v[N], x[N], Alocal[N][N], xlocal[N];
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    if(myid == 0)
        leer(A,v);

    k = N/nproc;
    MPI_Bcast(v,N,MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatter(A, k*N, MPI_INT, Alocal, k*N, MPI_INT, 0, MPI_COMM_WORLD);
    
    for(i = 0; i <k; i++){
        xlocal[i] = 0;
        for(j = 0; j < N; j++){ //la N se mantiene igual porque las columnas NO SE MODIFICAN
            xlocal[i] += Alocal[i][j]*v[j];
        }
    }

    MPI_Gather(xlocal, k, MPI_INT, x, k, MPI_INT, 0, MPI_COMM_WORLD);
    //si pidiera que todos recibieran los vectores, se usaria AllGather
    //MPI_ALLgather(xlocal, k, MPI_INT, x, k, MPI_INT, MPI_COMM_WORLD);

    if(myid == 0)
        escribir(x);
}



//ejercicio 2-11
/**
    apartado b)
    coste secuencial: sumatorio desde 0 hasta n-1 de 2 => 2n
    coste paralelo: 2 * coste de scatter (PORQUE HAY 2 SCATTER) + coste de reduce + coste de la operacion local
    coste de scatter: 2 * ((p-1)* (ts + n/p*n*tw))     => 2p*ts + 2n*tw
    coste de reduce: (p-1)* (ts + tw) + (p-1)          => p*ts + p*tw + p
    coste de la operacion local: 2n/p                  => 2n/p
*/
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
    //MPI_AllReduce(&sum_De_CADA_VECTOR, &global_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    return global_sum;
}



//ejercicio 2-17
int ejercicio2_17(int argc, char *argv[]){
    int i, j, myid, size,k;
    double a, b, A[M][N], B[M][N], C[M][N], Alocal[M][N], Blocal[M][N], Clocal[M][N];
    MPI_Init(&argc, &argv); //SI ESTAMOS EN UN MAIN, SE PONE ESTA LINEA
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    
    if(myid == 0){
        LeeOperandos(A,B, &a, &b);
    }
    
    MPI_Bcast(&a, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD); //para enviar un valor a todos los procesos
    MPI_Bcast(&b, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD); //para enviar un valor a todos los procesos

    k = M/size;
    for(i = 0; i < k; i){
        //el indice de envio es 0,3,6,9,.... el indice de recibo es 0,1,2,3,...
        MPI_Scatter(A[i*size], N, MPI_DOUBLE, Alocal[i], N, MPI_DOUBLE, 0, MPI_COMM_WORLD); //cada vez enviamos solo 1 fila, entonces no es k*N sino N
        MPI_Scatter(B[i*size], N, MPI_DOUBLE, Blocal[i], N, MPI_DOUBLE, 0, MPI_COMM_WORLD); //cada vez enviamos solo 1 fila, entonces no es k*N sino N
    }
    
    for(i = 0; i < k; i++){
        for(j = 0; j < N; j++){
            Clocal[i][j] = a*Alocal[i][j] + b*Blocal[i][j];
        }

    }

    for(i = 0; i < k; i++){
        MPI_Gather(Clocal[i], N, MPI_DOUBLE, C[i*size], N, MPI_DOUBLE, 0, MPI_COMM_WORLD); 
    }
    
    
    if(myid == 0)
        EscribirOperandos(C);
    
    MPI_Finalize();
    return 0;
}