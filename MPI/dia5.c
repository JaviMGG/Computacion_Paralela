#include <mpi.h>

/**
 * Esquemas de asignacion estatica
 * --------------------------------
 * 
 * Distribucion de matrices por bloques
 * ------------------------------------
 * asignacion hace corresponder porciones contiguas (bloques) del dominio de datos (matriz)
 distribuciones más usuales:
    - Unidimensionales por bloques:
         
    - Bidimensionales por bloques
        - 
*/

/**
 *  en reparto por bloques:
    paso del indice global a local => 
        I_global / k = I_proceso <=> I_global es de todo el bloque original. I_proceso es el indicie del proceso que le corresponde. k es numero de datos de cada bloque mas pequeño
        I_global % k = I_local   <=> I_global es de todo el bloque original. I_local es el indicie del bloque (más pequeño) que le corresponde. k es numero de datos de cada bloque mas pequeño
    en reparto ciclico:
    paso del indice global a local => I_global = p*I_local + I_proceso
        I_global % p = I_proceso <=> I_global es de todo el bloque original. I_proceso es el indicie del proceso que le corresponde. p es numero de procesos
        I_global / p = I_local   <=> I_global es de todo el bloque original. I_local es el indicie del bloque (más pequeño) que le corresponde. p es numero de procesos
*/
//ejercicio 1-12
void ejercicio1_12(double v[], int n){
    double max = v[0];
    int i, posmax = 0;
    for(i = 1; i < n; i++){
        if(v[i] > max){
            max = v[i];
            posmax = i;
        }
    }    
    printf("Maximo = %f, posmax = %d\n", max, posmax);
}

void ejercicio1_12_paralelo(double v[], int n, int rank, int nproc){
    int k, p, i; 
    double maxlocal = v[0];
    int posmaxlocal = 0;
    k = n/nproc; //numero de elementos por proceso
    if(rank==0){
        for(p = 1; p < nproc; p++){
            MPI_Send(&v[k*p], k, MPI_DOUBLE, p, 0, MPI_COMM_WORLD); //si envio desde el 5 elemento, el proceso 1 recibe en elemento 0
        }
    }
    else{
        MPI_Recv(v, k, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    for(i = 1; i < k; i++){
        if(v[i] > maxlocal){
            maxlocal = v[i];
            posmaxlocal = i;
        }
    }    
    if(rank==0){
        max = maxlocal; //para evitar que se machaque
        posmax = posmaxlocal; //para evitar que se machaque
        for(p = 1; p < nproc; p++){
            MPI_Recv(&maxlocal, 1, MPI_DOUBLE, p, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&posmaxlocal, 1, MPI_INT, p, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            if(maxlocal > max){
                max = maxlocal;
                posmax = posmaxlocal + k*p; //para transformar la posicion local en global
            }
        }   
    }
    else{
        MPI_Send(&maxlocal, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        MPI_Send(&posmaxlocal, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }

    if(rank==0){
        printf("Maximo = %f, posmax = %d\n", max, posmax);    
    }

}

/**
    * Otras funcionalidades:
    ----------------
    TIPOS DE DATOS
    ----------------
    * -MPI_CHAR
    * -MPI_SHORT
    * -MPI_INT
    * -MPI_LONG
    * -MPI_UNSIGNED_CHAR
    * -MPI_UNSIGNED_SHORT
    * -MPI_UNSIGNED
    * -MPI_UNSIGNED_LONG
    * -MPI_FLOAT
    * -MPI_DOUBLE
    * -MPI_LONG_DOUBLE

    ----------------
    Datos multiples
    ----------------
    MPI_Get_count(MPI_Status *status, MPI_Datatype datatype, int *count)
    status: estado de la comunicacion
    datatype: tipo de datos
    count: nº elementos contiguos en memoria, la funcion nos dice cuantos elementos se han recibido (lo almacena en count)

    ----------------------------------------
    Envio de una matriz que no tengan mismo numero de filas que de columnas:
    ----------------------------------------

    //IMPORTANTE Tipos de datos derivados regulares ENTRA SIEMPRE (es de los mas importantes)
    
    MPI_Type_vector(int count, int blocklength, int stride, MPI_Datatype oldtype, MPI_Datatype *newtype)
    count: numero de bloques que vas a enviar (bloque = elementos consecutivos en memoria)
    blocklength: numero de elementos por bloque (todos los bloques miden lo mismo)
    stride: numero de elementos entre bloques (distancia entre bloques: si pusieramos todos los elementos en una sola fila, distancia desde uno hasta otro: [1,2,3,4,5,6,1,2,3,4,5,6] si quisiera saber el stride del 1 al otro 1, seria 6)
    oldtype: tipo de datos
    newtype: nuevo tipo de datos (cambiar el tipo de datos)
    
*/
void codigo_ejemplo_datos_derivados(){
    double A[4][7];
    MPI_Datatype columna;
    MPI_Type_vector(4, 7, 10, MPI_DOUBLE, &columna);   
    MPI_Type_commit(&columna); //esencial, SIRVE PARA DEFINIR EL TIPO DE DATO NUEVO QUE SE USARÁ
    if(my_rank == 0){
        MPI_Send(&A[0][2], 1, columna, 1, 0, comm);   
    }
    else{
        //si quiero que reciba 
        MPI_Recv(&A[0][2], 1, columna, 0, 0, comm, &status);   
    }
    //MPI_Type_free(&columna); //esencial, SIRVE PARA LIBERAR EL TIPO DE DATO NUEVO QUE SE HA USADO
}

//ejercicio 3-08
void ejercicio3_08(int m, int n, double A[M][N], double v[MAX], MPI_Comm comm){
    MPI_Datatype submatriz;
    int myid; 
    MPI_Comm_rank(comm, &myid);
    
    if(myid == 0){
        MPI_Type_vector(m, n, N, MPI_DOUBLE, &submatriz);
        MPI_Type_commit(&submatriz);
        MPI_Send(&A[0][0], 1, submatriz, 1, 0, comm); //A / A[0] / &A[0][0] || ENVIO 1 SUBMATRIZ DE 6 elementos
        MPI_Type_free(&submatriz);
    }else if (myid == 1){
        MPI_Recv(v, m*n, MPI_DOUBLE, 0, 0, comm, MPI_STATUS_IGNORE); //SI FUERA MISMO TAMAÑO DE DATOS, DE PODRIA PONER 1
    }
}

//ejercicio 3-11
void ejercicio3_11(double A[N][N], double B[N/2][N]){
    MPI_Datatype submatriz;
    int myid;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    if(myid == 0){
        MPI_Type_vector(N/2, N, N*2, MPI_DOUBLE, &submatriz);
        MPI_Type_commit(&submatriz);
        MPI_Send(&A[0][0], 1, submatriz, 1, 0, MPI_COMM_WORLD);
        MPI_Send(&A[1][0], 1, submatriz, 2, 0, MPI_COMM_WORLD);
        MPI_Type_free(&submatriz);
    }
    else if (myid == 1 || myid == 2){
        MPI_Recv(&B[0][0], N*N/2, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
}

/**
 * Apartado b: costes
 * tc = 2* (ts+tw*(N*N/2))
*/

//ejercicio 2-17
int ejercicio2_17(int argc, char *argv[]){
    int i, j, myid, size,k;
    double a, b, A[M][N], B[M][N], C[M][N], Alocal[M][N], Blocal[M][N], Clocal[M][N];
    MPI_Init(&argc, &argv); //SI ESTAMOS EN UN MAIN, SE PONE ESTA LINEA
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    
    MPI_Datatype submatriz_envio;

    if(myid == 0){
        LeeOperandos(A,B, &a, &b);
    }
    
    MPI_Bcast(&a, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD); //para enviar un valor a todos los procesos
    MPI_Bcast(&b, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD); //para enviar un valor a todos los procesos

    k = M/size;

    MPI_Type_vector(k, N, size*N, MPI_DOUBLE, &submatriz_envio);
    MPI_Type_commit(&submatriz_envio);

    if(myid == 0){
        int p;
        //esto para A
        MPI_Sendrecv(A[0], 1, submatriz_envio, 0, 0, Alocal, k*N, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); //esto para enviar y recibir a mi mismo
        for(p=1; p < size; p++){
            MPI_Send(A[p], 1, submatriz_envio, p, 0, MPI_COMM_WORLD);
        }

        //esto para B
        MPI_Sendrecv(B[0], 1, submatriz_envio, 0, 0, Blocal, k*N, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); //esto para enviar y recibir a mi mismo
        for(p=1; p < size; p++){
            MPI_Send(B[p], 1, submatriz_envio, p, 0, MPI_COMM_WORLD);
        }
    }
    else{
        MPI_Recv(&Alocal[0][0], k*N, MPI_DOUBLE, 0 ,0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    
    for(i = 0; i < k; i++){
        for(j = 0; j < N; j++){
            Clocal[i][j] = a*Alocal[i][j] + b*Blocal[i][j];
        }

    }

    if(myid == 0){
        int p;
        MPI_Sendrecv(Clocal, k*N, MPI_DOUBLE, 0, 0, C[0], 1, submatriz_envio, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        for(p=1; p < size; p++){
            MPI_Recv(C[p], 1, submatriz_envio, p, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }
    else{
        MPI_Send(Clocal, k*N, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }

    
    if(myid == 0)
        EscribirOperandos(C);
    MPI_Type_free(&submatriz_envio);
    MPI_Finalize();
    return 0;
}