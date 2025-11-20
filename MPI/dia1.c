/**
 * En MPI se abren todos los procesos desde el minuto cero
 * 
 * Procesos organizados en grupos
 * Comunicador = grupo + contexto 
 * comunicacion en un comunicador NO PUEDE INFERIR EN LA DE OTRO (util para aislar la comunicacion dentro de un libreria)
 * 
 * Tipos de comunicaciones:
 * - Punto a punto:  1 proceso envia, 1 proceso recibe
 *      1) send(x, 1) => envio el contenido x al comunicador 1
 *      2) recv(y, 0) => recibo el contenido y del comunicador 0 (es un id)
 * - Colectivas: +2 procesos envian, +2 procesos reciben
 *      1 proceso raiz => proceso que envia a todos o recibe de todos (1 unico proceso que tiene una funcion distinta al resto)
 * - Operaciones avanzadas (MPI-2, MPI-3): 
 * 
 * 
 * 
 *  * Envio con sincronizacion: si el emisor quiere enviar un mensaje antes que el receptor ejecute la "comunicacion", se espera
 *    ejemplo bastante bueno: dos personas quedan en un sitio concreto para intercambiar informacion, si una llega antes que la otra, se espera. Cuando ambos estan, se intercambian la informacion y se van
 * 
 *  * Envio con buffer (sincrono): buffer almacena una copia temporal del mensaje
 *    ejemplo: una persona espera a que la carta que le tiene que llegar esté en el buzon 
 * 
 *  ** Operacion bloqueante: se espera a que el receptor active el recv o el buffer esté lleno para poder recoger los datos
 *     ejemplo: en lugar de ir hasta el buzon del vecino para darle la carta, será mas eficiente poner la carta en el buzon de correos y que alguien (mientras nosotros seguimos con nuestro codigo) recoja la carta y la entregue 
 * 
 *  ** Operacion no bloqueante: 
 * 
 * 
 * recv(z, any_src, any_tag, status)
 *    - z = que recibo
 *    - any_src = de quien recibo
 *    - any_tag = como ha sido etiquetado
 * 
*/


/**
 * Programar con MPI
 *  
 * www.mpi-forum.org
 * 
 * comunicador = grupo + contexto
 * 
 * - ejecutar codigo paralelo: mpiexec -n p programa [argumentos]
 *  * -n 
 *  * p = nº copias  
 * 
 * Enviar: 
 *      MPI_Send(buf, count, datatype, dest, tag, comm)
 *      buf         = buffer (direccion de memoria de la variable donde queremos enviar los datos, indicamos solo posicion inicial)
 *      count       = nº datos que enviamos / longitud del vector que enviamos
 *      datatype    = tipo de dato que enviamos (entero, float, char...)
 *      dest        = proceso destino
 *      tag         = etiqueta para los mensajes para poder diferenciarlos (tipo entero)
 *      comm        = comunicador donde se situa el envio
 * 
 *  
 * Recibir:
 *      MPI_Recv(buf, count, datatype, src, tag, comm, stat)
 *      buf         = Buffer (direccion de memoria de la variable donde queremos recibir los datos, indicamos solo posicion inicial)
 *      count       = Nº datos que recibimos / longitud del vector que recibimos
 *      datatype    = Tipo de dato que recibimos (entero, float, char...)
 *      src         = Proceso origen
 *      tag         = Etiqueta para los mensajes para poder diferenciarlos (tipo entero). debe coincidir con la que envia
 *      comm        = Comunicador donde se situa el envio
 *      stat        = Si queremos saber de muchos emisores quien ha enviado el mensaje, dice quien lo ha enviado y con que etiqueta. 
 * 
 *      MPI_Irecv(buf, count, datatype, src, tag, comm, stat)
 * 
 * Enviar y Recibir:
 *      MPI_Sendrecv(sendbuf, sendcount, sendtype, dest, sendtag, recvbuf, recvcount, recvtype, source, recvtag, comm, status)
 *      MPI_Sendrecv_replace(buf, count, datatype, dest, sendtag, source, recvtag, comm, status) //tenemos un único buffer.
 *      
 * 
*/

#include <stdio.h>

int main2(int argc, char* argv[]){
    int k; /*rango del proceso*/
    int p; /*numero de procesos*/
    MPI_Init(&argc, &argv);
    //MPI_Comm_rank(MPI_COMM_WORLD, &k); //"soy el proceso k"
    //MPI_Comm_size(MPI_COMM_WORLD, &p); //numero de procesos totales
    printf("Soy el proceso %d de %d\n", k, p);
    MPI_Finalize();
    return 0;
}


//ejercicio 1
void ejercicio1(int argc, char **argv){
    double numero;
    int myid;
    MPI_Status = status;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    if(myid == 0){
        printf("dame un numero: ");
        scanf("%lf",&numero);
        printf("el valor es %f\n", numero);
        MPI_Send(&numero, 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
    }
    else if (myid ==1){
        //MPI_Recv(&numero, 1, MPI_DOUBLE, 0,0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&numero, 1, MPI_DOUBLE,MPI_ANY_SOURCE ,MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        
        
        printf("el valor es %f\n", numero);
    }
    MPI_Finalize();
}

//ejercicio 2
#define N 5
void ejercicio2(int argc, char **argv){
    int a[N], i,myid;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    if(myid == 0){
        for (int i = 0; i < N; i++)
        {
            a[i] = i;
        }
        MPI_Send(a, N, MPI_INT, 1, 0, MPI_COMM_WORLD); //le envias la direccion de memoria del primero elemento del vector //tambien se puede poner en lugar de a => &a[0]
    }
    else if(myid == 1){
        //MPI_Recv(&a[0],N, MPI_INT, 0,0, MPI_COMM_WORLD); 
        MPI_Recv(&a[0],N, MPI_INT, MPI_ANY_SOURCE, MPY_ANY_TAG, MPI_COMM_WORLD, &status); 
        
        printf("los elementos del vector son; \n");
        for (int i = 0; i < N; i++)
        {
            printf("%d\n", a[i]);
        }
    }
    MPI_Finalize();
    
    
}