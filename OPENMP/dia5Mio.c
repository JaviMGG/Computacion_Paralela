/**
 * Comunicacion
 * ---------------------------------------------------------
 * Enviando mensajes: (envia y recibe mensajes explicitos) => Envia(P, dato) Recibe(P0, dato)
 * Ventajas: universalidad, se comprende más facil, es más expresivo y eficiente
 * Desventajas: programacion compleja, control total de las comunicaciones
 * ---------------------------------------------------------
 * MPI: especificacion propuesta por un comité.
 * ---------------------------------------------------------
 * Caracteristicas:
 * Portable a cualquier plataforma paralela
 * Simple (con 6 funciones puedes hacer cualquier programa)
 * Potente (mas de 300 funciones)
 * Estandar principal para Fortran y C
 * ---------------------------------------------------------
 * Basada en funciones de biblioteca => #include <mpi.h>
 * OBLIGATORIO:
 *              Llamar a MPI_Init();
 *              MPI_Finalize();
 * ---------------------------------------------------------
 * Tipos de operaciones (la mayoria opera sobre comunicadores):
 * Punto a Punto                => Intercambio de informacion entre 2 procesos
 * Comunicacion Colectiva       => Intercambio de informacion entre conjuntos de procesos
 * Gestión de datos             => Tipos de datos derivados
 * Operaciones de Alto Nivel    => Grupos, comunicadores, atributos...
 * Operaciones Avanzadas        => E/S, Creacion de procesos...
 * Utilidades                   => Interaccion con el entorno
 * ---------------------------------------------------------
 * COMUNICADORES: Abstraccion que engloba los siguientes conceptos:
 * GRUPO: Conjunto de procesos
 * CONTEXTO: Evita interferencias entre mensajes distintos
 * ---------------------------------------------------------
 * Comunicador agrupa P procesos
 * int MPI_Comm_size(MPI_Comm comm, int *size)
 * Cada proceso tiene un rango(identificador) desde 0 a P-1
 * int MPI_Comm_rank(MPI_Comm comm, int *rank)
 * 
*/

// #include <mpi.h> //necesita instalar en el portatil MPI para que no de error
// void ejemplo_funcion_MPI(){
//     int k; /* rango del proceso */
//     int p; /* número de procesos */
//     MPI_Init(&argc, &argv); --------------------> IMPORTANTE PONERLO
//     MPI_Comm_rank(MPI_COMM_WORLD, &k); ---------> CADA PROCESO CON SU ID
//     MPI_Comm_size(MPI_COMM_WORLD, &p); ---------> AGRUPACION DE PROCESOS
//     printf("Soy el proceso %d de %d\n", k, p);
//     MPI_Finalize(); ----------------------------> IMPORTANTE PONERLO
//     return 0;


/**
 * Comunicacion Punto a punto:
 * ---------------------------------------------------------
 * Mensaje enviado explicitamente por el emisor y recibidos por el receptor (logico)
 * ---------------------------------------------------------
 * ENVIO estándar:
 * MPI_Send(buf, count, datatype, dest, tag, comm)
 * 
 * RECEPCION estándar:
 * MPI_Recv(buf, count, datatype, src, tag, comm, stat)
 * ---------------------------------------------------------
 * LOS 3 PRIMEROS PARAMETROS - EL MENSAJE:
 * buf          => buffer de memoria donde se almacena la informacion
 * count        => nº elementos que componen el mensaje
 * datatype     => tipo de datos (MPI_INT, MPI_DOUBLE...)
 * ---------------------------------------------------------
 * LOS 5 SIGUIENTES PARAMETROS - EL SOBRE:
 * dest         => destino (indicado mediante ID de procesos)
 * src          => origen de la comunicacion (indicado mediante ID de procesos)
 * comm         => comunicación permitida solo dentro del mismo comunicador
 * tag          => numero entero para distinguir mensajes de tipos distintos
 * stat         => estado que contiene informacion (proceso emisor, etiqueta, longitud del mensaje)
 * (req)        => para reutilizar el buffer, asegurarse de que se ha completado el envio
 * ---------------------------------------------------------
 * Modos de Envío Punto a Punto:
 * SINCRONO:        MPI_Ssend(buf, count, datatype, dest, tag, comm)
 * BUFFER:          MPI_Bsend(buf, count, datatype, dest, tag, comm)
 * ESTANDAR:        MPI_Send (buf, count, datatype, dest, tag, comm)
 * NO BLOQUEANTE:   MPI_Isend(buf, count, datatype, dest, tag, comm, req)
 * 
 * Modos de Recepción Punto a Punto:
 * ESTANDAR:        MPI_Recv (buf, count, datatype, src, tag, comm, stat)
 * NO BLOQUEANTE:   MPI_Irecv(buf, count, type, src, tag, comm, req)
 * ---------------------------------------------------------
 * Operaciones Combinadas:
 * MPI_Sendrecv(sendbuf, sendcount, sendtype, dest, sendtag, recvbuf, recvcount, recvtype, source, recvtag, comm, status)
 * MPI_Sendrecv_replace(buf, count, datatype, dest, sendtag, source, recvtag, comm, status)
 * 
*/

/**
 * Comunicacion Colectiva
 * ---------------------------------------------------------
 * Involucran todos los procesos de un grupo(informados por un comunicador)
 * ---------------------------------------------------------
 * Operaciones disponibles:
 * BARRIER:     MPI_Bcast(buffer, count, datatype, root, comm)
 * BCAST:       MPI_Bcast(buffer, count, datatype, root, comm)
 * SCATTER:     MPI_Scatter(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm)
 * GATHER:      MPI_Gather (sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm)
 * ALLGATHER:   MPI_Allgather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm)
 * ALLTOALL:    MPI_Alltoall (sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm)
 * REDUCE:      MPI_Reduce(sendbuf, recvbuf, count, datatype, op, root, comm)
                MPI_Allreduce(sendbuf, recvbuf, count, type, op, comm)
 * SCAN:        MPI_Scan(sendbuf, recvbuf, count, datatype, op, comm)
*/

/**
 * Vectores:
 * MPI_Type_vector(count, length, stride, type newtype)
 * count    => cuantos bloques compone
 * length   => longitud de los bloques
 * stride   => separacion entre elemento de un bloque y mismo elemento
 * type     => tipo de elementos individuales
 * newtype  => (modificar el tipo supongo)
 */