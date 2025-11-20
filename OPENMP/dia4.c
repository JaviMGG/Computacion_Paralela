/**
 * Condiciones de Bernstein:
 * ----------------------------------------------------------------------------------
 * Dos tareas son independientes si:
 * - la segunda tarea se puede poner en el lugar de la primera y NO se produce fallo.
 * · Condiciones:
 *
 *      Input_I ----> I ----> Output_I ---> Input_J ----> J ----> Output_J (ESTE EJEMPLO ES UNA DEPENDENCIA DE FLUJO)
 *
 *      - lo que meto en la tarea I y lo que saco en la tarea J, NO ESTÁN RELACIONADOS (interseccion da valor nulo)
 *      - lo que meto en la tarea J y lo que saco en la tarea I, NO ESTÁN RELACIONADOS (interseccion da valor nulo)
 *      - lo que sale tanto en I y J, NO ESTÁN RELACIONADOS (interseccion da valor nulo)
 *
 */

/**
 * double a = 3, b = 5, c,d;
 * c = T1(a,b);
 * d = T2 (a,b,c);      //la tarea C tiene dependencia (similar a lo de ETC de lectura y escritura)
 */

// si hay dependencias, NO es paralelizable

/**
 * Grafo de dependencias:
 * ----------------------------------------------------------------------------------
 * Representacion grafica/Abstraccion acíclica de las dependencia de los datos (van de arriba a abajo).
 * Sirve para ver las dependencias y como asignar las tareas para evitar dependencias.
 * ----------------------------------------------------------------------------------
 * Cada nodo es una TAREA (puede tener un coste).
 * Cada arista representa la dependencia (puede tener un coste).
 *
 * 2 tareas al mismo nivel (un nodo al lado del otro): NO hay dependencia.
 * 1 tarea por debajo de otra: TIENE dependencia.
 * ----------------------------------------------------------------------------------
 * Longitud del camino: SUMA DE LOS COSTES de los nodos que componen UN CAMINO.
 * Camino crítico: camino más largo DESDE EL NODO INICIAL hasta el nodo final.
 * ----------------------------------------------------------------------------------
 * A nivel de concurrencia:
 *      · Grado máximo = nº tareas MÁXIMO que se ejecutan concurrentemente (concurrentemente = por niveles).
 *      · Grado medio = (formula en las transparencias).
 * ----------------------------------------------------------------------------------
 * Si hay un nodo que tiene coste 1 y otro que tiene coste 5 y se ejecutan CONCURRENTEMENTE, el que tiene coste 1 esperará a que termine el otro para avanzar
 */

/**
 * ejercicio 1-2:
 *      Longitud del camino crítico L = 6
 *      Grado máximo = 2
 *      Grado medio = (suma de los costes de todas las tareas) / camino crítico
 *
 *      Longitud del camino crítico L = 7 (hay que sumar los valores de cada nodo para calcularla)
 *      Grado máximo = 3
 */

/**
 * Repartir segun identificador de hilos
 * int myindex;
 * #pragma omp parallel private (myindex)
 *
 *
 * bucle compartido:
 * #pragma omp parallel{
 *
 * #pragma omp for
 * for(...){...}
 *
 * }
 *
 * CONSTRUCCION SECTIONS:
 * #pragma omp parallel sections{
 *      #pragma omp section
 *      Xaxis();
 *      #pragma omp section
 *      Yaxis();
 *      #pragma omp section
 *      Zaxis();
 * }
 *
 * CONSTRUCCION SINGLE (solo hace una, el resto se las salta)
 * #pragma omp parallel{
 *      #pragma omp single nowait   //solo hace el printf 1, el resto se espera a que termine
 *          printf("holamundo");
 *      work1();
 *
 *      #pragma omp single
 *      printf("holamundo");        //solo hace el printf 1, el resto avanzan
 *
 *      #pragma omp single nowait
 *      printf("holamundo");
 * }
 *
 */

// ejercicios boletines
//----------------------------------------------------------------------------------

#include <omp.h>

int busqueda(int x[], int n, int valor)
{
    int encontrado = 0, i;
#pragma omp parallel private(i)
    {
        int i = omp_get_thread_num();
        int nhilos = omp_get_num_threads();
        while (!encontrado && i < n)
        {
            if (x[i] == valor)
            {
                encontrado = 1;
            }
            i += nhilos; // porque hay que ir de nhilos en nhilos
        }
    }
    return encontrado;
}

// depende de la cantidad de hilos que pongas

#include <stdio.h>
// ejercicio 12
int main()
{
#pragma omp parallel private(a, niteraciones, miID)
    {
    int niteraciones = 0;
    int miID = omp_thread_thread_num();
    #pragma omp for //el private(a) se puede poner a continuacion del for
        for (int i = 0; i < n; i++)
        {
            a = funcion(i * 4);
            v[i] = a + 2;
            niteraciones++;
        }
    }
    printf("soy el hilo %d y he ejecutado %d iteraciones", miID, niteraciones);
}