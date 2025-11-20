void sumatorio();
double producto_escalar(double x[], double y[], int n);
double factorial(int n);

int main()
{
    int a = 5; // valor a del padre
#pragma omp parallel private(a)
    {
        a = 8;                                                              // como hemos declarado una variable a, utilizamos ESTA, NO LA DEL PADRE
        printf("Soy hilo %d. Valor de a = %d\n", omp_get_thread_num(), a);  // soy hilo 0, valor de a = 8, //soy hilo 1, valor de a = 8
    }
    printf("Soy hilo %d. Valor de a = %d\n", omp_get_thread_num(), a);      // soy hilo 0, valor de a = 5
#pragma omp parallel private(a)
    {
        // como no hay inicializacion de la variable a, tiene valor indefinido
        printf("Soy hilo %d. Valor de a = %d\n", omp_get_thread_num(), a);  // soy hilo 0, valor de a = ¿?,//soy hilo 1, valor de a = ¿?
    }
#pragma omp parallel // como no hemos puesto private(a), utiliza el valor de a del PADRE
    {
        printf("Soy hilo %d. Valor de a = %d\n", omp_get_thread_num(), a);  // soy hilo 1, valor de a = 5 //soy hilo 0, valor de a = 5
    }

    // el orden de los hilos padre e hijo no tiene por que ser 0 y luego 1, puede que llegue primero 1 y luego 0
}

void sumatorio()
{
    int suma = 0;
/**
 * #pragma omp parallel
 * si ponemos private, DARA MAL porque la variable no existe dentro del pragma y no se podria modificar
 * si ponemos shared/nada, cada hilo utiliza un hilo pero no la comparten (no hacen la suma correctamente)
 * -----------------------------------------------------------------------------------------------------------
 * reduccion: operacion colectiva entre diferentes hilos (+,-,*,/,...)
 * Reduccion OCURRE si TODOS los hilos colaboran al valor de una variable, por ejemplo: suma+=suma, f*=i...
 *
 */

// se crea una variable privada dentro, se mantiene el valor de la variable inicial y cada hilo hace la suma y la acumula
#pragma omp parallel for reduction(+ : suma)
}

double producto_escalar(double x[], double y[], int n)
{
    int i;
    double suma;
    suma = 0;
#pragma omp parallel for reduction(+ : suma) // comparte suma sin de manera que actua como si fuera parallel y como si fuera private
    for (int i = 0; i < n; i++)
    {
        suma = suma + x[i] * y[i];
    }
    return suma;
}

double factorial(int n)
{
    int i;
    double f;
    f = 1;
#pragma omp parallel for reduction(* : f) // f tiene condicion de carrera si no ponemos nada
    for (i = 2; i <= n; i++)
    {
        f *= i;
    }

    return f;
}

/**
 * firstprivate:
 * lastprivate: es como si ejecutaras el hilo de manera secuencial
 *
 * clausula if
 * #pragma omp parallel for if(expresion)
 * si la expresion es falsa: bucle se ejecuta de manera secuencial
 * si la expresion es cierta: bucle se ejecuta en paralelo con los hilos que hayamos creado/decidido
 *
 * diferencia entre poner el pragma antes de 2 bucles for y en medio de ambos:
 *
 * si lo ponemos antes de ambos bucles, reduce tiempo de ejecucion (cuesta menos)
 * si lo ponemos entre ambos bucles, va a crear varias regiones paralelas (activa, desactiva, activa, desactiva...), el de fuera será secuencial y el de dentro paralelo
 *
 */
// ANOTACION: todas las variables de control de los bucles son privadas (la j)

int M = 2, N = 3;
double funcion(double A[M][N], double x[N], double y[M])
{
    int i, j;
    double suma, suma2;
    suma2 = 0;
    //    #pragma omp parallel for private(suma,j) reduction(+:suma2)  //reparte los elementos del vector y para cada hilo, //para cada hilo, suma tiene que valer 0
    for (i = 0; i < M; i++)
    {
        suma = 0;
#pragma omp parallel for reduction(+ : suma)
        for (j = 0; j < N; j++)
        {
            suma = suma + A[i][j] * x[j];
        }
        y[i] = suma;
        suma2 = suma2 + y[i];
    }
    return suma2;
}

/**
 *
 */

// en los examenes hay que hacerlo eficientemente INDEPENDIENTEMENTE DE SI LO PONE O NO
void ejercicio_Diez(double A[M][N])
{
    int i, j;
    double prod;
#pragma omp parallel for private(i) // es mejor que el recorrido se de izquierda a derecha que de arriba a abajo
    for (j = 0; j < N; j++)
    {
        for (i = 1; i < M; i++)
        {
            A[i][j] = 2.0 * A[i][j];
        }
    }
    prod = 1.0;

#pragma omp parallel for private(j) reduction(*:prod)
    for (i = 0; i < M; i++)
    {
        for (j = 0; j < N; j++)
        {
            prod = prod * A[i][j];
        }
    }
    return prod;
}


/**
 * planificar de manera estatica: cada hilo sabe quien le va a tocar 
 * planificar de manera dinamica: cada hilo se adapta a la ejecucion actual 
 * 
*/