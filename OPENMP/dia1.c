#include <stdio.h>  //LIBRERIA USADA PARA ENTRADA-SALIDA: printf, scanf, fopen...

//& = direccion de memoria de una variable, NO EL CONTENIDO
//* = accede al dato que hay detras de la direccion de memoria

void primerEjemplo(){
  double a[4] = {1.1,2.2,3.3,4.4};
  double *p,x;
  p = &a[2]; 
  x = *p;    
  *p = 0.0;
  p = a;  //es lo mismo que poner &a[0]
}

void segundoEjemplo(){
  char s[] = "Comput. Paralela";
  
  char *p = s;
  while (*p != "P")
  {
    p++; //si fueran float y entero, se incrementa en 4, si es double en 8
  }
}

void tercerEjemplo(){
  void *p;          //no indicas el tipo de dato que se va a guardar
  double x = 10.0, z;
  p = &x;
  z = *(double*)p;  //hay que hacer una conversion porque p es generico y no sabe cuantos bits necesita para almacenar (no sabe si es entero o double)
}

double valorSuma(double x){  //como es double, tiene que devolver un double, SOLAMENTE UN DATO
  return x+1.44;
}

void mensaje(int x){         //como es void y no devuelve nada, no hace falta return, PUEDO DEVOLVER MAS DE UN OBJETO
  printf("valor %d", x);
}

float funcionResta(float a, float b){
  return a-b;
}

float funcionRestaPorDiferencia(float *a, float *b){
  return *a - *b; //si se hace un paso por referencia, SIEMPRE QUE APAREZCA ALGUNO DE LOS PARAMETROS, HAY QUE PONER *
}

#include <math.h>
void calculaRacies(double v[], int m){
  int i;
  for (i=0; i< m; i++){
    v[i] = sqrt(i+1);
  }
}

#include <stdlib.h>
#include <stdio.h>



void ejemplo(int n, double a, double *x, double *y, double *z){
  int i;
  #pragma omp parallel for  //significa que lo hará el padre, for: las iteraciones del bucle for se reparten entre los hilos
  for (i = 0; i < n; i++)   //i se comparte, para que todos hagan lo mismo, tiene que ser private
    z[i] = a*x[i]+y[i];
  //cuando los hilos llegan a una parte secuencial, se quedan suspendidos

  /**
   * siempre empiezan asi: #pragma omp <directiva> [clausula[]]
   * en el compilador hay que poner gcc -fopenmp (si no ponemos eso, hara el codigo secuencial/ignorará los pragmas)
   * si ponemos num_threads(numero_de_hilos), se haran un numero concreto de hilos /el padre + los hijos | tambien se puede poner omp_set_num_threads() ANTES de region paralela
   * en terminal se puede poner OMP_NUM_THREADS = 4 ./ejecutable.c => se pondrá en marcha, cuando llegue a region paralela, se ejecutara en este caso con 4 hilos.
   * con get_num_threads() => si es zona secuencial, 1, si es paralela, X hilos
  */

 /**
  * #pragma omp parallel [clause[clause ...]]{
  * //codigo
  * }
  */

}

void main(){
  /**
   * dia 1:
   * float w = funcionResta(6.0, 3.5);
   * float x = 3.9, y = 1.4;
   * float z = funcionRestaPorDiferencia(&x, &y);
  */

  /**
   * dia 2
    int i, n;
    double *x;  //es puntero
    printf("dime el numero de elementos");
    scanf("%d", &n);                          //pasamos la direccion de memoria de n
    x = (double *) malloc(n*sizeof(double));  //malloc reserva un espacio de valor 8*sizeof(). x será un puntero a double
    calculaRacies(x,n);
    for ( i = 0; i < n; i++){
      printf("sqrt(%d) = %f\n", i+1, x[i]);
    }
    free(x);                                  //accion contraria a un malloc
    //cuando pasemos una matriz estatica, hay que pasarle la direccion de memoria del primer elemento
    //cuando pasemos una matriz dinamica(no se sabe el tamaño/el usuario lo dice en tiempo de ejecucion), hay que reservar el espacio con malloc

   * reservar espacio para vector de punteros a double:
   * A = (double **)malloc(m*sizeof(double *))  EJEMPLO: A vale #300 que apunta a #500 que es un puntero a un double (es puntero a un puntero de un double)
   * lo que hay en el sizeof, indica que lo que se reservan son espacios para punteros double.
   * 
  
  */

}