#include <stdio.h>
#include <omp.h>

void parallelo(){
    #pragma omp parallel for [clausula[]]
    //SOLO SE PUEDEN PARA BUCLES FOR DONDE SE SEPAN EL NUMERO DE ITERACIONES QUE SE HACEN
    for (int i = 0; i < 10; i++){/*...*/}
    
}

void combinacionlineal(double x[], double y[], double z[], double a, double b, int n){
    int i;
    #pragma omp parallel for
    for (i = 0; i < n; i++)     //por defecto es privada
        z[i] = a*x[i] + b*y[i]; //NO HAY condicion de carrera porque Z es un VECTOR y cada hilo accede a una posicion DISTINTA
    
    
}

void combinacionlineal2(double x[], double y[], double z[], double a, double b, int n){
    int i;
    double temp;
    #pragma omp parallel for private(temp)
    for (i = 0; i < n; i++){     //por defecto es privada
        //double temp;             si se declarara aqui, seria tambien privada         
        temp = a*x[i] + b*y[i];  //NO HAY condicion de carrera porque Z es un VECTOR y cada hilo accede a una posicion DISTINTA
        z[i] = temp;
    }
}


int main(){
    #pragma omp parallel
    /**alternativa:
     * #pragma omp parallel{
           printf("soy el hilo %d\n", omp_get_thread_num());
     * }
     * 
     */
    printf("soy el hilo %d\n", omp_get_thread_num());
    
    /**
    #pragma omp parallel private (omp_get_num_threads()){
        if(omp_get_thread_num() == 0){
            printf("somos %d hilos\n", omp_get_num_threads());
        }
    }
    */ 

    return 0;
}