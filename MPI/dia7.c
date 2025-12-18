#include <mpi.h>
#include <stdio.h>

//ejercicio 1.1
/**
 * coste secuencial:
 t1 = 2n^2  
 t2 = 2n^2
 t3 = 2n^2
 t4 = 3n^2
 t5 = 3n^2
 t6 = 2n + 1
 * 3*2n^2 + 2*3n^2 + 2n + 1 => 12n^2 +2n +1 ==> 12n^2
 *----------------------------------
 * coste paralelo:
 * taritmetico = 2n^2 + 3n^2 + 2n => 5n^2 + 2n ==> 5n^2
 * tcomunicaciones = 2(ts + tw) + 2(ts + n*tw) + 2(ts + n*tw) + (ts + tw) + 2(ts + tw) => 9ts + (5 + 4n)t ==> 9ts + 4nt (el 5 se elimina por ser menos relevante)
 * ttotal = 5n^2 + 9ts + 4nt
 *----------------------------------
 * speedup = 12n^2 / (5n^2 + 9ts + 4nt) ==> 12/5 = 2.4
*/
double funcion()
{
    int i, j, n, myid;
    double *v, *w, *z, sv, sw, x, res = 0.0;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    if(myid == 0)
    {
        leer(&n, &v, &w, &z);
        MPI_Send(&n, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
        MPI_Send(&n, 1, MPI_INT, 2, 0, MPI_COMM_WORLD);
        MPI_Send(w, n, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
        MPI_Send(z, n, MPI_DOUBLE, 2, 0, MPI_COMM_WORLD);
        calcula_v(n, v); /*tarea 1*/
        MPI_Recv(&w, n, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        /*tarea 4 */
        for(j = 0; j < n; j++)
        {
            sv = 0;
            for(i = 0; i < n; i++) sv += v[i]*w[i];
            for(i = 0; i < n; i++) v[i] = sv*v[i];
        }
        MPI_Send(&sv, 1, MPI_DOUBLE, 2, 0, MPI_COMM_WORLD);
        MPI_Recv(&res, 1, MPI_DOUBLE, 2, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    else if (myid == 1)
    {
        MPI_Recv(&n, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        w = (double *) malloc(n*sizeof(double)); //usado para reservar memoria dinámica (*w o *z)
        MPI_Recv(w, n, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        calcula_w(n, w); /*tarea 2*/
        MPI_Send(w, n, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        MPI_Send(w, n, MPI_DOUBLE, 2, 0, MPI_COMM_WORLD);
        MPI_Recv(&res, 1, MPI_DOUBLE, 2, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    }
    else if (myid == 2)
    {
        MPI_Recv(&n, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        z = (double *) malloc(n*sizeof(double)); //usado para reservar memoria dinámica (*z o *w)
        MPI_Recv(z, n, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        calcula_z(n, z); /*tarea 3*/
        MPI_Recv(&w, n, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        /*tarea 5*/
        for(j = 0; j < n; j++)
        {
            sw = 0;
            for(i = 0; i < n; i++) sw += w[i]*z[i];
            for(i = 0; i < n; i++)  z[i] = sw*z[i];
        }
        MPI_Recv(&sv, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        /*tarea 6*/
        x  = sv+sw;
        for(j = 0; j < n; j++) res += x*z[j];
        MPI_Send(&res, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        MPI_Send(&res, 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
    }
    
    printf("%f\n", res);
    return res;    
}

//ejercicio 2.12
/**
 * coste secuencial:
 * N^2
 * coste paralelo:
 * taritmeticas = N^2/p 
 * tcomunicaciones = 2(p-1)(ts + N^2/p*tw)
 * ttotal = N^2/p + 2(p-1)(ts + N^2/p*tw)
 * speedup = N^2 / (N^2/p + 2(p-1)(ts + N^2/p*tw))
 * eficiencia = speedup/p ==> (N^2 / (N^2 + 2(p-1)(ts + N^2/p*tw))
*/
void ejercicio2_12()
{
    int i, j, myid, nprocs;
    double A[N][N], B[N][N];
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    int k = N/nprocs;


    if(myid == 0)
    {
        leer(A);   
    }
    MPI_Scatter(A, k*N, MPI_DOUBLE, B, k*N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    for(i = 0; i < k; i++)
    {
        for(j = 0; j < N; j++)
        {
            B[i][j] = B[i][j] * B[i][j];
        }
    }
    MPI_Gather(B, k*N, MPI_DOUBLE, A, k*N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}


//Ejercicio 2.07
/**
 * coste secuencial:
 * N^2
 * coste paralelo:
 * taritmeticas = N^2/p 
 * tcomunicaciones = (p-1)(ts + N^2/p*tw) + (p-1)(ts + tw)
 * ttotal = N^2/p + (p-1)(ts + N^2/p*tw) + (p-1)(ts + tw) ~= 2pts + N^2tw
 * speedup = N^2 / (N^2/p + 2pts + N^2tw)
 * eficiencia = speedup/p ==> (N^2 / (N^2 + 2pts + N^2tw)) ==> 1/1+ptw
*/
#include <math.h>
#define N 800
double ejercicio2_07(double A [][N], double Alocal [][N])
{
    int i, j, myid, nprocs;
    double s, nrm = 0.0, nrm_local;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    int k = N/nprocs; //es N por el bucle de las filas del codigo original

    MPI_Scatter(A, k*N, MPI_DOUBLE, Alocal, k*N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    for(i = 0; i < k; i++)
    {
        s = 0.0;
        for(j = 0; j < N; j++) s += fabs(Alocal[i][j]);
        if(s > nrm_local) nrm_local = s;
    }

    MPI_Reduce(&nrm_local, &nrm, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    return nrm;
}

//Ejercicio 3-16
int ejercicio3_16(float mat[M][N], float lmat[M][N/4])
{
    int myid;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Datatype col_cic;

    MPI_Type_vector(M*N/4, 1, 4, MPI_FLOAT, &col_cic);
    MPI_Type_commit(&col_cic);

    if(myid == 0)
    {
        MPI_Sendrecv(&mat[0][0], 1, col_cic, 0, 0, lmat, M*N/4, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        for(int i = 1; i < 4; i++)
        {
            MPI_Send(&mat[0][i], 1, col_cic, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);    
        }
    }
    else
    {
        MPI_Recv(lmat, M*N/4, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    MPI_Type_free(&col_cic);
    
    return 0;
}

//Ejercicio 3-15
void ejercicio3_15(double A[N][N])
{
    MPI_Datatype diagonal, antidiagonal;
    MPI_Type_vector(N, 1, N+1, MPI_DOUBLE, &diagonal);
    MPI_Type_vector(N, 1, N-1, MPI_DOUBLE, &antidiagonal);
    MPI_Type_commit(&diagonal);
    MPI_Type_commit(&antidiagonal);
    
    //apartado a
    MPI_Bcast(&A[0][0], 1, diagonal, 0, MPI_COMM_WORLD);
    MPI_Bcast(&A[0][N-1], 1, antidiagonal, 0, MPI_COMM_WORLD);

    MPI_Type_free(&diagonal);
    MPI_Type_free(&antidiagonal);
}

//apartado b
void ejercicio3_15_b(double A[N][N], double prin[N], double anti[N])
{
    int myid, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Datatype diagonal, antidiagonal;
    MPI_Type_vector(N, 1, N+1, MPI_DOUBLE, &diagonal);
    MPI_Type_vector(N, 1, N-1, MPI_DOUBLE, &antidiagonal);
    MPI_Type_commit(&diagonal);
    MPI_Type_commit(&antidiagonal);
    
    if(myid == 0)
    {
        MPI_Sendrecv(&A[0][0], 1, diagonal, 0, 0, prin, N, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Sendrecv(&A[0][N-1], 1, antidiagonal, 0, 0, anti, N, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        for(int i = 1; i < nprocs; i++)
        {
            MPI_Send(&A[0][i], 1, diagonal, i, 0, MPI_COMM_WORLD);    
            MPI_Send(&A[0][N-1], 1, antidiagonal, i, 0, MPI_COMM_WORLD);    
        }
    }
    else
    {
        MPI_Recv(prin, N, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(anti, N, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    
    

    MPI_Type_free(&diagonal);
    MPI_Type_free(&antidiagonal);
}