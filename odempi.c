#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <mpi.h>
//#define num 1000 //Number of asteroids
#define WORKTAG 1
#define DIETAG 2

void initialize(double array[], int num)
{
   int i;
   for( i = 0; i < num; i++)
   {
      array[i] = 0.;
   }
}

double arrprint(double array[], int num)
{
    int k;
    for(k = 0; k < num; k++)
    {
        printf("arr of %d %g\n",k, array[k]); 
    }
}

double RHS_J_1 (double m3, //Dimensionless Newton's 2st Law equation for Jupiter
            double x1,  double y1, double ux1, double uy1)
{
    return ux1;
}

double RHS_J_2 (double m3, 
            double x1,  double y1, double ux1, double uy1)
{
    return -m3*x1/(pow(x1*x1 + y1*y1, 1.5));
}

double RHS_J_3 (double m3, 
            double x1,  double y1, double ux1, double uy1)
{
    return uy1;
}

double RHS_J_4 (double m3, 
            double x1,  double y1, double ux1, double uy1)
{
    return -m3*y1/(pow(x1*x1 + y1*y1, 1.5));
}

double RHS_a_1 (double m1,  double m3, //Dimensionless Newton's 2st Law equation for asteroids
            double x1,  double x2,
            double y1,  double y2,
            double ux2, double uy2)
{
        return ux2;
}

double RHS_a_2 (double m1,  double m3,
            double x1,  double x2,
            double y1,  double y2,
            double ux2, double uy2)
{
        return m1*(x1 - x2)/(pow((x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1), 1.5)) - m3*x2/(pow(x2*x2 + y2*y2, 1.5));
}

double RHS_a_3 (double m1,  double m3,
            double x1,  double x2,
            double y1,  double y2,
            double ux2, double uy2)
{
        return uy2;
}

double RHS_a_4 (double m1,  double m3,
            double x1,  double x2,
            double y1,  double y2,
            double ux2, double uy2)
{
        return m1*(y1 - y2)/(pow((x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1), 1.5)) - m3*y2/(pow(x2*x2 + y2*y2, 1.5));
}

double R (double x, double y) //Radial distance of asteroids
{
    return sqrt(x*x + y*y);
}


double * RK4_Jupiter(int N, int m_Sun, double dt,
    double k1_x_J[],   double k2_x_J[],   double k3_x_J[],   double k4_x_J[],
    double k1_y_J[],   double k2_y_J[],   double k3_y_J[],   double k4_y_J[],
    double x_J[], double y_J[], double ux_J[], double uy_J[])
{
    int i;
    double k1_ux_J;  double k2_ux_J;  double k3_ux_J;  double k4_ux_J;
    double k1_uy_J;  double k2_uy_J;  double k3_uy_J;  double k4_uy_J;
    for(i = 0; i < N; i++)
    {

        k1_x_J[i] = dt * RHS_J_1(m_Sun,  x_J[i],  y_J[i],  ux_J[i],  uy_J[i]);
        k1_ux_J   = dt * RHS_J_2(m_Sun,  x_J[i],  y_J[i],  ux_J[i],  uy_J[i]);
        k1_y_J[i] = dt * RHS_J_3(m_Sun,  x_J[i],  y_J[i],  ux_J[i],  uy_J[i]);
        k1_uy_J   = dt * RHS_J_4(m_Sun,  x_J[i],  y_J[i],  ux_J[i],  uy_J[i]);


        k2_x_J[i] = dt * RHS_J_1(m_Sun,  x_J[i] + 0.5 * k1_x_J[i],  y_J[i] + 0.5 * k1_y_J[i],  ux_J[i] + 0.5 * k1_ux_J,  uy_J[i] + 0.5 * k1_uy_J);
        k2_ux_J   = dt * RHS_J_2(m_Sun,  x_J[i] + 0.5 * k1_x_J[i],  y_J[i] + 0.5 * k1_y_J[i],  ux_J[i] + 0.5 * k1_ux_J,  uy_J[i] + 0.5 * k1_uy_J);
        k2_y_J[i] = dt * RHS_J_3(m_Sun,  x_J[i] + 0.5 * k1_x_J[i],  y_J[i] + 0.5 * k1_y_J[i],  ux_J[i] + 0.5 * k1_ux_J,  uy_J[i] + 0.5 * k1_uy_J);
        k2_uy_J   = dt * RHS_J_4(m_Sun,  x_J[i] + 0.5 * k1_x_J[i],  y_J[i] + 0.5 * k1_y_J[i],  ux_J[i] + 0.5 * k1_ux_J,  uy_J[i] + 0.5 * k1_uy_J);
        

        k3_x_J[i] = dt * RHS_J_1(m_Sun,  x_J[i] + 0.5 * k2_x_J[i],  y_J[i] + 0.5 * k2_y_J[i],  ux_J[i] + 0.5 * k2_ux_J,  uy_J[i] + 0.5 * k2_uy_J);
        k3_ux_J   = dt * RHS_J_2(m_Sun,  x_J[i] + 0.5 * k2_x_J[i],  y_J[i] + 0.5 * k2_y_J[i],  ux_J[i] + 0.5 * k2_ux_J,  uy_J[i] + 0.5 * k2_uy_J);
        k3_y_J[i] = dt * RHS_J_3(m_Sun,  x_J[i] + 0.5 * k2_x_J[i],  y_J[i] + 0.5 * k2_y_J[i],  ux_J[i] + 0.5 * k2_ux_J,  uy_J[i] + 0.5 * k2_uy_J);
        k3_uy_J   = dt * RHS_J_4(m_Sun,  x_J[i] + 0.5 * k2_x_J[i],  y_J[i] + 0.5 * k2_y_J[i],  ux_J[i] + 0.5 * k2_ux_J,  uy_J[i] + 0.5 * k2_uy_J);


        k4_x_J[i] = dt * RHS_J_1(m_Sun,  x_J[i] + k3_x_J[i],  y_J[i] + k3_y_J[i],  ux_J[i] + k3_ux_J,  uy_J[i] + k3_uy_J);
        k4_ux_J   = dt * RHS_J_2(m_Sun,  x_J[i] + k3_x_J[i],  y_J[i] + k3_y_J[i],  ux_J[i] + k3_ux_J,  uy_J[i] + k3_uy_J);
        k4_y_J[i] = dt * RHS_J_3(m_Sun,  x_J[i] + k3_x_J[i],  y_J[i] + k3_y_J[i],  ux_J[i] + k3_ux_J,  uy_J[i] + k3_uy_J);
        k4_uy_J   = dt * RHS_J_4(m_Sun,  x_J[i] + k3_x_J[i],  y_J[i] + k3_y_J[i],  ux_J[i] + k3_ux_J,  uy_J[i] + k3_uy_J);


        x_J[i+1]  = x_J[i]  + (1./6.) * k1_x_J[i] + (1./3.) * k2_x_J[i] + (1./3.) * k3_x_J[i] + (1./6.) * k4_x_J[i];
        ux_J[i+1] = ux_J[i] + (1./6.) * k1_ux_J   + (1./3.) * k2_ux_J   + (1./3.) * k3_ux_J   + (1./6.) * k4_ux_J;
        y_J[i+1]  = y_J[i]  + (1./6.) * k1_y_J[i] + (1./3.) * k2_y_J[i] + (1./3.) * k3_y_J[i] + (1./6.) * k4_y_J[i];
        uy_J[i+1] = uy_J[i] + (1./6.) * k1_uy_J   + (1./3.) * k2_uy_J   + (1./3.) * k3_uy_J   + (1./6.) * k4_uy_J;
        //Runge-Kutta method solving Jupiter's motion
    }
    return x_J, ux_J, y_J, uy_J;
}

double RK4_asteroid(int N, double m_Sun, double m_J, double dt, int k,
    double x_J[], double y_J[N],
    double k1_x_J[], double k1_y_J[], 
    double k2_x_J[], double k2_y_J[], 
    double k3_x_J[], double k3_y_J[],
    double x_a[], double y_a[], 
    double ux_a[], double uy_a[])
{
    int i;
            double k1_x_A;  double k2_x_A;  double k3_x_A;  double k4_x_A;
            double k1_ux_A; double k2_ux_A; double k3_ux_A; double k4_ux_A;
            double k1_y_A;  double k2_y_A;  double k3_y_A;  double k4_y_A;
            double k1_uy_A; double k2_uy_A; double k3_uy_A; double k4_uy_A;
            double x[N]; double y[N];
            double ux[N]; double uy[N];
            double x_0; double y_0;
         double R_a; double R_p; double R_0; double R_1; double R_2;

        x[0]  = x_a[k]; 
        y[0]  = y_a[k]; 
        ux[0] = ux_a[k]; 
        uy[0] = uy_a[k];
    for(i = 0; i < N; i++)
    {
        k1_x_A  = dt * RHS_a_1(m_J,  m_Sun,  x_J[i],  x[i],  y_J[i],  y[i],  ux[i],  uy[i]);
        k1_ux_A = dt * RHS_a_2(m_J,  m_Sun,  x_J[i],  x[i],  y_J[i],  y[i],  ux[i],  uy[i]);
        k1_y_A  = dt * RHS_a_3(m_J,  m_Sun,  x_J[i],  x[i],  y_J[i],  y[i],  ux[i],  uy[i]);
        k1_uy_A = dt * RHS_a_4(m_J,  m_Sun,  x_J[i],  x[i],  y_J[i],  y[i],  ux[i],  uy[i]);
            
            
        k2_x_A  = dt * RHS_a_1(m_J,  m_Sun,  x_J[i] + 0.5 * k1_x_J[i],  x[i] + 0.5 * k1_x_A,  y_J[i] + 0.5 * k1_y_J[i],  y[i] + 0.5 * k1_y_A,  ux[i] + 0.5 * k1_ux_A,  uy[i] + 0.5 * k1_uy_A);
        k2_ux_A = dt * RHS_a_2(m_J,  m_Sun,  x_J[i] + 0.5 * k1_x_J[i],  x[i] + 0.5 * k1_x_A,  y_J[i] + 0.5 * k1_y_J[i],  y[i] + 0.5 * k1_y_A,  ux[i] + 0.5 * k1_ux_A,  uy[i] + 0.5 * k1_uy_A);
        k2_y_A  = dt * RHS_a_3(m_J,  m_Sun,  x_J[i] + 0.5 * k1_x_J[i],  x[i] + 0.5 * k1_x_A,  y_J[i] + 0.5 * k1_y_J[i],  y[i] + 0.5 * k1_y_A,  ux[i] + 0.5 * k1_ux_A,  uy[i] + 0.5 * k1_uy_A);
        k2_uy_A = dt * RHS_a_4(m_J,  m_Sun,  x_J[i] + 0.5 * k1_x_J[i],  x[i] + 0.5 * k1_x_A,  y_J[i] + 0.5 * k1_y_J[i],  y[i] + 0.5 * k1_y_A,  ux[i] + 0.5 * k1_ux_A,  uy[i] + 0.5 * k1_uy_A);

            
        k3_x_A  = dt * RHS_a_1(m_J,  m_Sun,  x_J[i] + 0.5 * k2_x_J[i],  x[i] + 0.5 * k2_x_A,  y_J[i] + 0.5 * k2_y_J[i],  y[i] + 0.5 * k2_y_A,  ux[i] + 0.5 * k2_ux_A,  uy[i] + 0.5 * k2_uy_A);
        k3_ux_A = dt * RHS_a_2(m_J,  m_Sun,  x_J[i] + 0.5 * k2_x_J[i],  x[i] + 0.5 * k2_x_A,  y_J[i] + 0.5 * k2_y_J[i],  y[i] + 0.5 * k2_y_A,  ux[i] + 0.5 * k2_ux_A,  uy[i] + 0.5 * k2_uy_A);
        k3_y_A  = dt * RHS_a_3(m_J,  m_Sun,  x_J[i] + 0.5 * k2_x_J[i],  x[i] + 0.5 * k2_x_A,  y_J[i] + 0.5 * k2_y_J[i],  y[i] + 0.5 * k2_y_A,  ux[i] + 0.5 * k2_ux_A,  uy[i] + 0.5 * k2_uy_A);
        k3_uy_A = dt * RHS_a_4(m_J,  m_Sun,  x_J[i] + 0.5 * k2_x_J[i],  x[i] + 0.5 * k2_x_A,  y_J[i] + 0.5 * k2_y_J[i],  y[i] + 0.5 * k2_y_A,  ux[i] + 0.5 * k2_ux_A,  uy[i] + 0.5 * k2_uy_A);


        k4_x_A  = dt * RHS_a_1(m_J,  m_Sun,  x_J[i] + k3_x_J[i],  x[i] + k3_x_A,  y_J[i] + k3_y_J[i],  y[i] + k3_y_A,  ux[i] + k3_ux_A,  uy[i] + k3_uy_A);
        k4_ux_A = dt * RHS_a_2(m_J,  m_Sun,  x_J[i] + k3_x_J[i],  x[i] + k3_x_A,  y_J[i] + k3_y_J[i],  y[i] + k3_y_A,  ux[i] + k3_ux_A,  uy[i] + k3_uy_A);
        k4_y_A  = dt * RHS_a_3(m_J,  m_Sun,  x_J[i] + k3_x_J[i],  x[i] + k3_x_A,  y_J[i] + k3_y_J[i],  y[i] + k3_y_A,  ux[i] + k3_ux_A,  uy[i] + k3_uy_A);
        k4_uy_A = dt * RHS_a_4(m_J,  m_Sun,  x_J[i] + k3_x_J[i],  x[i] + k3_x_A,  y_J[i] + k3_y_J[i],  y[i] + k3_y_A,  ux[i] + k3_ux_A,  uy[i] + k3_uy_A);

                
        R_0 = R(x_0, y_0);
        R_1 = R(x[i], y[i]);
        x_0 = x[i]; 
        y_0 = y[i];

        x[i+1]  = x[i]  + (1./6.) * k1_x_A  + (1./3.) * k2_x_A  + (1./3.) * k3_x_A  + (1./6.) * k4_x_A;
        ux[i+1] = ux[i] + (1./6.) * k1_ux_A + (1./3.) * k2_ux_A + (1./3.) * k3_ux_A + (1./6.) * k4_ux_A;
        y[i+1]  = y[i]  + (1./6.) * k1_y_A  + (1./3.) * k2_y_A  + (1./3.) * k3_y_A  + (1./6.) * k4_y_A;
        uy[i+1] = uy[i] + (1./6.) * k1_uy_A + (1./3.) * k2_uy_A + (1./3.) * k3_uy_A + (1./6.) * k4_uy_A;

        //Runge-Kutta method solving for asteroids' motion


        R_2 = R(x[i+1], y[i+1]);
            
        if(R_0 < R_1 && R_2 < R_1)
            R_a = R_1;
        //finding aphelion
        else if(R_0 > R_1 && R_2 > R_1)
            R_p = R_1;
        //finding perihelion
    }
    return 0.5*(R_a + R_p);
}

void master(double R_output[], int num)
{
      int numProcs;
	
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs); 
  
    MPI_Status status;
  
      int kcount = 0;
      int workerRank;
      
      for(workerRank = 1; workerRank < numProcs; workerRank++)
      {   
         MPI_Send(&kcount, 1, MPI_INT, workerRank, WORKTAG, MPI_COMM_WORLD); 
         kcount++;
      }
      // receive the work done by workers
      double termAnswer;
      double finalAnswer = 0;
      
      while(kcount < num)
      {    
      	  double term[2];
          // Receive the results from workers and accumulate the result
          MPI_Recv(&term, 2, MPI_DOUBLE, MPI_ANY_SOURCE, WORKTAG, MPI_COMM_WORLD, &status);
          int rcount  = (int)(term[1]);
         
          R_output[rcount] = term[0];
          workerRank = status.MPI_SOURCE;

          MPI_Send(&kcount, 1, MPI_INT, workerRank, WORKTAG, MPI_COMM_WORLD); 
          kcount++;
      }
        
      // Send (tag = DIETAG) to all workers
      for(workerRank = 1; workerRank < numProcs; workerRank++)
      {
        // sending garbage value because it will be ignored at worker
        kcount = -111;  
        MPI_Send(&kcount, 1, MPI_INT, workerRank, DIETAG, MPI_COMM_WORLD); 
      }
     
      // Do pending receives for outstanding messages from workers
      for(workerRank = 1; workerRank < numProcs; workerRank++)
      {
        double term[2];
        MPI_Recv(&term, 2, MPI_DOUBLE, workerRank, WORKTAG, MPI_COMM_WORLD, &status);
        int rcount  = (int)(term[1]);
        R_output[rcount] = term[0];
      }
}
  
void worker(int N, double m_Sun, double m_J, double dt, int k,
    double x_J[], double y_J[],
    double k1_x_J[], double k1_y_J[], 
    double k2_x_J[], double k2_y_J[], 
    double k3_x_J[], double k3_y_J[],
    double x_a[], double y_a[], 
    double ux_a[], double uy_a[])
{
  MPI_Status status;
    int val; 
        
    // worker keeps looping; waiting for work items from master process, 
    // unless it gets termination message
    while(1)
    { 
       int kcount;
       MPI_Recv(&kcount, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

       if(status.MPI_TAG == DIETAG)
       {
           //printf("TERMINATING. BYE \n");
           fflush(stdout);
           return;
       }
       else    // (status.MPI_TAG = WORKTAG)
       {
         // evaluateTerm(int coefficient, int degree, double x)
         double answer = RK4_asteroid(N,  m_Sun, m_J,  dt,  kcount,
                x_J, y_J,
                k1_x_J,  k1_y_J, 
                k2_x_J,  k2_y_J, 
                k3_x_J,  k3_y_J,
                x_a,  y_a, 
                ux_a,  uy_a);
         double term[2];
         term[0] = answer;
         term[1] = kcount;
         MPI_Send(&term, 2, MPI_DOUBLE, 0, WORKTAG, MPI_COMM_WORLD);
         fflush(stdout);
       }
    } 
}

int main(int argc, char **argv)
{
    int num;;
    int AllTime;

    FILE *fp;
    fp = fopen("num_of_astr_config.txt","r");
    if(fp == NULL)  //If file failed to open
    {
        printf("Opening the file failed.Exiting...");
        return -1;
    }

    fscanf(fp, "%d", &num);
    fclose(fp);

    FILE *fpp;
    fpp = fopen("last_time_config.txt","r");
    if(fpp == NULL)  //If file failed to open
    {
        printf("Opening the file failed.Exiting...");
        return -1;
    }

    fscanf(fpp, "%d", &AllTime);
    fclose(fpp);

    //printf("%d %d\n", num, AllTime);


    double m_Sun = 4*pow(M_PI, 2); //Dimensionless mass of the Sun
    double m_J = 0.0376807; //Dimensionless mass of Jupiter
    double t = 0;
    int i, k, T;    
    double dt = 1./128.;
    int N = AllTime/dt; //time for all

    double k1_x_J[N],   k2_x_J[N],   k3_x_J[N],   k4_x_J[N];
    double k1_y_J[N],   k2_y_J[N],   k3_y_J[N],   k4_y_J[N];

    double m_a[num], x_a[num], y_a[num], ux_a[num], uy_a[num];
    double x_0[num], y_0[num];
    double x_J[N], y_J[N];
    double ux_J[N], uy_J[N];
//    double R_output;

    x_J[0] = 5.4588;
    y_J[0] = 0;
    ux_J[0] = 0;
    uy_J[0] = 2.62267;

    FILE *file;
    file = fopen("IV.txt", "r");
    
    if(file == NULL)
    {
        printf("Unable to open file\n");
        return -1;
    }
    
    for(i = 0; i < num; i++)
    {
        fscanf(file, "%lg %lg %lg %lg %lg", &m_a[i], &x_a[i], &y_a[i], &ux_a[i], &uy_a[i]);
        //R_output[i] = 0;
    }
    fclose(file);
    
    //reading initial values from the outcome of IV

    double *p;
    int temp_i;

    p = RK4_Jupiter( N,  m_Sun, dt,
                k1_x_J,   k2_x_J,   k3_x_J,   k4_x_J,
                k1_y_J,   k2_y_J,   k3_y_J,   k4_y_J,
                x_J, y_J,
                ux_J, uy_J);
    
    double *R_output = (double *)malloc(sizeof(double) * num);
    //initialize(R_output, num);
    
    
    MPI_Init(&argc, &argv);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        
    if (rank == 0) 
    {
       master(R_output,num);
    } 
    else 
    {
       worker(N,  m_Sun, m_J,  dt,  k,
                x_J, y_J,
                k1_x_J,  k1_y_J, 
                k2_x_J,  k2_y_J, 
                k3_x_J,  k3_y_J,
                x_a,  y_a, 
                ux_a,  uy_a);
    }
    
    MPI_Finalize();
    FILE *fpWrite=fopen("odempi.txt","w");
  
    for(k = 0; k < num; k++)
    {
        fprintf(fpWrite,"%lg\n", R_output[k]); 
    }   
	  free(R_output);
    fclose(fpWrite);

    printf("final");
    
    return 0;
}
