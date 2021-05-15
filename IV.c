#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <string.h>
//#define n 100
//n is the number of the asteroid

int main()
{
    FILE *fp;
    int n;
    fp = fopen("num_of_astr_config.txt","r");
    if(fp==NULL)  //If file failed to open
    {
        printf("Opening the file failed.Exiting...");
        return -1;
    }

    fscanf(fp, "%d", &n);
    fclose(fp);

    int i;
    double R[n], theta[n], v[n], m[n], x[n], y[n], ux[n], uy[n];
    double e[n];
    //srand(time(NULL));
    
    for(i = 0; i <= n - 1; i++)
    {
        m[i] = 1.8648*1e-8*(rand()/(RAND_MAX + 1.0)) + 2.7794*1e-18;
        //random mass of asteroid

        R[i] = 2 + i*(1.5/(n - 1));
        //random initial radial distance
        
        theta[i] = 2*M_PI*(rand()/(RAND_MAX + 1.0));
        //random initial angular position
        
        theta[i] = 0*M_PI;
        //set an identical initial angular position (if possible)

        e[i] = 0;
    
    }
    //R[0]=3.1;R[1]=3.2;R[2]=3.3;R[3]=3.4;R[4]=3.5;
    
    for(i = 0; i <= n - 1; i++)
    {
        v[i] = sqrt((4*pow(M_PI, 2))*(1 - e[i])/R[i]);
    }
    
    //calculate initial velocity for aphelion
    for(i = 0; i <= n - 1; i++)
    {
        x[i] = R[i]*cos(theta[i]);
        //initial x position
        
        y[i] = R[i]*sin(theta[i]);
        //initial y position
        
        ux[i] = -v[i]*sin(theta[i]);
        //initial velocity in x direction
        
        uy[i] = v[i]*cos(theta[i]);
        //initial velocity in y direction
    }
    
    FILE *fpWrite=fopen("IV.txt","w");

    for(i = 0; i <= n - 1; i++)
    {
        fprintf(fpWrite,"%lg %lg %lg %lg %lg\n", m[i], x[i], y[i], ux[i], uy[i]); 

        //printf("%g %g %g %g %g\n", m[i], x[i], y[i], ux[i], uy[i]);
    }
    fclose(fpWrite);
    printf("IV done \n");
    //print mass, initial x,y position, and initial velocity in x,y direction
    return 0;
}
