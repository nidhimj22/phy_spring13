#include<stdio.h>
#include<stdlib.h>
#include<math.h>

int main()
{
    int i,ns,nsi,si;

    double dt,ke,pe,r,r2,ri6,forceij;


    /* Initialize all variable to
     * appropriate values
     */

    double	x1=0.00;//Initial position
    double	x2=2.00;//Initial position

    double	v1=0.0;//velocity
    double	v2=0.0;//velocity
    double	a1=0.0;//acceleration
    double	a2=0.0;//acceleration
    double	ao1=0.0;//old acceleration
    double	ao2=0.0;//old acceleration
    double energy = 0.0;
    float p=0.0;
    
        dt=0.001;//time step of integration
    
        //Time steps for integration - total time is unit of 10
        ns=(int)((1.0/dt)*10);

    //Sampling interval at which data
    //is recorded from the simulations, the interval is 0.1 time units
    si=(int)((1.0/dt)*0.01);

    FILE *inp=fopen("traj.txt","w");

    FILE* fp;
    FILE* fp1;
    FILE* fp2;
    FILE* fp3;
    FILE* fp4;
    FILE* fp5;
    FILE* fp6;

    fp = fopen ("velocity.txt", "w");
    fp2 = fopen ("energy.txt", "w");
    fp1 = fopen ("position.txt", "w");
    fp3 = fopen ("momentum.txt", "w");
    fp4 = fopen ("ke.txt", "w");
    fp5 = fopen ("pe.txt", "w");
    fp6 = fopen ("r.txt", "w");

    for(i=0;i<ns;i++){

        //Calculation of Lennard-Jones forces - potential energy
        r=x2-x1;
        r2=r*r;
        ri6=1.0/r2/r2/r2;
        forceij=(48.0*ri6-24.0)*ri6/r;

        //Force on first particle
        a1=-forceij;
        a2=+forceij;

        //L-J Potential energy
        pe=(ri6-1.00)*ri6*4.0;

        //Integrating the equation of motion
        //

        //Update velocities
        v1=v1+0.5*(a1+ao1)*dt;
        v2=v2+0.5*(a2+ao2)*dt;

        //Update the positions
        x1=x1+v1*dt+0.5*a1*dt*dt;
        x2=x2+v2*dt+0.5*a2*dt*dt;

        //Update the old accelerations with new
        ao1=a1;
        ao2=a2;

        //Compute the kinetic energy
        ke=0.5*v1*v1+0.5*v2*v2;

        //Collecting the output at specified interval

        if(i%si==0)fprintf(inp,"%f %lf %lf %lf %lf %lf %lf\n",dt*(float)i,x1,x2,ke,pe,ke+pe,v1);

        energy = pe+ ke;
        fprintf(fp2,"%3.4f\n",energy);
        fprintf(fp4,"%f\n",ke);
        fprintf(fp5,"%3.4f\n",pe);
        fprintf(fp3,"%f\n",p);
        fprintf(fp,"%f, %f\n",v1,v2);
        fprintf(fp1,"%f, %f\n",x1,x2);
        fprintf(fp6,"%20.18f,\n",r);

    }

    fclose(inp);

    return i;
}
