#include<stdio.h>
#include<stdlib.h>
#include<math.h>

int main()
{
    int i,ns,nsi,si,j,k,l;

    int n=10;    //no of particles
    double r[n],a[n],ao[n],v[n];
    double dt,ke=0.0,pe=0.0,ds=0.0,ds2=0.0,dsi6=0.0,fij=0.0,fjk=0.0,p=0.0,energy=0.0;

    /* Initialize all variable to	 * appropriate values	 */


    for (j=0;j<n;j++)
    {
        v[j]= 0.0;
        a[j]=0.0;
        ao[j] =0.0;
        r[j]= j*1.0+0.5;

    }

    float rho =0.6;
    float scaling_factor = 1/rho;
    float L = scaling_factor * n;


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

    fp = fopen ("velocity.txt", "w");
    fp2 = fopen ("energy.txt", "w");
    fp1 = fopen ("position.txt", "w");
    fp3 = fopen ("ke.txt", "w");
    fp4 = fopen ("accerelation.txt", "w");

    for(i=0;i<2955;i++){

        p =0.0;
        pe = 0.0;
        ke =0.0;
        for (l=0;l<n;l++)
        {
            ao[l]=a[l];
            a[l]=0.0;
        }
        for (j=0;j<n;j++)
        {
            for (k=j+1;k<n;k++)
            {

                ds=r[k]-r[j];

                if (ds<0)
                {ds=-ds;}

                if (ds> 0.5 * L)
                {ds = L-ds;}

                if (ds<0)
                {ds=-ds;}

                ds2=ds*ds;
                dsi6=1.0/ds2/ds2/ds2;

                fjk=(48.0*dsi6-24.0)*dsi6/ds;
                printf("dsi is %f\n",dsi6);
                //Force on first particle
                a[j]+=-fjk;
                a[k]+=+fjk;

                //L-J Potential energy
                pe+=(dsi6-1.00)*dsi6*4.0;


            }
        }

        //Compute the kinetic energy

        for (l=0;l<n;l++)
        {
            v[l]=v[l]+0.5*(a[l]+ao[l])*dt;

            //Update the positions

            r[l]=r[l]+v[l]*dt+0.5*a[l]*dt*dt;

            if (r[l]<0)
            {r[l]=r[l]+L;}

            if (r[l]>L)
            {r[l]=r[l]-L;}


            ao[l]=a[l];

            ke+=0.5*v[l]*v[l];
            p+=v[l];

        }

        //Collecting the output at specified interval

        if(i%si==0)fprintf(inp,"%f %lf %lf %lf %lf %lf %lf\n",dt*(float)i,r[0],r[1],ke,pe,ke+pe,v[0]);


        for (l=0;l<n;++l)
        {
            fprintf(fp,"%3.4f, ",v[l]);
            fprintf(fp1,"%3.4f,  ",r[l]);
            fprintf(fp4,"%3.4f,  ",a[l]);
        }

        fprintf(fp,"\n ");
        fprintf(fp1,"\n");
        fprintf(fp4,"\n");

        energy=pe+ke;

        fprintf(fp2,"%3.4f\n",energy);
        fprintf(fp3,"%3.4f\n",ke);

    }

    fclose(inp);

    return i;
}
