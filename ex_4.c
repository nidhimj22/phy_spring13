#include<stdio.h>
#include<stdlib.h>
#include<math.h>


int main()
{
	int i,ns,nsi,si,j,k,l;

	int n=100;    //no of particles
	double r[n],a[n],ao[n],v[n];double dt,ke=0.0,pe=0.0,ds,ds2,dsi6,fij,fjk,p,energy;
	float T_sum=0.0,T=0.0;
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
//    ns=20;
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

	fp = fopen ("velocity.txt", "w");
	fp2 = fopen ("energy.txt", "w");
	fp1 = fopen ("position.txt", "w");
	fp3 = fopen ("momentum.txt", "w");
	fp4 = fopen ("accerelation.txt", "w");
	fp5 = fopen ("temperature.txt", "w");

	for(i=0;i<10000;i++){
		if ((i<2950)||(i>3000))
{	        p =0.0;
	        pe = 0.0;
		ke =0.0;
//		force1f =0.0;
//		force2f = 0.0;
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
		               // ds= abs(ds);
				if (ds<0)
			       {printf("r is %f \n",ds);}
	    				ds2=ds*ds;
				dsi6=1.0/ds2/ds2/ds2;
				fjk=(48.0*dsi6-24.0)*dsi6/ds;

				//Force on first particle
				a[j]+=-fjk;
				a[k]+=+fjk;

				//L-J Potential energy
				pe+=(dsi6-1.00)*dsi6*4.0;

				//Integrating the equation of motion

				//Update velocities

			}
		}

		//Compute the kinetic energy

		for (l=0;l<n;l++)
		{
			v[l]=v[l]+0.5*(a[l]+ao[l])*dt;
			//v[k]=v[k]+0.5*(a[k]+ao[k])*dt;
			//Update the positions
			r[l]=r[l]+v[l]*dt+0.5*a[l]*dt*dt;
			if (r[l]<0)
			{		r[l]=r[l]+L;}

			if (r[l]>L)
			{	r[l]=r[l]-L;}

			//r[k]=r[k]+v[k]*dt+0.5*a[k]*dt*dt;
			ao[l]=a[l];
			//ao[k]=a[k];

			ke+=0.5*v[l]*v[l];
			p+=v[l];
			if ((i<1000)&&(i>=1))
			{v[l] = v[l]/sqrt(T);}

		}

		T_sum = T_sum +(ke/n);
		T = ke/n;


		//Collecting the output at specified interval

		if(i%si==0)fprintf(inp,"%f %lf %lf %lf %lf %lf %lf\n",dt*(float)i,r[0],r[1],ke,pe,ke+pe,v[0]);


		for (l=0;l<n;++l)
		{
	            fprintf(fp,"%6.2f, ",v[l]);
                    fprintf(fp1,"%6.2f,  ",r[l]);
		    fprintf(fp4,"%6.2f,  ",a[l]);
		}

                    fprintf(fp,"\n ");
        	    fprintf(fp1,"\n");
		    fprintf(fp4,"\n");

		    energy=pe+ke;

		    fprintf(fp2,"%f\n",energy);
		    fprintf(fp3,"%f\n",p);
		    fprintf(fp5,"%f\n",T);
		}
	}

	fclose(inp);

	return 0;
}
