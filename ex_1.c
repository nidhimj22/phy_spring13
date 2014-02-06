#include<stdio.h>
int main(void)
{
    FILE* fp;
    FILE* fp1;
    FILE* fp2;
    FILE* fp3;
    fp = fopen ("velocity.txt", "w");
    fp2 = fopen ("energy.txt", "w");  
    fp1 = fopen ("position.txt", "w");  
    fp3 = fopen ("momentum.txt", "w");  

    int n=0; // the number of time steps involved , t=n*dt
    float x1=0,x2=1,v1=0,v2=0;//x1, x2, v1, v2 are  the positions and velocities of the 2 particles
    int  m1=1,m2=1; // masses of the 2 bodies=unity
    float force1p=0.0,force2p=0.0,force1f=0.0,force2f=0.0;

    int k=100;
    float x0=0.5;
    //they are connected with a spring of force constant k and its natural length is x0
    int time=0;
    float dt=0.001;

    float ke=0,pe=0,energy=0,p=0;
    // energys and momentum defined
    //the spring is initially elongated 
    force1p=k*((x2-x1)-x0);
    force2p=- force1p;
    // present force on particle 1 and particle 2
    // Using Newton's third law , force on the 2 masses are opposite to each other in direction and equal in magnitude    


    for(n=0,time=0;time<5;n++,time=time+dt)
    {
        //using verlet equation to find the positions of the particles
        // as masses is unity , forces= acceleration 
        x1=x1+v1*dt+ (0.5*dt*dt*force1p);
        x2=x2+v2*dt+ (0.5*dt*dt*force2p);
        //future forces 
        force1f=k*((x2-x1)-x0);
        force2f=-force1f;
        //calculating velocity using verlet equation
        v1=v1+0.5*dt*(force1p+force1f);
        v2=v2+0.5*dt*(force2p+force2f);
        // mementum
        p=m1*v1+m2*v2;

        //energy
        ke=0.5*((v1*v1)+(v2*v2));
        pe=0.5*k*((x2-x1)-x0)*((x2-x1)-x0);
        energy=pe+ke;
        fprintf(fp2,"%f\n",energy);
        fprintf(fp3,"%f\n",p);
        fprintf(fp,"%f, %f\n",v1,v2);
        fprintf(fp1,"%f, %f\n",x1,x2);

        // present forces = future forces 
        force1p=force1f;
        force2p=force2f;
    }
    return 0;

}
