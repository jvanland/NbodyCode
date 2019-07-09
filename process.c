/* Take a list of positions and velocities for the two
10 Msun black holes.  Using this, compute a and e
for the motion of the center of mass around the 10^5 Msun
black holes, a and e for the binary, and the relative
inclination of the binary to its superorbital plane. 
Requires fulllist.dat
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI 3.141592653589793
#define G 39.47841760435743 /* Newton's constant for AU,yr,Msun*/

void aei(double Mtot, double x, double y, double z, double vx, 
double vy, double vz, double *a, double *e, double *Lx, double *Ly, 
double *Lz);

int main(int argc, char *argv[ ])
{
  if(argc != 3) {
    printf("Usage: %s M_SMBH data_file\n", argv[0]);
    return 1;
  }

   double M = atof(argv[1]);
   char *Initfile = argv[2];
   int i,N=0;
   double m1,m2,mtot,t,t1,t2,dt,temp;
   double x1,y1,z1,vx1,vy1,vz1;
   double x2,y2,z2,vx2,vy2,vz2;
   double xcm,ycm,zcm,vxcm,vycm,vzcm;
   double abin,ebin,Lxbin,Lybin,Lzbin;
   double acm,ecm,Lxcm,Lycm,Lzcm;
   double mut_inc,com_inc;
   FILE *data;
   char buf[5000];

   data=fopen(Initfile,"r");
   while((fgets(buf,5000,data)) != NULL) {
     N++;
   }
   rewind(data);

   t1 = t2 = 0;
   for (i=0; i<N/2; i++)
   {
      fscanf(data,"%lg %lg %lg %lg %lg %lg %lg",
         &temp,&x1,&y1,&z1,&vx1,&vy1,&vz1);
      fscanf(data,"%lg %lg %lg %lg %lg %lg %lg",
         &temp,&x2,&y2,&z2,&vx2,&vy2,&vz2);
      dt = t2-t1;
      if (i == 0) t = temp;
      else if (i == 1) t = t1 = temp;
      else if (i == 2) t = t2 = temp;
      else t += dt; 
      m1 = m2 = 10.0;
      xcm=(m1*x1+m2*x2)/(m1+m2);
      ycm=(m1*y1+m2*y2)/(m1+m2);
      zcm=(m1*z1+m2*z2)/(m1+m2);
      vxcm=(m1*vx1+m2*vx2)/(m1+m2);
      vycm=(m1*vy1+m2*vy2)/(m1+m2);
      vzcm=(m1*vz1+m2*vz2)/(m1+m2);
      x1-=xcm;
      x2-=xcm;
      y1-=ycm;
      y2-=ycm;
      z1-=zcm;
      z2-=zcm;
      vx1-=vxcm;
      vx2-=vxcm;
      vy1-=vycm;
      vy2-=vycm;
      vz1-=vzcm;
      vz2-=vzcm;
      mtot=m1+m2;
      aei(mtot,x1-x2,y1-y2,z1-z2,vx1-vx2,vy1-vy2,vz1-vz2,&abin,&ebin,&Lxbin,&Lybin,&Lzbin);
      mtot=M;
      aei(mtot,xcm,ycm,zcm,vxcm,vycm,vzcm,&acm,&ecm,&Lxcm,&Lycm,&Lzcm);
      mut_inc=acos(Lxbin*Lxcm+Lybin*Lycm+Lzbin*Lzcm)*180.0/PI;
      com_inc=acos(Lzcm)*180.0/PI;
      printf("%.16g %.16g %.16g %.16g %.16g %.16g %.16g\n",t,abin,ebin,mut_inc,acm,ecm,com_inc);
   }
   fclose(data);
}

void aei(double Mtot, double x, double y, double z, double vx, double vy, 
double vz, double *a, double *e, double *Lx, double *Ly, double *Lz)
{
/* Given position and velocity, compute semimajor axis, eccentricity,
and inclination */

   double dist,vel,energy,L;

   dist=sqrt(x*x+y*y+z*z);
   vel=sqrt(vx*vx+vy*vy+vz*vz);
   energy=0.5*vel*vel-G*Mtot/dist;
   *a=-G*Mtot/(2.0*energy);
   *Lx=y*vz-z*vy;
   *Ly=z*vx-x*vz;
   *Lz=x*vy-y*vx;
   L=sqrt((*Lx)*(*Lx)+(*Ly)*(*Ly)+(*Lz)*(*Lz));
   *Lx/=L;
   *Ly/=L;
   *Lz/=L;
   *e=sqrt(1.0-L*L/(G*Mtot*(*a)));
}
