/* Reads in the output from OUT3 nbody6.
Follows orbital elements of 1 object
N_tot = first[1], Time = rec[0], N_pairs = rec[1] */

#include<stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

#define G 1.0 /* Newton's constant */
#define PI 3.141592653589793  /* pi to unnecessary precision */

typedef struct
{
   float x;     /* x component of vector. */
   float y;     /* y component of vector. */
   float z;     /* z component of vector. */
} vector;

void zerovec(vector *r);    /* Set components of vector to zero */
vector vecadd(vector a, vector b);   /* Vector addition a+b. */
vector vecsub(vector a, vector b);   /* Vector subtraction a-b. */
vector vecmult(float d, vector a);  /* Multiplication: da. */
vector vecdiv(vector a, float d);  /* Division: a/d. */
vector crossprod(vector a, vector b); /* Cross product, axb */
float magnitude(vector a);  /* Magnitude of vector a */

int main(int argc, char *argv[ ])
{
  if(argc != 8) {
    printf("Usage: %s Infile Nobj Nout Rconv Vconv Tconv Obj\n", argv[0]);
    return 1;
  }
        int i,j;
	int N_esc,N_pairs,N_tot;
	char *data = argv[1];
	int N_obj = atoi(argv[2]);
	int N_out = atoi(argv[3]);
	double Rconv = atof(argv[4]);
	double Vconv = atof(argv[5]);
	double Tconv = atof(argv[6]);
	int Obj = atoi(argv[7]);
	FILE *infile;
	int first[7];
	float rec[20];
	float *m,*pos,*vel;
	int *name;
	vector r1,r2,rcent,v1,v2,vcent;
	float mass1,mass2,dist,d,vmag;

	m = malloc(4);
	pos = malloc(4);
	vel = malloc(4);
	name = malloc(4);

       	infile=fopen(data,"rb");
       	if (!infile)
       	{
       		printf("Unable to open file!\n");
       		return 1;
       	}
	for (j=0;j<3;j++)
	{
	        // Read in particle independent data
		fread(&first,sizeof(first),1,infile);
		fread(&rec,sizeof(rec),1,infile);
		// Determine N_esc,N_pairs,N_tot
		N_tot = first[1];
		N_pairs = (int) rec[1];
		N_esc = N_obj - (N_tot - N_pairs);
		// Adjust size of arrays
		m = (float *) realloc(m,N_tot*sizeof(float));
		name = (int *) realloc(name,N_tot*sizeof(int));
		pos = (float *) realloc(pos,3*N_tot*sizeof(float));
		vel = (float *) realloc(vel,3*N_tot*sizeof(float));
		// Read in particle data
       		fread(m,N_tot*sizeof(float),1,infile);
       		fread(pos,3*N_tot*sizeof(float),1,infile);
       		fread(vel,3*N_tot*sizeof(float),1,infile);
       		fread(name,N_tot*sizeof(int),1,infile);
		fseek(infile,4,SEEK_CUR);
		zerovec(&r2);
		for (i=0;i<N_tot;i++)
		  {
		    if(name[i] == 1) {
		      mass1 = m[i];
		      r1.x = pos[3*i];
		      r1.y = pos[3*i+1];
		      r1.z = pos[3*i+2];
		      v1.x = vel[3*i];
		      v1.y = vel[3*i+1];
		      v1.z = vel[3*i+2];
		    }
		    if(name[i] == Obj) {
		      mass2 = m[i];
		      r2.x = pos[3*i];
		      r2.y = pos[3*i+1];
		      r2.z = pos[3*i+2];
		      v2.x = vel[3*i];
		      v2.y = vel[3*i+1];
		      v2.z = vel[3*i+2];
		    }
		  }
		if(magnitude(r2) == 0) {
		  printf("Object Escaped!\n");
		  return 1;
		}
		rcent = vecsub(r2,r1);
		vcent = vecsub(v2,v1);
		/* Conversions */
		rec[0] = rec[0]/(2.0*PI*Tconv);
		rcent = vecdiv(rcent,Rconv);
		vcent = vecdiv(vcent,Vconv);
		printf("%d %d %lg %lg %lg %lg %lg\n",j,first[1],rec[0],Rconv,Tconv,rec[3],rec[2]);
		//printf("%d %lg %lg %lg %lg %lg %lg\n",j,rcent.x,rcent.y,rcent.z,vcent.x,vcent.y,vcent.z);
	 }
       	fclose(infile);
	free(m);
	free(pos);
	free(vel);
	free(name);
       	return 0;
}

void zerovec(vector *r)
{
   (*r).x=0.0;
   (*r).y=0.0;
   (*r).z=0.0;
}

vector vecadd(vector a, vector b)
{
   vector sum;

   sum.x=a.x+b.x;
   sum.y=a.y+b.y;
   sum.z=a.z+b.z;
   return sum;
}

vector vecsub(vector a, vector b)
{
   vector sub;

   sub.x=a.x-b.x;
   sub.y=a.y-b.y;
   sub.z=a.z-b.z;
   return sub;
}

vector vecmult(float d, vector a)
{
   vector mult;

   mult.x=d*a.x;
   mult.y=d*a.y;
   mult.z=d*a.z;
   return mult;
}

vector vecdiv(vector a, float d)
{
   vector div;

   div.x=a.x/d;
   div.y=a.y/d;
   div.z=a.z/d;
   return div;
}

vector crossprod(vector a, vector b)
{
   vector c;

   c.x=a.y*b.z-a.z*b.y;
   c.y=a.z*b.x-a.x*b.z;
   c.z=a.x*b.y-a.y*b.x;
   return c;
}

float magnitude(vector a)
{
   float mag=0.0;

   mag=sqrt(a.x*a.x+a.y*a.y+a.z*a.z);
   return mag;
}
