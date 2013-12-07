#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <memory.h>

void calculatecumsum(unsigned char *,int, int, int**);
void boxfilter1(unsigned char *, int, int, int);


void guidedFilter(unsigned char * src, unsigned char * guidance,unsigned char * dest,int radius,float eps, int rows, int cols)
{
  int i,j, k, l;
  float mean_I, mean_Ip, cov_Ip, mean_II, var_I, a, b, mean_a, mean_b, mean_p, q;
  const float e=eps*eps; 
  float *distance1;
  float *distance2;
  float *distance3;
  float *distance4;
  float *distance5;
  float *distance6;
  float *N;
  int value = 1;
 //the size of each local patch; N=(2r+1)^2 except for boundary pixels
  unsigned char* unitfilter;
  unsigned char Ip;
  unsigned char* Ipp;
  unsigned char Ii;
  unsigned char* Iii;
  unsigned char ai;
  unsigned char* aii;
  unsigned char bi;
  unsigned char* bii;

  unitfilter = (unsigned char*) calloc(rows*cols, sizeof(unsigned char));
  
  for(k=0; k<rows; k++) 
	{
		for(int l=0; l<cols; l++) 
			{
				unitfilter[k*cols+l] = value; 
			}
	}

  
  
  boxfilter1(unitfilter, rows, cols, radius, N); //N = boxfilter(ones(hei, wid), r);
  boxfilter1(guidance, rows, cols, radius, distance1); 
  boxfilter1(src, rows, cols, radius, distance2); 
  
  Ip = (*src)*(*guidance);
  Ipp = &Ip;

  boxfilter1(Ipp, rows, cols, radius, distance3);
  

  mean_I = (*distance1)/(*N); 
  mean_p = (*distance2)/(*N);

  mean_Ip = (*distance3)/(*N);
  cov_Ip = mean_Ip - mean_I * mean_p;


  Ii = (*guidance)*(*guidance);
  Iii = &Ii;
  boxfilter1(Iii, rows, cols, radius, distance4);
  mean_II = (*distance4)/(*N);

  var_I = mean_II - mean_I * mean_I;
  a = cov_Ip / (var_I + eps);// Eqn. (5) in the paper;
  b = mean_p - a * mean_I;// Eqn. (6) in the paper;

  ai = (unsigned char)a;
  aii = &ai;

  bi = (unsigned char)b;
  bii = &bi;

  boxfilter1(aii, rows, cols, radius, distance5);
  boxfilter1(bii, rows, cols, radius, distance6);

  mean_a = (*distance5)/(*N);
  mean_b = (*distance6)/(*N);

  q = mean_a * (float)*src + mean_b;// Eqn. (8) in the paper;

}





void boxfilter1(unsigned char *image, int N, int M, int radius, float *distance)
{
  int m, n, i, j;
  int *cumsum; 
  //float *distance;

  (distance) = (float *) calloc(N*M, sizeof(float));
  calculatecumsum(image, M, N, &cumsum);
   
  
 // Cumulative SUM over Y axis 

   for(n = 0; n <N; n++)
      {
		  for(m = 0; m < M; m++)
         	{
				if(m==0)
				{
					distance[m,n]= distance[m,n];
				}
				else
				{
					distance[m,n]= distance[m,n] + cumsum[m-1, n];
				}
		  }
   }
         
   // Difference over Y axis 

   for(m = 0; m <=radius; m++)
      {
		  for(n = 0; n < N; n++)
			{
				distance[m,n]= cumsum[m+radius, n];
		  }
   }
   //imDst(1:r+1, :) = imCum(1+r:2*r+1, :);
      for(m = radius+1; m <= M-radius-1; m++)
      {
		  for(n = 0; n < N; n++)
			{
				distance[m,n]= cumsum[m+radius, n]-cumsum[m*radius-1, n];
		  }
   }
      
	  
	  for(m = M-radius; m < M; m++)
      {
		  for(n = 0; n < N; n++)
			{
				distance[m,n]= cumsum[M-1, n]-cumsum[m*radius-1, n];
		  }
   }					

// Cumulative SUM over X axis 
   for(m = 0; m <M; m++)
      {
		  for(n = 0; n < N; n++)
         	{
				if(n==0)
				{
					distance[m,n]= distance[m,n];
				}
				else
				{
					distance[m,n]= distance[m,n] + cumsum[m, n-1];
				}
		  }
   }


// Difference over X axis 

 for(m = 0; m <M; m++)
      {
		  for(n = 0; n <= radius; n++)
         	{
					distance[m,n]= cumsum[m,n+radius];
		  }	
   }
 free(distance);
}


void calculatecumsum(unsigned char *image, int rows, int cols, 
         int **cumsumim)
{
	int y,x;
	(*cumsumim) = (int *) calloc(rows*cols, sizeof(int));
	(*cumsumim) = (int *) image;
	for (y = 0 ; y < rows ; y++ )
	{
		for (x = 0; x < cols ; x++ )
			{
				(cumsumim)[y*cols+x] += (int)(cumsumim)[y*cols+ x];
			}
	}
}