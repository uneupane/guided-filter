/**************************************
 * Zeyun Yu (yuz@uwm.edu)             *
 * Department of Computer Science     *
 * University of Wisconsin-Milwaukee  *
 **************************************/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <malloc.h>  
#include <memory.h>

#define BOOSTBLURFACTOR 90.0
#define max(x, y) ((x>y) ? (x):(y))
#define min(x, y) ((x<y) ? (x):(y))
#define SIGMA         2


void gaussian_smooth(unsigned char *, int , int , short int **);
void make_gaussian_kernel(float , float **, int *);
void derrivative_x_y(short int *, int , int , short int **, short int **);
void magnitude_x_y(short int *, short int *, int , int , float *);


void EdgeMag(int rows, int cols, unsigned char *img, float *edge_mag)
{
  int i,j;
  float max_mag,min_mag;
  short int *delta_x, *delta_y;
  short int *smoothedim;
  
  gaussian_smooth(img, rows, cols, &smoothedim);
  derrivative_x_y(smoothedim, rows, cols, &delta_x, &delta_y);
  free(smoothedim);
  magnitude_x_y(delta_x, delta_y, rows, cols, edge_mag);
  free(delta_x);
  free(delta_y);
  
  max_mag = -99999;
  min_mag = 99999;
  for (j=0; j<rows; j++)   
    for (i=0; i<cols; i++) {
      if (edge_mag[j*cols+i] > max_mag)
	max_mag = edge_mag[j*cols+i];
      if (edge_mag[j*cols+i] < min_mag)
	min_mag = edge_mag[j*cols+i];
    }


  for (j=0; j<rows; j++)   
    for (i=0; i<cols; i++)
      edge_mag[j*cols+i] = 255*(edge_mag[j*cols+i]-min_mag)/
	                       (max_mag-min_mag);
 
}


void magnitude_x_y(short int *delta_x, short int *delta_y, int rows, int cols,
        float *magnitude)
{
  int r, c, pos, sq1, sq2;

   
  for(r=0,pos=0;r<rows;r++){
    for(c=0;c<cols;c++,pos++){
      sq1 = (int)delta_x[pos] * (int)delta_x[pos];
      sq2 = (int)delta_y[pos] * (int)delta_y[pos];
      magnitude[pos] = 0.5 + sqrt((float)sq1 + (float)sq2);
    }
  }

}



void derrivative_x_y(short int *smoothedim, int rows, int cols,
        short int **delta_x, short int **delta_y)
{
   int r, c, pos;

   (*delta_x) = (short *) calloc(rows*cols, sizeof(short));
   (*delta_y) = (short *) calloc(rows*cols, sizeof(short));
   
   for(r=0;r<rows;r++){
      pos = r * cols;
      (*delta_x)[pos] = smoothedim[pos+1] - smoothedim[pos];
      pos++;
      for(c=1;c<(cols-1);c++,pos++){
         (*delta_x)[pos] = smoothedim[pos+1] - smoothedim[pos-1];
      }
      (*delta_x)[pos] = smoothedim[pos] - smoothedim[pos-1];
   }

   for(c=0;c<cols;c++){
      pos = c;
      (*delta_y)[pos] = smoothedim[pos+cols] - smoothedim[pos];
      pos += cols;
      for(r=1;r<(rows-1);r++,pos+=cols){
         (*delta_y)[pos] = smoothedim[pos+cols] - smoothedim[pos-cols];
      }
      (*delta_y)[pos] = smoothedim[pos] - smoothedim[pos-cols];
   }
}


void gaussian_smooth(unsigned char *image, int rows, int cols, 
        short int **smoothedim)
{
   int r, c, rr, cc,     /* Counter variables. */
      windowsize,        /* Dimension of the gaussian kernel. */
      center;            /* Half of the windowsize. */
   float *tempim,        /* Buffer for separable filter gaussian smoothing. */
         *kernel,        /* A one dimensional gaussian kernel. */
         dot,            /* Dot product summing variable. */
         sum;            /* Sum of the kernel weights variable. */
   float sigma;


   sigma = SIGMA;

   make_gaussian_kernel(sigma, &kernel, &windowsize);
   center = windowsize / 2;

   tempim = (float *) calloc(rows*cols, sizeof(float));

   (*smoothedim) = (short int *) calloc(rows*cols, sizeof(short int));
  
   for(r=0;r<rows;r++){
      for(c=0;c<cols;c++){
         dot = 0.0;
         sum = 0.0;
         for(cc=(-center);cc<=center;cc++){
            if(((c+cc) >= 0) && ((c+cc) < cols)){
               dot += (float)image[r*cols+(c+cc)] * kernel[center+cc];
               sum += kernel[center+cc];
            }
         }
         tempim[r*cols+c] = dot/sum;
      }
   }

   for(c=0;c<cols;c++){
      for(r=0;r<rows;r++){
         sum = 0.0;
         dot = 0.0;
         for(rr=(-center);rr<=center;rr++){
            if(((r+rr) >= 0) && ((r+rr) < rows)){
               dot += tempim[(r+rr)*cols+c] * kernel[center+rr];
               sum += kernel[center+rr];
            }
         }
         (*smoothedim)[r*cols+c] = (short int)(dot*BOOSTBLURFACTOR/sum + 0.5);
      }
   }

   free(tempim);
   free(kernel);
}


void make_gaussian_kernel(float sigma, float **kernel, int *windowsize)
{
   int i, center;
   float x, fx, sum=0.0;

   *windowsize = 1 + 2 * ceil(2.5 * sigma);
   center = (*windowsize) / 2;

   (*kernel) = (float *) calloc((*windowsize), sizeof(float));

   for(i=0;i<(*windowsize);i++){
     x = (float)(i - center);
     fx = pow(2.71828, -0.5*x*x/(sigma*sigma)) / (sigma * sqrt(6.2831853));
     (*kernel)[i] = fx;
     sum += fx;
   }

   for(i=0;i<(*windowsize);i++) (*kernel)[i] /= sum;

   
}
