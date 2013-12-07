#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <memory.h>

#define nelem(a) (sizeof(a) / (sizeof(a[0])))
#define min(x,y) ((x) < (y) ? (x) : (y))
#define max(x,y) ((x) < (y) ? (y) : (x))
//typedef double element;


//   2D MEAN FILTER implementation
//     image  - input image
//     result - output image
//     N      - width of the image (columns)
//     M      - height of the image (rows)
void meanfilter1(unsigned char *image, double* result, int N, int M)
{
  int m, n, i, j;
	//   Move window through all elements of the image
   for(m = 1; m < M - 1; m++)
      {
		  for(n = 1; n < N - 1; n++)
         //   Take the average
			{
				result[(m - 1) * (N - 2) + n - 1] = (image[(m - 1) * N + n - 1] + 
            image[(m - 1) * N + n] + 
            image[(m - 1) * N + n + 1] +
            image[m * N + n - 1] + 
            image[m * N + n] + 
            image[m * N + n + 1] +
            image[(m + 1) * N + n - 1] + 
            image[(m + 1) * N + n] + 
            image[(m + 1) * N + n + 1]) / 9;

			 //printf("%d\t", result);
		  }
   }
	//for (j=1; j<M-1; j++)   {
	 //for (i=1; i<N-1; i++){
		//printf("%d\t", (int) result[j*N+i]);
		//printf("hello");
	//}
	//printf("\n");
 //}
}
//   2D MEAN FILTER wrapper//
//     image  - input image
//     result - output image
//     N      - width of the image (columns)
//     M      - height of the image (rows)
void meanfilter(unsigned char* image, double* result, int N, int M)
{
   unsigned char *extension;
	int i;
	//   Check arguments
   if (!image || N < 1 || M < 1)
      return;
   //   Allocate memory for signal extension
   extension = (unsigned char *) calloc((N+2)*(M+2), sizeof(unsigned char));
   //double* extension = new double[(N + 2) * (M + 2)];
   //   Check memory allocation
   if (!extension)
      return;
   
   //   Create image extension
   for(i = 0; i < M; i++)
   {
      memcpy(extension + (N + 2) * (i + 1) + 1,image + N * i, N * sizeof(unsigned char));
      extension[(N + 2) * (i + 1)] = image[N * i];
      extension[(N + 2) * (i + 2) - 1] = image[N * (i + 1) - 1];
   }
   //   Fill first line of image extension
   memcpy(extension, extension + N + 2, (N + 2) * sizeof(unsigned char));
   //   Fill last line of image extension
   memcpy(extension + (N + 2) * (M + 1), extension + (N + 2) * M, (N + 2) * sizeof(unsigned char));
   //   Call mean filter implementation
   meanfilter1(extension, result ? result : (double *)image, N + 2, M + 2);
   //   Free memory
   free(extension);
}

void boxfilter1(unsigned char *image, double* distance, int N, int M, int radius)
{
  int m, n, i, j;
  int *cumsum; 
  calculatecumsum(image, M, N, &cumsum);
  distance[M,N] = 0;
  //   Move window through all elements of the image

  
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

 for(m = 0; m <M; m++)
      {
		  for(n = radius+1; n <= N-radius-1; n++)
         	{
					distance[m,n]= cumsum[m,n+radius]-cumsum[m,n*radius-1];
		  }	
   }
  
 for(m = 0; m <M; m++)
      {
		  for(n = radius+1; n <= N-radius-1; n++)
         	{
					distance[m,n]= cumsum[m,n+radius]-cumsum[m,n*radius-1];
		  }	
   }
 
 for(m = 0; m <M; m++)
      {
		  for(n = N-radius; n < N; n++)
         	{
					distance[m,n]= cumsum[m,N-1]-cumsum[m,n*radius-1];
		  }	
   }
 //return distance;
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


void guidedfilter(unsigned char *srcimg, unsigned char *guideimg, float *result, int radius, float eps)
{
  int i,j;
  float N,mean_I, mean_Ip, cov_Ip, mean_II, var_I, a, b, mean_a, mean_b;
   
 //the size of each local patch; N=(2r+1)^2 except for boundary pixels
[hei, wid] = size(I);
  center = (*windowsize) / 2;
N = boxfilter(ones(hei, wid), r); //the size of each local patch; N=(2r+1)^2 except for boundary pixels.

mean_I = boxfilter(I, r) ./ N;
mean_p = boxfilter(p, r) ./ N;
mean_Ip = boxfilter(I.*p, r) ./ N;
cov_Ip = mean_Ip - mean_I .* mean_p; //this is the covariance of (I, p) in each local patch.

mean_II = boxfilter(I.*I, r) ./ N;
var_I = mean_II - mean_I .* mean_I;


// calculating the linear coefficients a and b
a = cov_Ip ./ (var_I + eps); // Eqn. (5) in the paper;
b = mean_p - a .* mean_I; // Eqn. (6) in the paper;

mean_a = boxfilter(a, r) ./ N;
mean_b = boxfilter(b, r) ./ N;

q = mean_a .* I + mean_b; // Eqn. (8) in the paper;



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
