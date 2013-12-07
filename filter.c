#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <memory.h>

void boxfilter(float *, int, int, int, float *);
void unitfilter(int, int, float *);
float* matDiff(float* one, float* other, int rows, int cols);
float* repmat(float* image, int rows, int cols, int repRows, int repCols);
void printMat(float* mat, int rows, int cols);
float* mean(float* value, float* total,int rows, int cols);
float* imagemultiplication(float* one, float* other, int rows, int cols);
float* calculatevariance(float* mean1, float* mean2,float* mean3,int rows, int cols);
void calculatea(float* value1, float* value2, float eps, int rows, int cols, float* result);
void calculateoutput(float* meancoeffa, float* image, float* meancoeffb, int rows, int cols, float* result);
void printMat2(unsigned char* mat, int rows, int cols);

void guidedFilter(float * guidance, float * src,float * dest,int radius,float eps, int rows, int cols)
{
  int i,j, k, l;
  float *mean_I, *mean_p, *mean_Ip, *cov_Ip, *mean_II, *var_I; 
  float *distance; // distance of local patch
  float *guidanceimdistance; // distance of guidance image
  float *filterimdistance; // distance of guidance image
  float *Ipimdistance; // distance of I*p image
  float *IIimdistance; // distance of I*I
  float *coeffadistance; // distance of coeffa
  float *coeffbdistance; // distance of coeffb
  float *coeffa, *coeffb, *mean_coeffa, *mean_coeffb, *output;
  

  //unsigned char* unf;
  float* unf;
  float* improduct;
  float* imsquareproduct;
  
  //unf = (unsigned char*)malloc(sizeof(unsigned char)*rows*cols);
  unf = (float*)malloc(sizeof(float)*rows*cols);
  unitfilter(rows, cols, unf);
  //printf("\Printing unit filter\n");
  //printMat(unf,rows,cols);

  /*N = boxfilter(ones(hei, wid), r);*/
  distance = (float*)malloc(sizeof(float)*rows*cols);
  boxfilter(unf, rows, cols, radius, distance);
  //printf("\Printing patch distance\n");
  //printMat(distance,rows,cols);

  /* Mean of guidance image: mean_I = boxfilter(I, r) ./ N;*/
  guidanceimdistance = (float*)malloc(sizeof(float)*rows*cols);
  boxfilter(guidance, rows, cols, radius, guidanceimdistance);	
  //printf("\Printing guidanceimdistance\n");
  //printMat(guidanceimdistance,rows,cols);

  mean_I = (float*)malloc(sizeof(float)*rows*cols);
  mean_I = mean(guidanceimdistance,distance,rows,cols);
  //printf("\nPrinting mean\n");
  //printMat(mean_I,rows,cols);

  /* Mean of filtering input image: mean_p = boxfilter(p, r) ./ N;*/
  filterimdistance = (float*)malloc(sizeof(float)*rows*cols);
  boxfilter(src, rows, cols, radius, filterimdistance);	
  //printf("\Printing filterimdistance\n");
  //printMat(filterimdistance,rows,cols);
  
  mean_p = (float*)malloc(sizeof(float)*rows*cols);
  mean_p = mean(filterimdistance,distance,rows,cols);
  //printf("\nPrinting mean p\n");
  //printMat(mean_p,rows,cols);
  
  /* Mean of I*p: mean_Ip = boxfilter(I.*p, r) ./ N;*/
  improduct = (float*)malloc(sizeof(float)*rows*cols);
  improduct = imagemultiplication(guidance,guidance,rows,cols); 
  //printf("\nPrinting Ip\n");
  //printMat(improduct,rows,cols);  

  Ipimdistance = (float*)malloc(sizeof(float)*rows*cols);
  boxfilter(improduct, rows, cols, radius, Ipimdistance);	
  //printf("\nPrinting Ipimdistance\n");
  //printMat(Ipimdistance,rows,cols);  

  mean_Ip = (float*)malloc(sizeof(float)*rows*cols);
  mean_Ip = mean(Ipimdistance,distance,rows,cols);
  //printf("\nPrinting mean\n");
  //printMat(mean_Ip,rows,cols);
  
  //cov_Ip = mean_Ip - mean_I .* mean_p;% this is the covariance of (I, p) in each local 
  cov_Ip = (float*)malloc(sizeof(float)*rows*cols);
  cov_Ip = calculatevariance(mean_Ip,mean_I,mean_p,rows,cols);
  //printf("\nPrinting cov_Ip\n");
  //printMat(cov_Ip,rows,cols);

  /* mean_II = boxfilter(I.*I, r) ./ N;*/
  imsquareproduct = (float*)malloc(sizeof(float)*rows*cols);
  imsquareproduct = imagemultiplication(guidance,guidance,rows,cols);
  //printf("\nPrinting I2\n");
  //printMat(imsquareproduct,rows,cols); 
  
  IIimdistance = (float*)malloc(sizeof(float)*rows*cols);
  boxfilter(imsquareproduct, rows, cols, radius, IIimdistance);
  
  //printf("\nPrinting IIimdistance\n");
  //printMat(IIimdistance,rows,cols);
  
  mean_II = (float*)malloc(sizeof(float)*rows*cols);
  mean_II = mean(IIimdistance,distance,rows,cols);
  //printf("\nPrinting mean\n");
  //printMat(mean_II,rows,cols);
  
  /*var_I = mean_II - mean_I .* mean_I;*/
  var_I = (float*)malloc(sizeof(float)*rows*cols);
  var_I = calculatevariance(mean_II,mean_I,mean_I,rows,cols);
  //printf("\nPrinting var_I\n");
  //printMat(var_I,rows,cols);

  /*a = cov_Ip ./ (var_I + eps); % Eqn. (5) in the paper;*/
  coeffa = (float*)malloc(sizeof(float)*rows*cols);
  calculatea(cov_Ip,var_I,eps,rows,cols,coeffa);

  /*b = mean_p - a .* mean_I; % Eqn. (6) in the paper;*/
  coeffb = (float*)malloc(sizeof(float)*rows*cols);
  coeffb = calculatevariance(mean_p,coeffa,mean_I,rows,cols);

  /*mean_a = boxfilter(a, r) ./ N;*/

  coeffadistance = (float*)malloc(sizeof(float)*rows*cols);
  boxfilter(coeffa, rows, cols, radius, coeffadistance);	
  
  mean_coeffa = (float*)malloc(sizeof(float)*rows*cols);
  mean_coeffa = mean(coeffadistance,distance,rows,cols);

  
  /*mean_b = boxfilter(b, r) ./ N;*/
  coeffbdistance = (float*)malloc(sizeof(float)*rows*cols);
  boxfilter(coeffb, rows, cols, radius, coeffbdistance);	
  
  mean_coeffb = (float*)malloc(sizeof(float)*rows*cols);
  mean_coeffb = mean(coeffbdistance,distance,rows,cols);

  /*q = mean_a .* I + mean_b; % Eqn. (8) in the paper;*/
  //output = (float*)malloc(sizeof(float)*rows*cols);
  calculateoutput(mean_coeffa,guidance,mean_coeffb,rows,cols,dest);


   
  free(unf);
  free(distance);
  free(guidanceimdistance);
  free(mean_I);
  free(filterimdistance);
  free(mean_p);
  free(improduct);
  free(Ipimdistance);
  free(mean_Ip);
  free(cov_Ip);
  free(imsquareproduct);
  free(IIimdistance);
  free(mean_II);
  free(var_I);
  free(coeffa);
  free(coeffb);
  free(coeffadistance);
  free(coeffbdistance);
  free(mean_coeffa);
  free(mean_coeffb);

}



void unitfilter(int rows, int cols, float *unf)
{
  int k,l;
  float value = 1.0;

  for(k=0; k<rows; k++) 
	{
		for(l=0; l<cols; l++) 
			{
				unf[k*cols+l] = value; 
			}
	}
}


void boxfilter(float *image,  int M, int N, int radius, float *distance)
{
  int m, n, i, j;
  float *cumdistance;
  float *cumheightvector; //1-D, After repmat will become 2-D
  float *cumheightvectorcol;
  float *cumheightvector2d;
  float *cumheightvector2dcol;
  float *cumdistancevector;
  float *cumdistancevectorcol;
  float *cumheightdiffdistance; //cumheightvector2d - cumheightvector
  float *cumheightdiffdistancecol;
  float *printcumheight;
  float *chvc;
 cumdistance = (float*)malloc(sizeof(float)*M*N);
  // Cumulative SUM over Y axis 
   for(n = 0; n <N; n++)
      {
		  for(m = 0; m < M; m++)
         	{
				if(m==0)
				{
					cumdistance[m*N+ n]= (float) image[m*N + n];
				}
				else
				{
					cumdistance[m*N+ n] = (float) image[m*N + n] + (float) cumdistance[(m-1)*N+n];
				}
		  }
   }
         
 
   // Difference over Y axis 

   for(m = 0; m <=radius; m++)
      {
		  for(n = 0; n < N; n++)
			{
				distance[m*N+ n] = (float) cumdistance[(m+radius)*N+ n];
		  }
   }

   for(m = radius+1; m <= M-(radius+1); m++)
      {
		  for(n = 0; n < N; n++)
			{
				//imDst(r+2:hei-r, :) = imCum(2*r+2:hei, :) - imCum(1:hei-2*r-1, :);
				distance[m*N+ n] = (float) cumdistance[(m+radius)*N + n]-(float) cumdistance[(m-radius-1)*N+ n];
		  }
   }
 
 //imDst(hei-r+1:hei, :) = repmat(imCum(hei, :), [r, 1]) - imCum(hei-2*r:hei-r-1, :);
  cumheightvector = (float*) calloc(N,sizeof(float*));
   for(n=0; n<N; n++) {
	   cumheightvector[n] = cumdistance[(M-1)*N + n];
   }
  
   cumdistancevector = (float*) calloc(N*M,sizeof(float*));
   for (i=0,m=M-1-2*radius; m<=M-1-radius-1; m++,i++) {
	   for (j=0,n=0; n<N; n++,j++) {
			cumdistancevector[i*N+j] = cumdistance[m*N+n];
	   }
   }

   cumheightvector2d = repmat(cumheightvector,1,N,radius,1);
   cumheightdiffdistance = matDiff(cumheightvector2d, cumdistancevector, radius, N);
   //printf("\nPrinting Cumheightdiffdistance\n");
   //printMat(cumheightdiffdistance,M,N);
  for(i=0,m = M-radius; m < M; m++,i++)
      {
		  for(j=0,n = 0; n < N; n++,j++)
			{
				distance[m*N+n] = cumheightdiffdistance[i*N+j];
		  }
   }					

  //Cumulative SUM over X axis 
	for(m = 0; m <M; m++)
      {
		  for(n = 0; n < N; n++)
         	{
				if(n==0)
				{
					cumdistance[m*N+ n]= (float) distance[m*N + n];
				}
				else
				{
					cumdistance[m*N+ n] = (float) distance[m*N+ n] + cumdistance[m*N+(n-1)];
				}
		  }
   }

// Difference over X axis 
//imDst(:, 1:r+1) = imCum(:, 1+r:2*r+1)
 for(m = 0; m <M; m++)
      {
		  for(n = 0; n <= radius; n++)
         	{
					distance[m*N+ n] = (float) cumdistance[m*N+(n+radius)];
		  }	
	}

 //imDst(:, r+2:wid-r) = imCum(:, 2*r+2:wid) - imCum(:, 1:wid-2*r-1);

 for(m = 0; m <M; m++)
      {
		  for(n = radius+1; n <= N-radius-1; n++)
         	{
					distance[m*N+ n] = (float) cumdistance[m*N+(n+radius)]-(float) cumdistance[m*N+(n-radius-1)];
		  }	
   }

//imDst(:, wid-r+1:wid) = repmat(imCum(:, wid), [1, r]) - imCum(:, wid-2*r:wid-r-1);

  cumheightvectorcol = (float*) calloc(M,sizeof(float*));
  for(m=0; m<M; m++) {
	  	   cumheightvectorcol[m] = cumdistance[m*N + N-1];
  }
  //printf("\n Printing cumheightvectorcol \n");
  //printMat(cumheightvectorcol,M,1);

 cumdistancevectorcol = (float*)malloc(sizeof(float)*M*N);
 cumdistancevector = (float*) calloc(radius*M,sizeof(float*));
   for (i=0,m=0; m<M; m++,i++) {
	   for (j=0,n=N-1-2*radius; n<=N-1-radius-1; n++,j++) {
					cumdistancevectorcol[i*radius+j] = (float) cumdistance[m*N+n];
		  }	
   }

   cumheightvector2dcol = repmat(cumheightvectorcol,M,1,1,radius); //M*radius
   cumheightdiffdistancecol = matDiff(cumheightvector2dcol, cumdistancevectorcol, M, radius);

   for(i=0,m = 0; m <M; m++,i++)
      {
		  for(j=0,n = N-radius; n <N; n++,j++)
         	{
				distance[m*N+n] = cumheightdiffdistancecol[i*radius+j];	
		  }	
   }

 free(cumdistance);
 free(cumheightvector);
 free(cumheightvectorcol);
 free(cumheightvector2d);
 free(cumheightvector2dcol);
 free(cumdistancevector);
 free(cumdistancevectorcol);
 free(cumheightdiffdistance);
 free(cumheightdiffdistancecol);
}

//								M		  1			1			 radius=2
 float* repmat(float* image, int rows, int cols, int repRows, int repCols) 
  {
	int totRows = rows * repRows;//M
	int totCols = cols * repCols;//2
	int r,c, ro, co;
	int j,i;
	float val;
	float* newImage = (float*) calloc(totRows * totCols, sizeof(float*));
	//Loop thru rows and cols of given image
	for (r=0; r<rows; r++) {
		for (c=0; c<cols; c++) {
			val = image[r*cols+c];
			//Loop thru rows and cols of the new image and replicate the value
			for (ro=r; ro<totRows; ro=ro+rows){
				for (co=c; co<totCols; co=co+cols){
					newImage[ro*totCols+co] = val;
				}
			}
		}
	}
	return newImage;
 }

float* imagemultiplication(float* one, float* other, int rows, int cols) {
	int i,j,m,n;
	float* B = (float*)malloc(sizeof(float)*rows*cols);
	for (j=0; j<rows; j++) {  
		for (i=0; i<cols; i++){
			B[j*cols+i] = one[j*cols+i]*other[j*cols+i]; 
		}
	}
	return B;
}

float* matDiff(float* one, float* other, int rows, int cols) {
	int i,j;
	float* result = (float*) calloc (rows*cols, sizeof(float*));
	for (j=0; j<rows; j++) {  
		for (i=0; i<cols; i++){
			result[j*cols+i] = one[j*cols+i] - other[j*cols+i];
		}
	}
	return result;
}

float* mean(float* value, float* total,int rows, int cols) {
	int i,j;
	float sumtotal;
	float* result = (float*) calloc (rows*cols, sizeof(float*));
	for (j=0; j<rows; j++) {  
		for (i=0; i<cols; i++){
			result[j*cols+i] = value[j*cols+i] /total[j*cols+i];
		}
	}
	return result;
}

float* calculatevariance(float* mean1, float* mean2,float* mean3,int rows, int cols) {
	int i,j;
	float sumtotal;
	float* result = (float*) calloc (rows*cols, sizeof(float*));
	for (j=0; j<rows; j++) {  
		for (i=0; i<cols; i++){
			result[j*cols+i] = mean1[j*cols+i] - mean2[j*cols+i]*mean3[j*cols+i];
		}
	}
	return result;
}

void calculatea(float* value1, float* value2, float eps, int rows, int cols, float* result){
	int i,j;
	for (j=0; j<rows; j++) {  
		for (i=0; i<cols; i++){
			result[j*cols+i] = value1[j*cols+i]/(value2[j*cols+i]+eps);
		}
	}
}

void calculateoutput(float* meancoeffa, float* image, float* meancoeffb, int rows, int cols, float* result){
	int i,j;
	for (j=0; j<rows; j++) {  
		for (i=0; i<cols; i++){
			result[j*cols+i] = meancoeffa[j*cols+i]*(float)image[j*cols+i]+meancoeffb[j*cols+i];
		}
	}
}


void printMat2(unsigned char* mat, int rows, int cols) {
	int i,j;
	for (j=0; j<rows; j++) {  
		for (i=0; i<cols; i++){
			printf("%5d  ", (int) mat[j*cols+i]); 
		}
		printf("\n"); 
	}
}
