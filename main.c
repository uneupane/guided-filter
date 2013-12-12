
/* Guided Image Filtering by by Kaiming He, Jian Sun, and Xiaoou Tang, in ECCV 2010 (Oral)

%   GUIDEDFILTER   O(1) time implementation of guided filter
%
%   - guidance image: I (should be a gray-scale/single channel image)
%   - filtering input image: p (should be a gray-scale/single channel image)
%   - local window radius: r
%   - regularization parameter: eps

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


void ReadPGM(FILE* , unsigned char** , int* , int* );
void WritePGM(int, int, unsigned char *, unsigned char *, FILE*);
void printMat2(unsigned char* mat, int rows, int cols);
void printMat(float* mat, int rows, int cols);
void guidedFilter(float* guidance, float * src,float * dest,int radius,float eps, int rows, int cols);

main(int argc, char *argv[])
{
  int rows, cols, radius;
  float eps;
  unsigned char *image1;
  unsigned char *image2;
  unsigned char *output;

  float *fimage1;
  float *fimage2;
  float *distance;
  float *filtered_image;

  FILE * fp;
  int i,j;

  if (argc != 6){
    printf("Usage: Test <guidance_image_filename><input_filtering_image_filename><radius><eps><output_filename> \n");
    printf("       <guidance_image_filename>:   PGM file \n");
	printf("       <input_filtering_image_filename>:   PGM file \n");
    printf("       <local window radius>:   r\n");
	printf("       <regularization parameter>:  eps\n");
	printf("       <output_filename>:  PGM file \n");

    exit(0);              
  }


  printf("Reading PGM.... \n");
  if ((fp=fopen(argv[1], "rb"))==NULL){
    printf("reading error...\n");
    exit(0); 
  }
  ReadPGM(fp,&image1,&rows,&cols);
  //printMat2(image1,rows,cols);
  
  if ((fp=fopen(argv[2], "rb"))==NULL){
    printf("reading error...\n");
    exit(0); 
  }
  ReadPGM(fp,&image2,&rows,&cols);

  if (argv[3]==NULL){
    printf("Please enter the local window radius...\n");
    exit(0); 
  }
  else
  radius = atoi(argv[3]);
  printf("\n Radius = %d\n", radius);
  eps = atof(argv[4]);
  printf("\n eps = %f\n", eps);

  /* you may replace your own applications here */

  fimage1 = (float*)malloc(sizeof(float)*rows*cols);
  for (j=0; j<rows; j++) {  
    for (i=0; i<cols; i++){
		fimage1[j*cols+i] = (float)image1[j*cols+i]/255;
	}
  }

  fimage2 = (float*)malloc(sizeof(float)*rows*cols);
  for (j=0; j<rows; j++) {  
    for (i=0; i<cols; i++){
		fimage2[j*cols+i] = (float)image2[j*cols+i]/255;
	}
  }
  filtered_image = (float*)malloc(sizeof(float)*rows*cols);
  printf("begin calculting filtered_output.... \n");
  guidedFilter(fimage1, fimage2, filtered_image, radius, eps, rows, cols);
  //printf("\nPrinting result\n");
  //printMat(filtered_image,rows,cols);
  /* end of your application */

  output = (unsigned char*)malloc(sizeof(unsigned char)*rows*cols);
  for (j=0; j<rows; j++){   
    for (i=0; i<cols; i++){
		output[j*cols+i] = (unsigned char)255*filtered_image[j*cols+i];
	}
  }
  //printf("\nPrinting result\n");
  //printMat2(image,rows,cols);

  //printf("\nPrinting result\n");
  //printMat2(output,rows,cols);
  if ((fp=fopen(argv[5], "wb"))==NULL){
    printf("reading error...\n");
    exit(0); 
  }
  printf("Writing PGM....\n");
  WritePGM(rows, cols, image2, output, fp);

 printf("\nFreeing memory\n");
 free(filtered_image);
 free(output);
}


void printMat(float* mat, int rows, int cols) {
	int i,j;
	for (j=0; j<rows; j++) {  
		for (i=0; i<cols; i++){
			printf("%5f  ",  mat[j*cols+i]); 
		}
		printf("\n"); 
	}
}




