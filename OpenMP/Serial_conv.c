//========================================================
//				SERIAL CONVERGENCE OF THE GRAYSCALE AND RGB IMAGE
//	 	
//Developers: Papantonis.X , Tzimas.G
//Supervisor: Cotronis.I
//
//=========================================================
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

//image attributes
#define image_width 1920
#define image_height 2520
#define image_colour 1

//I/O filenames
#define i_image "initial_image.raw"
#define o_image "conv_image.raw"



void main (){
	
	const unsigned char filter_mask[9];
	int i , j, pass; 
	clock_t clk_start, clk_stop; // variable fto time
	FILE *filePointer //image pointer file 
	
	
	short **image;
	
	for(pass = 0 ; pass <= 300; pass++){
		array = (short*)malloc(image_height*sizeof(short*));//Dynamic memory allocation 
		for (i = 0; i < image_height ; i++){
			*(array +i) = (short)malloc(image_width*sizeof(short));
			for (j = 0; i < image_width ; j++){
					for 
			
			}
			
			
			
		}
	
	
		
		
		
		
	}
	
	
	
	
	
}
