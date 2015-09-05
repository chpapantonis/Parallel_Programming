#include <stdio.h>
#include <cuda.h>

#define Width 1920
#define Height 2520
#define iterations 100


__global__ void convolution_kernel(unsigned char* in_device_buffer,unsigned char* out_device_buffer);
void swap_images(unsigned char **in_image,unsigned char **out_image);


int main (){
	
	unsigned char *in_image ;
	unsigned char *in_device_buffer;
	unsigned char *out_image ;
	unsigned char *out_device_buffer ;
	int grid, i;
	grid = Width* Height;
	in_image =(unsigned char*)malloc(grid);
	out_image =(unsigned char*)malloc(grid);
	
	
	FILE *fp;
	fp = fopen("images/grey_X1.raw","rb");
	if (fp == NULL){
		printf("Empty file ... Exiting");
		exit(1);
	}else {
		fread(in_image,grid,1,fp);
		fclose(fp);
	}
	 
	
	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	
	dim3 threadNum(32,32);
	dim3 blockNum(Width/((threadNum.x)-2),Height/((threadNum.y)-2));
	
	
	
	
	
	
	cudaMalloc(&in_device_buffer,grid);
	cudaMalloc(&out_device_buffer,grid);
	cudaMemcpy(in_device_buffer,in_image,grid,cudaMemcpyHostToDevice);
	
	cudaEventRecord(start, 0);
	for (i = 0 ; i < iterations; i++){
		convolution_kernel <<<blockNum,threadNum>>> (in_device_buffer, out_device_buffer);
		swap_images(&in_device_buffer,&out_device_buffer);	
	}
	cudaEventRecord(stop, 0);


	cudaMemcpy(out_image, out_device_buffer, grid, cudaMemcpyDeviceToHost);
	cudaFree(out_device_buffer);
	cudaFree(in_device_buffer);
	
	
	cudaEventSynchronize(stop);
	
	fp = fopen("out.raw","w+");
	fwrite(out_image,grid,1,fp);
	fclose(fp);
	
return 0;

}


void swap_images(unsigned char **in_image,unsigned char **out_image)
{
	unsigned char* temp = *in_image;
	*in_image = *out_image;
	*out_image = temp;
}

__global__ void convolution_kernel(unsigned char* A, unsigned char* B)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	int j = blockIdx.y * blockDim.y + threadIdx.y;
	
	int x = i-2*blockIdx.x-1;
	int y = j-2*blockIdx.y-1;
	
	__shared__ unsigned char As[32][32];
	
	//Copy from global memory to shared memory	
		
		if (x<0) {
			x=0;
		} else if (x==Width) {
			x=Width-1;
		}
		if (y<0) {
			y=0;
		} else if (y == Height) {
			y = Height-1;
		}
		As[threadIdx.x][threadIdx.y] = A[Width*y + x];
		
		__syncthreads();
	
	// Computations
	
		if (threadIdx.x!=0 && threadIdx.x!=31 && threadIdx.y!=0 && threadIdx.y!=31) {
			B[Width*y + x] =     (As[threadIdx.x-1][threadIdx.y-1]  +
										As[threadIdx.x  ][threadIdx.y-1] * 2 +
										As[threadIdx.x+1][threadIdx.y-1]  +
										As[threadIdx.x-1][threadIdx.y  ] *2 +
										As[threadIdx.x  ][threadIdx.y  ] *4 +
										As[threadIdx.x+1][threadIdx.y  ] * 2 +
										As[threadIdx.x-1][threadIdx.y+1] * 1 +
										As[threadIdx.x  ][threadIdx.y+1] * 2 +
										As[threadIdx.x+1][threadIdx.y+1] * 1)/16;
		}
}
