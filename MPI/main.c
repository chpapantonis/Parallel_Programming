#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <string.h>

#define Height 2520
#define Width 1920



//filter mask
//filter mask
const float h[3][3] = {{0.0625,0.125,0.0625}, {0.125 ,0.25,0.125},{0.0625,0.125,0.0625}} ;

//Setting the height and width of the sub-images global
int Ph, Pw;
//number of dimensions (grey or rgb)
int d;
//functions
void filteringInnerImage (float *input_imag , float *output_imag);
void filteringOutterImage (float *input_imag , float *output_imag );



int main(int argc, char** argv) {


    MPI_Comm cart_comm;
    MPI_Datatype Filetype;
    MPI_File* fp;

    // total number of processes
    int world_size;
    // process rank in MPI_COMM_WORLD
    int world_rank;
    // number of processors in each dimension
    int dims;
    int proc_per_dim[2];

    // specify periodicity in dimensions
    int periods[2] = {0, 0};
    // Cartesian coordinates of process
    int coords[2];
    int temp_coords[2];

    // ranks of source and destination processes for local communication
    int up_rank, down_rank, left_rank, right_rank;
    int up_left_rank, up_right_rank, down_left_rank, down_right_rank;

    //
    double time0, time1, time2, time3;
    double comp_time, max_comp_time;
    double tot_time, max_tot_time;

    // number of iteration
    int Iter;
    // total number of iterations
    int TotIter;



    int file_offset;

    //
    unsigned char* image;

    // data before filtering
    float *input_imag;
    // filtered data
    float *output_imag;
    // change pointers to data to avoid copy of matrices
    float *temp_pointer;


    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    int a[world_size];

    //getting the #of pixels on each subimage acording to processors we use
    Ph = Height/sqrt(world_size);
    Pw = Width/sqrt(world_size);
    dims = sqrt(world_size);
    proc_per_dim[0] = dims;
    proc_per_dim[1] = dims;

    /* CREATE CASRTESIAN TOPOLOGY */
    MPI_Cart_create(MPI_COMM_WORLD, 2, proc_per_dim, periods, 1, &cart_comm);
    MPI_Cart_coords(cart_comm, world_rank, 2, coords);
    MPI_Cart_shift(cart_comm, 0, 1, &up_rank, &down_rank);
    MPI_Cart_shift(cart_comm, 1, 1, &left_rank, &right_rank);


    //Setting the corner pixels communication
    //upper left
    if (coords[0] != 0 && coords[1] != 0){
        temp_coords[0] = coords[0]-1;
        temp_coords[1] = coords[1]-1;
        MPI_Cart_rank(cart_comm, temp_coords, &up_left_rank);
    }
    else{
        up_left_rank = -1;
    }
    //upper right
    if (coords[0] != 0 && coords[1] != dims-1){
        temp_coords[0] = coords[0]-1;
        temp_coords[1] = coords[1]+1;
        MPI_Cart_rank(cart_comm, temp_coords, &up_right_rank);
    }
    else{
        up_right_rank = -1;
    }
    //bottom left
    if (coords[0] != dims-1 && coords[1] != 0){
        temp_coords[0] = coords[0]+1;
        temp_coords[1] = coords[1]-1;
        MPI_Cart_rank(cart_comm, temp_coords, &down_left_rank);
    }
    else{
        down_left_rank = -1;
    }
    //bottom right
    if (coords[0] != dims-1 && coords[1] != dims-1){
        temp_coords[0] = coords[0]+1;
        temp_coords[1] = coords[1]+1;
        MPI_Cart_rank(cart_comm, temp_coords, &down_right_rank);
    }
    else{
        down_right_rank = -1;
    }


    if (world_rank==0) {
        printf("number of argc:%d\n", argc);
        int i;
        for (i=0; i<argc; i++) {
            printf("argv[%d]: %s\n", i, argv[i]);
        }
    }
    
    TotIter = atoi(argv[2]); // number of iteration from arguments
    
    
	 fp = malloc(sizeof(MPI_File));
    if (atoi(argv[1])==1) {
        d = 1;
        MPI_File_open(cart_comm, "images/grey_X1.raw", MPI_MODE_RDONLY, MPI_INFO_NULL, fp);
    } else if (atoi(argv[1])==3) {
		MPI_File_open(cart_comm, "images/waterfall_1920_2520.raw", MPI_MODE_RDONLY, MPI_INFO_NULL, fp);
        d = 3;
    }
    

    time0 = MPI_Wtime();
    image = malloc(d*Pw*Ph);
    MPI_Type_vector(Ph, d*Pw, d*Width, MPI_UNSIGNED_CHAR, &Filetype);
    MPI_Type_commit(&Filetype);

    file_offset = d*Pw*coords[1] + d*Width*Ph*coords[0];

    /* READ IMAGE DATA */
   
    if (fp != NULL) {
		MPI_File_set_view (*fp, file_offset, MPI_UNSIGNED_CHAR, Filetype, "native", MPI_INFO_NULL);
		MPI_File_read(*fp, image, d*Pw*Ph, MPI_UNSIGNED_CHAR, MPI_STATUSES_IGNORE);
		MPI_File_close(fp);
	}else{
		printf("Empty file...The programm will exit");
		exit(1);
	}
		
    input_imag = (float*)malloc(sizeof(float)*d*(Ph+2)*(Pw+2));
    output_imag = (float*)malloc(sizeof(float)*d*(Ph+2)*(Pw+2));

    /* INITIALIZE MATRICES */
    //assign -1 to the outer matrix and the proper value from initial processor
    int i,j,k,m;
    for (i=0; i<Ph+2; i++){
        for (j=0; j<Pw+2; j++){
            if (i==0 || i==Ph+1 || j==0 || j==Pw+1){
                for (m=0; m<d; m++) {
                    input_imag[i*d*(Pw+2)+d*j+m]=-1;
                }
            }
            else{
                for (m=0; m<d; m++) {
                    input_imag[i*d*(Pw+2)+d*j+m] = (float)image[(i-1)*d*Pw+d*(j-1)+m];
                }
            }
            for (m=0; m<d; m++) {
                output_imag[i*d*(Pw+2)+d*j+m]=-1;
            }
        }
    }

    MPI_Barrier(cart_comm);

    time1 = MPI_Wtime();
	
    for (Iter=1; Iter<=TotIter; Iter++){

        if (world_rank == 0)
            printf("iteration:%d\n", Iter);


        /* LOCAL COMMUNICATION */
        MPI_Datatype vector;
        MPI_Type_vector(Ph, d, d*(Pw+2), MPI_FLOAT, &vector);
        MPI_Type_commit(&vector);

        MPI_Request request[16];
        MPI_Status status[16];
        //up-down local communication
        MPI_Irecv(input_imag+d, d*Pw, MPI_FLOAT, up_rank, 0, cart_comm, &request[0]);
        MPI_Isend(input_imag+d*(Pw+2)*Ph+d, d*Pw, MPI_FLOAT, down_rank, 0, cart_comm, &request[1]);
        MPI_Irecv(input_imag+d*(Pw+2)*(Ph+1)+d, d*Pw, MPI_FLOAT, down_rank, 1, cart_comm, &request[2]);
        MPI_Isend(input_imag+d*(Pw+2)+d, d*Pw, MPI_FLOAT, up_rank, 1, cart_comm, &request[3]);
        //left right local commnication
        MPI_Irecv(input_imag+d*(Pw+2), 1, vector, left_rank, 2, cart_comm, &request[4]);
        MPI_Isend(input_imag+2*d*(Pw+2)-2*d, 1, vector, right_rank, 2, cart_comm, &request[5]);
        MPI_Irecv(input_imag+2*d*(Pw+2)-d, 1, vector, right_rank, 3, cart_comm, &request[6]);
        MPI_Isend(input_imag+d*(Pw+2)+d, 1, vector, left_rank, 3, cart_comm, &request[7]);
        //corner communication
        MPI_Irecv(input_imag, d, MPI_FLOAT, up_left_rank, 4, cart_comm, &request[8]);
        MPI_Isend(input_imag+d*(Pw+2)*(Ph+1)-2*d, d, MPI_FLOAT, down_right_rank, 4, cart_comm, &request[9]);
        MPI_Irecv(input_imag+d*(Pw+2)-d, d, MPI_FLOAT, up_right_rank, 5, cart_comm, &request[10]);
        MPI_Isend(input_imag+d*(Pw+2)*Ph+d, d, MPI_FLOAT, down_left_rank, 5, cart_comm, &request[11]);
        MPI_Irecv(input_imag+d*(Pw+2)*(Ph+1), d, MPI_FLOAT, down_left_rank, 6, cart_comm, &request[12]);
        MPI_Isend(input_imag+2*d*(Pw+2)-2*d, d, MPI_FLOAT, up_right_rank, 6, cart_comm, &request[13]);
        MPI_Irecv(input_imag+d*(Pw+2)*(Ph+2)-d, d, MPI_FLOAT, down_right_rank, 7, cart_comm, &request[14]);
        MPI_Isend(input_imag+d*(Pw+2)+d, d, MPI_FLOAT, up_left_rank, 7, cart_comm, &request[15]);


        /* Filter inner pixels before local communication finish */
       // filteringInnerImage(input_imag, output_imag);
		int i,j,n,m,l,k;
      /* COMPUTE INNER DATA */
      /* inner data can be computated before communication finish */
		for (i=2; i<Ph; i++) {
			for (j=2; j<Pw; j++) {
				for (l=0; l<d; l++) {
					 output_imag[ i*d*(Pw+2)+d*j+l ] = input_imag[(i+0-1)*d*(Pw+2)+d*(j+0-1)+l] * h[0][0]
													 + input_imag[(i+0-1)*d*(Pw+2)+d*(j+1-1)+l] * h[0][1]
													 + input_imag[(i+0-1)*d*(Pw+2)+d*(j+2-1)+l] * h[0][2]
													 + input_imag[(i+1-1)*d*(Pw+2)+d*(j+0-1)+l] * h[1][0]
													 + input_imag[(i+1-1)*d*(Pw+2)+d*(j+1-1)+l] * h[1][1]
													 + input_imag[(i+1-1)*d*(Pw+2)+d*(j+2-1)+l] * h[1][2]
													 + input_imag[(i+2-1)*d*(Pw+2)+d*(j+0-1)+l] * h[2][0]
													 + input_imag[(i+2-1)*d*(Pw+2)+d*(j+1-1)+l] * h[2][1]
													 + input_imag[(i+2-1)*d*(Pw+2)+d*(j+2-1)+l] * h[2][2];
				}
			}
		}

        /* wait for local communication to finish */
        MPI_Waitall(16, request, status);

        /* Filter outter pixels after local communication finish */
       // filteringOutterImage(input_imag, output_imag);
                
		float x[3][3], temp;

		for (i=0; i<=1; i++) {
			for (j=1; j<=Pw; j++) {

				for (l=0; l<3; l++) {
					k = d*(Pw+2+j)+i*d*(Pw+2)*(Ph-1)+l;
					temp = 0;
					for (n=0; n<3; n++) {
						for (m=0; m<3; m++) {

							x[n][m] = input_imag[k+(n-1)*d*(Pw+2)+d*(m-1)];
							if (x[n][m]==-1) {
								x[n][m] =  input_imag[k];
							}
							temp += x[n][m] * h[n][m];
						}
					}
					output_imag[k] = temp;
				}
			}
		}


		/* filter columns */
		for (i=2; i<Ph; i++) {
			for (j=0; j<=1; j++) {
				for (l=0; l<3; l++) {
					k = i*d*(Pw+2)+d+j*d*(Pw-1)+l;
					temp = 0;
					for (n=0; n<3; n++) {
						for (m=0; m<3; m++) {
							x[n][m] = input_imag[k+(n-1)*d*(Pw+2)+d*(m-1)];
							if (x[n][m]==-1) {
								x[n][m] =  input_imag[k];
							}
							temp += x[n][m] * h[n][m];
						}
					}
					output_imag[k] = temp;
				}
			}
		}
       

        /* swap pointers of new and old sub-images */
        temp_pointer = input_imag;
        input_imag = output_imag;
        output_imag = temp_pointer;


        /* GLOBAL COMMUNICATION */
        if ((Iter) % 10 == 0){
            int check;
            check = 0;
            i=0; j=0;
            while (check==0 && i<Ph && j<Pw) {
                k = Pw+3+j+i*(Pw+2);
                if ((int)*(input_imag+k)!=(int)*(output_imag+k)){
                    check = 1;
                }
                j++;
                if (j==Pw){
                    i++; j=0;
                }
            }
//	            printf("n:%d rank:%d check:%d\n", n, world_rank, check);
            int tot_check;
            MPI_Reduce(&check, &tot_check, 1, MPI_INT, MPI_MAX, 0, cart_comm);
            if (world_rank == 0) {
                printf("Iter:%d total_check:%d\n", Iter, tot_check);
            }
        }
    }

    time2 = MPI_Wtime();


    for (i=0; i<Ph; i++) {
        for (j=0; j<Pw; j++) {
            for (m=0; m<d; m++) {
                *(image+i*d*Pw+d*j+m) = (unsigned char)*(input_imag+d*(Pw+3)+i*d*(Pw+2)+d*j+m);
            }
        }
    }

    /* PRINT OUTPUT TO FILE */
    if (atoi(argv[1]) == 3)
		MPI_File_open(cart_comm, "images/image_out_RGB.raw", MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, fp);
	else 
		MPI_File_open(cart_comm, "images/image_out_gray.raw", MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, fp);
		
    MPI_File_set_view (*fp, file_offset, MPI_UNSIGNED_CHAR, Filetype, "native", MPI_INFO_NULL);
    MPI_File_write(*fp, image, d*Pw*Ph, MPI_UNSIGNED_CHAR, MPI_STATUSES_IGNORE);
    MPI_File_close(fp);


	//Calculation of  Total and Computation time, 
    time3 = MPI_Wtime();
    comp_time = time2-time1;
    tot_time = time3-time0;
    MPI_Barrier(cart_comm);
    printf("computation time =% lf\n", comp_time);
    MPI_Barrier(cart_comm);
    printf("total time = %lf\n", tot_time);


    MPI_Reduce(&comp_time, &max_comp_time, 1, MPI_DOUBLE, MPI_MAX, 0, cart_comm);
    MPI_Reduce(&tot_time, &max_tot_time, 1, MPI_DOUBLE, MPI_MAX, 0, cart_comm);
    if (world_rank == 0) {
        printf("maximum computation time = %lf\n", max_comp_time);
        printf("maximum total time = %lf\n", max_tot_time);
    }

    MPI_Finalize();

    return 0;
}
//=====================================================================================
// The below functions has been implemented for better readabilty of the code, but since using them we 
//	achive worse computation time than before , we decided to calculate the convolution (3 nested for) inside 
//	of the main programm 
//=====================================================================================

/*
void filteringInnerImage (float *input_imag, float *output_imag) {
    int i,j,l;
      // COMPUTE INNER DATA 
      // inner data can be computated before communication finish 
    for (i=2; i<Ph; i++) {
        for (j=2; j<Pw; j++) {
            for (l=0; l<d; l++) {
                 output_imag[ i*d*(Pw+2)+d*j+l ] = input_imag[(i+0-1)*d*(Pw+2)+d*(j+0-1)+l] * h[0][0]
                                                 + input_imag[(i+0-1)*d*(Pw+2)+d*(j+1-1)+l] * h[0][1]
                                                 + input_imag[(i+0-1)*d*(Pw+2)+d*(j+2-1)+l] * h[0][2]
                                                 + input_imag[(i+1-1)*d*(Pw+2)+d*(j+0-1)+l] * h[1][0]
                                                 + input_imag[(i+1-1)*d*(Pw+2)+d*(j+1-1)+l] * h[1][1]
                                                 + input_imag[(i+1-1)*d*(Pw+2)+d*(j+2-1)+l] * h[1][2]
                                                 + input_imag[(i+2-1)*d*(Pw+2)+d*(j+0-1)+l] * h[2][0]
                                                 + input_imag[(i+2-1)*d*(Pw+2)+d*(j+1-1)+l] * h[2][1]
                                                 + input_imag[(i+2-1)*d*(Pw+2)+d*(j+2-1)+l] * h[2][2];
            }
        }
    }
}

void filteringOutterImage(float *input_imag , float *output_imag) {
    // filter rows and corner pixels
    int i,j,n,m,l,k;
    float x[3][3], temp;

    for (i=0; i<=1; i++) {
        for (j=1; j<=Pw; j++) {

            for (l=0; l<3; l++) {
                k = d*(Pw+2+j)+i*d*(Pw+2)*(Ph-1)+l;
                temp = 0;
                for (n=0; n<3; n++) {
                    for (m=0; m<3; m++) {

                        x[n][m] = input_imag[k+(n-1)*d*(Pw+2)+d*(m-1)];
                        if (x[n][m]==-1) {
                            x[n][m] =  input_imag[k];
                        }

                        temp += x[n][m] * h[n][m];
                    }

                }
                output_imag[k] = temp;
            }

        }
    }


    // filter columns 
    for (i=2; i<Ph; i++) {
        for (j=0; j<=1; j++) {
            for (l=0; l<3; l++) {
                k = i*d*(Pw+2)+d+j*d*(Pw-1)+l;
                temp = 0;
                for (n=0; n<3; n++) {
                    for (m=0; m<3; m++) {

                        x[n][m] = input_imag[k+(n-1)*d*(Pw+2)+d*(m-1)];
                        if (x[n][m]==-1) {
                            x[n][m] =  input_imag[k];
                        }

                        temp += x[n][m] * h[n][m];
                    }

                }
                output_imag[k] = temp;
            }

        }
    }

}
*/
