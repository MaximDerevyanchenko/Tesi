/****************************************************************************
 *
 ****************************************************************************/

#include "hpc.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#define BLKDIM 32

typedef struct{
	int h,w;
	int connectivity;
	int *image;
	int *result;
}bitmap;

/**
 * Reads the input bitmap.
 */
void read_input( FILE *f, bitmap *bm )
{
    int i, con, h, w, pixel;
	int size;
    int *pixels;
	int *positions;
    
    if ( 1 != fscanf(f, "%d", &con) ) {
        fprintf(stderr, "Error: cannot read connectivity\n");
        exit(EXIT_FAILURE);
    }
    if ( con != 4 && con != 8 ) {
        fprintf(stderr, "Error: connectivity cannot be %d. The only acceptable values are 4 and 8\n", con);
        exit(EXIT_FAILURE);
    }
    bm->connectivity = con;
    if ( 2 != fscanf(f, "%d %d", &h, &w) ){
		fprintf(stderr, "Error: cannot to read image sizes\n");
		exit(EXIT_FAILURE);
	}
    assert(h > 0);
    assert(w > 0);
    bm->h = h;
    bm->w = w;

	size = h * w;

    pixels = (int*)malloc( size * sizeof(int) );
	positions = (int*)malloc( size * sizeof(int) );
    assert(pixels);
	assert(positions);

    for (i=0; i<size; i++) {
		if (1 != fscanf(f, "%d ", &pixel)) {
			fprintf(stderr, "Error: cannot read the value of pixel at %d, %d\n", i, size % i);
			exit(EXIT_FAILURE);
		}
		assert(pixel == 0 || pixel == 1);
		pixels[i] = pixel;
		positions[i] = i;
    }
    bm->image = pixels;
    bm->result = positions;
}


/**
 * Reads the bitmap and creates the equivalences found in rows 
 */
__global__ void row_equivalences( int *input, int *res, int w, int h ){

    int i = blockIdx.y * blockDim.y + threadIdx.y;
	int j = blockIdx.x * blockDim.x + threadIdx.x;

	//The position is increased by one, so the position 0 is reserved to the background value and it's easier to access the other values later.
    int pos = (j + i * w) + 1;
	int left = (j - 1) + i * w;

	if ( i == 0 && j == 0 ){
		res[0] = 0;
	}

	if ( i < h && j < w ){
    	if ( j == 0 ){
        	res[pos] = pos * input[pos - 1];
		} else {
			if ( input[pos - 1] ){
				if ( input[left] ){
					res[pos] = pos-1;
				} else {
					res[pos] = pos;
				}
			} else {
				res[pos] = 0;
			}
		}
	}
}

/**
 * Reads the bitmap and creates the equivalences found in coulmns
 */
__global__ void col_equivalences( int *input, int *res, int w, int h ){

    int i = blockIdx.y * blockDim.y + threadIdx.y;
	int j = blockIdx.x * blockDim.x + threadIdx.x;

	//The position is increased by one, so the position 0 is reserved to the background value and it's easier to access the other values later.
    int pos = (j + i * w) + 1;
	int upper = j + (i - 1) * w;

	if ( i == 0 && j == 0 ){
		res[0] = 0;
	}

	if ( i < h && j < w ){
		if ( i == 0 ){
			res[pos] = pos * input[pos - 1];
		} else {
			if ( input[pos - 1] ){
				if ( input[upper] ){
					res[pos] = pos-(w);
				} else {
					res[pos] = pos;
				}
			} else {
				res[pos] = 0;
			}
		}
	}
}

/**
 * Sets the labels of the result from the vector labels
 */
__global__ void set_labels( int *input, int *result, int *labels, int w, int h ){
	int i = blockIdx.y * blockDim.y + threadIdx.y;
	int j = blockIdx.x * blockDim.x + threadIdx.x;
	//The position is increased by one, so the position 0 is reserved to the background value and it's easier to access the other values later.
    int pos = (j + i * w) + 1;
	int root = labels[pos];

		while ( root != labels[root] ){
			root = labels[root];
		}
		labels[pos] = root;
		result[pos - 1] = root;
		
}

/**
 * Merging column values and row values, getting the provvisional label values
 */
__global__ void merge( int *input, int *result, int *row, int *col, int *label, int w, int h){

	int i = blockIdx.y * blockDim.y + threadIdx.y;
	int j = blockIdx.x * blockDim.x + threadIdx.x;
	//The position is increased by one, so the position 0 is reserved to the background value and it's easier to access the other values later.
    int pos = (j + i * w) + 1;
	if ( i < h && j < w ){
		label[pos] = 0;
		if ( input[pos - 1] ){
			if ( col[pos] == pos && row[pos] == pos ){
				label[pos] = pos;
			} else if ( col[pos] && col[pos] != pos){
				label[pos] = col[pos];
			} else if ( row[pos] && row[pos] != pos ){
				label[pos] = row[pos];
			}
		} else {
			result[pos - 1] = 0;
		}
	}
}

/**
 * Checks if root leads to test_root
 */
int check_cycles( int *labels, int test_root, int root ){

	while ( labels[root] != root && labels[root] != test_root ){
		root = labels[root];
	}

	return labels[root] != test_root;
}

/**
 * Checks the roots and updates if no cycles are found
 */
void update_label( int root_to_update, int root_to_set, int *labels ){
	int is_acyclic;

	//Searching the root of the label to update
	while ( root_to_update != labels[root_to_update] && root_to_set != root_to_update ){
		root_to_update = labels[root_to_update];
	}

	is_acyclic = check_cycles( labels, root_to_update, root_to_set );
	if ( labels[root_to_set] != root_to_update && root_to_update != root_to_set && is_acyclic){
		labels[root_to_update] = root_to_set;
	}
}

/**
 * Creates the result matrix with labels from an initial binary matrix in connectivity 4
 */
void label( bitmap *bm ){
	int i, j;
	int pos;
	int upper;
	int left;
	int upper_left;
	int root, my_root;
	int size;

	int *label;

	int *d_input, *d_result;
	int *d_row, *d_col, *d_label;

	dim3 block( BLKDIM, BLKDIM );
	dim3 grid( (bm->w + BLKDIM - 1)/BLKDIM, (bm->h + BLKDIM - 1)/BLKDIM );

	size = bm->h * bm->w;

	size_t byte_size = size * sizeof(int);
	size_t label_size = (size + 1) * sizeof(int);

	label = (int*)malloc( (size + 1) * sizeof(int) );
	
	//Input and result contain the matrix as it is while row, col and label have an offset to contain the background value (0) at position 0
	cudaSafeCall( cudaMalloc((int**) &d_input, byte_size) );
	cudaSafeCall( cudaMalloc((int**) &d_result, byte_size) );
	//The size is increased by 1, so position 0 is reserved to the 0 value (background) and all the values are easily reached later
	cudaSafeCall( cudaMalloc((int**) &d_row, label_size) );
	cudaSafeCall( cudaMalloc((int**) &d_col, label_size) );
	cudaSafeCall( cudaMalloc((int**) &d_label, label_size) );

	cudaSafeCall( cudaMemcpy(d_input, bm->image, byte_size, cudaMemcpyHostToDevice) );

	row_equivalences<<<grid, block>>>(d_input, d_row, bm->w, bm->h);
	cudaCheckError();
	col_equivalences<<<grid, block>>>(d_input, d_col, bm->w, bm->h);
	cudaCheckError();

	//Merging column values and row values, getting the provvisional label values
	merge<<<grid, block>>>(d_input, d_result, d_row, d_col, d_label, bm->w, bm->h);
	cudaCheckError();

	cudaSafeCall( cudaFree( d_row ) );
	cudaSafeCall( cudaFree( d_col ) );

	//Setting provisional labels
	set_labels<<<grid, block>>>(d_input, d_result, d_label, bm->w, bm->h);
	cudaCheckError();

	cudaSafeCall( cudaMemcpy( label, d_label, label_size, cudaMemcpyDeviceToHost ) );

	/**
	 * Refining of the local equivalences
	 * In particular removing the equivalences of the following type
	 * 
	 * 			* * * * *
	 *			* * a * *
	 *			* b a * *
	 *			* * * * *
	 * 
	 * This is the only type of the equivalences that could be found after the previous passes of the algorithm.
	 */
	for (i = 1; i < bm->h; i++){
		for (j = 1; j < bm->w; j++){
			pos = (j + i * bm->w) + 1;
			upper = j + (i-1) * bm->w;
			left = (j - 1) + i * bm->w;
			upper_left = (j-1) + (i-1) * bm->w;
			if ( bm->image[pos - 1] ){
				my_root = label[pos];
				if ( bm->image[upper] && bm->image[left] && !bm->image[upper_left] ){
					root = label[pos-1];
					update_label(my_root, root, label);
				}
			}
		}
	}

	cudaSafeCall( cudaMemcpy( d_label, label, label_size, cudaMemcpyHostToDevice ) );
	
	//Setting final labels
	set_labels<<<grid, block>>>(d_input, d_result, d_label, bm->w, bm->h);
	cudaCheckError();

	cudaSafeCall( cudaMemcpy(bm->result, d_result, byte_size, cudaMemcpyDeviceToHost) );

	free(label);
	cudaSafeCall( cudaFree( d_label ) );
}

/**
 * Free the bitmap data structure
 */
void free_bitmap( bitmap *bm )
{
	free( bm->image );
	free( bm->result );
}

int main( void )
{
    bitmap bm;
    double tstart, elapsed;
	int x = 0, y = 0;
	int pos = 0, left = 0, upper = 0;
	int val = 0;
	int isCorrect = 1;
    
    read_input(stdin, &bm);

	tstart = hpc_gettime();
	label(&bm);
	elapsed = hpc_gettime() - tstart;

    fprintf(stderr, "Elapsed time = %f sec\n", elapsed);
	
	//Printing the result
	printf("Result:\n");
    for ( int i = 0; i < bm.h; i++ ){
		for ( int j = 0; j < bm.w; j++ ){
			pos = j + i * bm.w;
			left = (j - 1) + i * bm.w;
			upper = j + (i - 1) * bm.w;
			printf("%9d ", bm.result[pos]);
			if ( i && j ){
				//Checking if the result provided by the algorithm is correct
				if ( bm.image[pos] ){
					if ( (bm.image[pos] == 1 && bm.result[pos] == 0) || (bm.image[pos] == 0 && bm.result[pos] != 0) ){
						x = j;
						y = i;
						val = bm.result[pos];
						isCorrect = 0;
					}
					if ( bm.result[left] && bm.result[left] != bm.result[pos]){
						x = j;
						y = i;
						val = bm.result[pos];
						isCorrect = 0;
					}
					if ( bm.result[upper] && bm.result[upper] != bm.result[pos]){
						x = j;
						y = i;
						val = bm.result[pos];
						isCorrect = 0;
					}
				}
			}
		}
		printf("\n");
	}

	if ( isCorrect ){
		fprintf(stderr, "Correct\n");
	} else {
		fprintf(stderr, "Result Wrong. Last wrong value found: X = %d, Y = %d, val = %d\n", x, y, val);
	}

    free_bitmap(&bm);
    return EXIT_SUCCESS;    
}
