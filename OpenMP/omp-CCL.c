/****************************************************************************
 *
 ****************************************************************************/

#include "hpc.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

typedef struct{
	int h,w;
	int connectivity;
	int **image;
	int **result;
}bitmap;

/**
 * Reads the input bitmap.
 */
void read_input( FILE *f, bitmap *bm )
{
    int i, j, con, h, w, pixel;
    int **pixels;
	int **positions;
    
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
    pixels = (int**)malloc( h * sizeof(int*) );
	positions = (int**)malloc( h * sizeof(int*) );
    assert(pixels);
	assert(positions);
    for (i=0; i<h; i++) {
		pixels[i] = (int*)malloc( w * sizeof(int) );
		positions[i] = (int*)malloc( w * sizeof(int) );
		for (j=0; j<w; j++){
			if (1 != fscanf(f, "%d ", &pixel)) {
				fprintf(stderr, "Error: cannot read the value of pixel at %d, %d\n", i, j);
				exit(EXIT_FAILURE);
			}
			assert(pixel == 0 || pixel == 1);
			pixels[i][j] = pixel;
			positions[i][j] = j + i * w;
		}
    }
    bm->image = pixels;
    bm->result = positions;
}


/**
 * Reads the bitmap and creates the equivalences found in rows 
 */
void row_equivalences( bitmap *bm, int *eq, int from, int to ){
    int i, j;
    int pos;

    for (i = from; i < to; i++){
		//The position is increased by one, so the position 0 is reserved to the background value and it's easier to access the other values later.
        eq[i * bm->w + 1] = (i * bm->w + 1) * bm->image[i][0];
		for (j = 1; j < bm->w; j++){
			pos = (j + i * bm->w) + 1;
			if ( bm->image[i][j] ){
				if ( bm->image[i][j-1] ){
					eq[pos] = pos-1;
				} else {
					eq[pos] = pos;
				}
			} else {
				bm->result[i][j] = 0;
				eq[pos] = 0;
			}
		}
	}
}

/**
 * Reads the bitmap and creates the equivalences found in coulmns
 */
void col_equivalences( bitmap *bm, int *eq, int from, int to, int *starts, int thread_num ){
    int i, j;
    int pos;
	int flag = 0;
	int special_row = 0;
	int count = 0;

	for( count = 0; count < thread_num; count++ ){
		if (starts[count] >= from){
			special_row = starts[count];
			count++;
			break;
		}
	}

    for (i = from; i < to; i++){
		for (j = 0; j < bm->w; j++){
			
			if ( i == special_row ){
				eq[i * bm->w + j + 1] = (i * bm->w + j + 1) * bm->image[i][j];
				flag = 1;
			} else {
				//The position is increased by one, so the position 0 is reserved to the background value and it's easier to access the other values later.
				pos = (j + i * bm->w) + 1;
				if ( bm->image[i][j] ){
					if ( bm->image[i-1][j] ){
						eq[pos] = pos-(bm->w);
					} else {
						eq[pos] = pos;
					}
				} else {
					eq[pos] = 0;
				}
			}
		}
		if (flag && count < thread_num){
			special_row = starts[count];
			count++;
			flag = 0;
		}
	}
}

/**
 * Sets the labels of the result from the vector labels
 */
void set_labels( bitmap *bm, int *labels ){
	int i, j;
	int pos, root;

#pragma omp for collapse(2)
	for (i = 0; i < bm->h; i++){
		for (j = 0; j < bm->w; j++){
            pos = (j + i * bm->w) + 1;
			root = labels[pos];
			while ( root != labels[root] ){
				root = labels[root];
			}
			labels[pos] = root;
			bm->result[i][j] = root;
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
void update_label( bitmap *bm, int root_to_update, int root_to_set, int *labels ){
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
	int *row, *col, *label;
	int size;
	int pos;
	int root, my_root;

	size = bm->h * bm->w;
	
	//the size is increased by 1, so position 0 is reserved to the 0 value (background) and all the values are easily reached later
	row = (int*)malloc( (size + 1) * sizeof(int) );
    col = (int*)malloc( (size + 1) * sizeof(int) );
    label = (int*)malloc( (size + 1) * sizeof(int) );
	row[0] = col[0] = label[0] = 0;

	int max_thread_num = omp_get_max_threads();
	int starts[max_thread_num];
	max_thread_num = bm->h > max_thread_num ? max_thread_num : bm->h;

#pragma omp parallel default(none) shared(bm, row, col, label, starts, stderr) private(my_root, root, pos, i, j) num_threads(max_thread_num)
{
	int my_id = omp_get_thread_num();
	int my_group = my_id % 2;
	int my_group_id = my_id / 2;
    int thread_num = omp_get_num_threads();
	int my_group_size = thread_num / 2;
	int my_group_start;
	int my_group_end;
	int my_start;
	int my_end;

	if (thread_num % 2 == 1 && my_group == 1){
		my_group_size++;
	}

	/**
	 * Dividing the matrix based on the columns and every thread executes the CCL algorithm on its section.
	 * The main problem will be with wide pictures with very low amount of rows (less than thread amount).
	 * This case is hardly met, because the CCL algorithm is used on images not lines, so it makes no problem for our situation.
	 * 
	 */
	my_start = bm->h * my_id / thread_num;
	my_end = bm->h * (my_id + 1) / thread_num;
	
	starts[my_id] = my_start;
	
#pragma omp barrier

	if ( thread_num > 1 ){
		my_group_start = bm->h * my_group_id / my_group_size;
		my_group_end = bm->h * (my_group_id + 1) / my_group_size;
		if ( my_group == 0 ){
			row_equivalences(bm, row, my_group_start, my_group_end);
		} else {
			col_equivalences(bm, col, my_group_start, my_group_end, starts, thread_num);
		}
	} else {
		row_equivalences(bm, row, my_start, my_end);
		col_equivalences(bm, col, my_start, my_end, starts, thread_num);
	}

#pragma omp barrier
	//Merging column values and row values, getting the provvisional label values
#pragma omp for collapse(2)
	for (i = 0; i < bm->h; i++){
		for (j = 0; j < bm->w; j++){
			pos = (j + i * bm->w) + 1;
            label[pos] = 0;
            if ( bm->image[i][j] ){
                if ( col[pos] == pos && row[pos] == pos ){
                    label[pos] = pos;
                } else if ( col[pos] && col[pos] != pos){
                    label[pos] = col[pos];
                } else if ( row[pos] && row[pos] != pos ){
                    label[pos] = row[pos];
                }
            }
		}
	}

#pragma omp barrier
#pragma omp master
{
	free(row);
	free(col);
}
	//Setting provisional labels
	set_labels(bm, label);

#pragma omp barrier

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
	for (i = my_start + 1; i < my_end; i++){
		for (j = 1; j < bm->w; j++){
			if (my_id == 11){
			}
			if ( bm->image[i][j] ){
				pos = (j + i * bm->w) + 1;
				my_root = label[pos];
				if ( bm->image[i-1][j] && bm->image[i][j-1] && !bm->image[i-1][j-1] ){
					root = label[pos-1];
					update_label(bm, root, my_root, label);
				}
			}
		}
	}

#pragma omp barrier

	/**
	 * Resolving the equivalences of the borders of the sections of the image and calcolating the final label values
	 * 
	 */
	for ( i = 1; i < thread_num; i++ ){
#pragma omp for
		for ( j = 0; j < bm->w; j++ ){
			if ( bm->image[starts[i]][j] && starts[i] > 0 ){
				pos = (j + starts[i] * bm->w) + 1;
				my_root = label[pos];
				if ( bm->image[starts[i]][j] && bm->image[starts[i] - 1][j] ){
					root = label[pos - bm->w];
					#pragma omp critical
					{
					update_label(bm, my_root, root, label);
					}
				}
			}
		}
	}

	//Setting final labels
	set_labels(bm, label);
}
	free(label);
}

/**
 * Free the bitmap data structure
 */
void free_bitmap( bitmap *bm )
{
	int i;
	for ( i = 0; i < bm->h; i++){
		free( bm->image[i] );
		free( bm->result[i] );
	}
	free( bm->image );
	free( bm->result );
}

int main( void )
{
    bitmap bm;
    double tstart, elapsed;
	int x = 0, y = 0;
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
			printf("%9d ", bm.result[i][j]);
			if ( i && j ){
				//Checking if the result provided by the algorithm is correct
				if ( bm.result[i][j] ){
					if ( (bm.image[i][j] == 1 && bm.result[i][j] == 0) || (bm.image[i][j] == 0 && bm.result[i][j] != 0) ){
						x = j;
						y = i;
						val = bm.result[i][j];
						isCorrect = 0;
					}
					if ( bm.result[i-1][j] && bm.result[i-1][j] != bm.result[i][j]){
						x = j;
						y = i;
						val = bm.result[i][j];
						isCorrect = 0;
					}
					if ( bm.result[i][j-1] && bm.result[i][j-1] != bm.result[i][j]){
						x = j;
						y = i;
						val = bm.result[i][j];
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
		fprintf(stderr, "Result Wrong. Last wrong value found: X = %d, Y = %d, val = %d", x, y, val);
	}

    free_bitmap(&bm);
    return EXIT_SUCCESS;    
}
