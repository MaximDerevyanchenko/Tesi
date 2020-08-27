
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <math.h>

typedef struct{
	int h,w;
	int **image;
}bitmap;

void create( bitmap *bm )
{
    int i, j, h, w, pixel;
    int **pixels;
	int val, count;
    
	h = 4320;
	w = 7680;

    bm->h = h;
    bm->w = w;
    pixels = (int**)malloc( h * sizeof(int*) );
	count = 0;
    assert(pixels);
	srand(time(NULL));
    for (i=0; i<h; i++) {
		pixels[i] = (int*)malloc( w * sizeof(int) );
		for (j=0; j<w; j++){
			val = rand() % 10;
			if ( val < 8 ){
				count++;
				pixel = 1;
			} else {
				pixel = 0;
			}
			
			pixels[i][j] = pixel;
		}
    }
	fprintf(stderr, "%d\n", count);
    bm->image = pixels;
}

void free_bitmap( bitmap *bm )
{
	int i;
	for ( i = 0; i < bm->h; i++){
		free( bm->image[i] );
	}
	free( bm->image );
}

int main( void )
{
    bitmap bm;
    
    create(&bm);
	printf("4\n%d %d\n", bm.h, bm.w);
	
	for ( int i = 0; i < bm.h; i++){
		for ( int j = 0; j < bm.w; j++){
			printf("%d ", bm.image[i][j]);
		}
		printf("\n");
	}

    free_bitmap(&bm);
    return EXIT_SUCCESS;    
}
