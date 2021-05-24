#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include <unistd.h>

int readin_grid(char* filename, double **data, int *ra_index, int *r0_index, int *rb_index){
	
    hid_t       file_id;
    hsize_t     dims[1];
    herr_t      status;

	if( access( filename, F_OK ) != -1 ) {
    	/* open file from ex_lite1.c */
		printf("Loading data from %s\n", filename);
    	file_id = H5Fopen (filename, H5F_ACC_RDONLY, H5P_DEFAULT);
	}else{
		printf("Required file %s not found\n", filename);
		exit(0);
	}

    /* get the dimensions of the dataset */
    status = H5LTget_dataset_info(file_id, "/r", dims, NULL, NULL);

    // Allocate  memory for the data
    *data = malloc(dims[0]*sizeof(double));
 
    /* read in grid */
    status = H5LTread_dataset_double(file_id, "/r", *data);
	
	// Read in important indexes
	int ii[3];
	status = H5LTread_dataset_int(file_id, "ImportantIndexes", ii);
	*ra_index = ii[0] - 1;
	*r0_index = ii[1] - 1;
	*rb_index = ii[2] - 1;

    /* close file */
    status = H5Fclose (file_id);
	
	return dims[0];
}
