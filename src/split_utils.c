/*
 * FILE: split_utils.c
 */

#include "flowFP.h"

/*
 * function: split_bin
 *
 *   This function divides a bin into two parts.
 *   It records which event belongs to which bin
 *   by setting values in the tags array.
 *
 */

void split_bin(CMATRIX *fcs, CIVECTOR *tags, CIVECTOR *lut, int split_axis, 
               double split_val, int low_bin_num) {
	int i;
	int col_idx = split_axis - 1; // Check this!!
	int ofs;
	int row_idx;
	int num_low = 0;
	int num_hi = 0;

	for(i=0; i < lut->len; i++) { 
		row_idx = lut->data[i];
		ofs = OFS(fcs, row_idx, col_idx);
		
		if (fcs->data[ofs] > split_val) {
			tags->data[row_idx] = low_bin_num + 1;
			num_hi++;
		} else {
			tags->data[row_idx] = low_bin_num;
			num_low++;
		}
	}
	
}

/*
 * function: find_split_val
 *
 *  Given a split_axis (a.k.a a column in the fcs data matrix)
 *  this routine extracts all of the points in the index look up table
 *  and returns the median of that column, this is the split value for
 *  this bin.
 */

double find_split_val(CMATRIX *fcs, int split_axis, CIVECTOR *lut) {

	int col_idx = split_axis - 1;
	if (col_idx < 0 || split_axis > fcs->ncols)
		printf("find_split_val: split_axis = %d\n", split_axis);
	double *vals = extract_col(NULL, fcs, col_idx, lut);
	double med = median(vals, lut->len);
	free(vals);
	
	return (med);

}

/*
 * Function: find_split_axis
 *
 *  This routine finds the axis (a.k.a column of the fcs matrix)
 *  with the largest variance. It only considers the points indexed
 *  by the index look-up vector.
 */
int find_split_axis(CMATRIX *fcs, CIVECTOR *param_idx, CIVECTOR *lut) {

	if (lut->len == 0)
		error("No events in this parent bin, you are probably trying to sub-divide the data into too many levels\n");
	int i;
	int split_axis = -1;
	double max_var = DBL_MIN;
	double *vals = malloc(sizeof(double) * lut->len);
	int col_idx;
	double tmp;
	
	for(i=0; i < param_idx->len; i++) { // loop over each parameter
		col_idx = param_idx->data[i] - 1;
		extract_col(vals, fcs, col_idx, lut);
		tmp = variance(vals, lut->len);
		if (max_var < tmp) {
			split_axis = col_idx + 1;
			max_var = tmp;
		}
	}
	free(vals);
	if (split_axis == -1)
		error("No events in this bin, you are probably trying to sub-divide the data into too many levels\n");
	return (split_axis);
}

/*
 * Function: extract_col
 *
 *  This routine, extracts a subset of a column out of the data matrix.
 *  and returns it as a array of doubles.So the fcs matrix contains
 *  all of the events in a matrix, and col_idx secifies which column,
 *  and row_idx is a vect that specifies which rows.
 */
double * extract_col(double *buf, CMATRIX *fcs, int col_idx, CIVECTOR *row_idx) {

	double *data;
	if (buf == NULL) {
		data = malloc(sizeof(double) * row_idx->len);
	} else {
		data = buf;
	}
	int i;
	int ofs;
	for(i=0; i < row_idx->len; i++) {
		ofs = OFS(fcs, row_idx->data[i], col_idx);
		data[i] = fcs->data[ofs];
	}
	return data;
}


/*
 * Function: variance
 *
 * This function calculates the variance of a set of data.
 * num_pts is the length of the data pointer.
 */
double variance(double *data, int num_pts) {
	if (data == NULL || num_pts == 0)
		return NAN;

	double val;
	double sum = 0.0, sum_sq = 0.0;
	int i;
	
	for(i=0; i < num_pts; i++) {
		val = (double)*data++;
		sum += val;
		sum_sq += val * val;
	}
		
	return ( (sum_sq - (sum * sum) / num_pts) / (num_pts-1) );
	
}

/*
 * Function: compare
 *
 * This is the compare function only used by qsort
 * called by median.
 */
int compare(const void *i, const void *j)
{ 
  double *a = (double *)i;
  double *b = (double *)j;

  if (*a > *b) return( 1);
  if (*a < *b) return(-1);
  return(0);
}

/*
 * Function: median
 *
 *  This routine calculates the meadian of a set of data.
 *  num_pts it the length of the data pointer.
 */
double median(double *data, int num_pts) {
	
	if (data == NULL || num_pts == 0)
		return NAN;
		
	qsort(data, num_pts, sizeof(double), compare);
	
	double  median;
	if (num_pts % 2 != 0) 
		median = data[num_pts / 2];
    else
		median = (data[num_pts / 2 - 1] + data[num_pts / 2]) / 2.0;
	
	return median;
}

/*
 * Function: destroy_idx_lut
 *
 * This function free's all of the memory allocated by
 * create_idx_lut
 */
void destroy_idx_lut(CIVECTOR **lut, int num_splits) {
	int i;
	for(i=0; i < num_splits; i++) {
		free(lut[i]->data);
		free(lut[i]);
	}
	free(lut);
}

/*
 * Function create_idx_lut
 *
 * This routine creates an array of CVECTORS (num_splits of them)
 * each vectors length corresponds to the number of events of that
 * number in the tags array. The tags vector has numbers in ranging
 * from 1 to num_splits. These values are arranged in event order
 * and this function simply re-organizes them putting all of the
 * 1's together and 2's etc. This enables us to race down a vector
 * of indices and extract the values for a particular bin.
 */

CIVECTOR **create_idx_lut(CIVECTOR *tags, int num_splits) {
	int i;
	int *lengths;
	int tmp;
	

	lengths = (int *)calloc(num_splits, sizeof(int));
	for(i=0; i < tags->len; i++) {
		tmp = tags->data[i] - 1;
		if (!INRANGE(tmp, 0, num_splits - 1))
			error("This tag (%d)[%d] has values either < 0 or >= to %d", tmp, i, num_splits);
		lengths[tmp]++;
	}
	
	CIVECTOR **lut;
	lut = malloc(sizeof(CIVECTOR *) * num_splits);
	for(i=0; i < num_splits; i++) {
		lut[i] = malloc(sizeof(CIVECTOR));
		lut[i]->len = lengths[i];
		lut[i]->data = malloc(sizeof(int) * lengths[i]);
		lengths[i] = 0;
	}
	for(i=0; i < tags->len; i++) {
		tmp = tags->data[i] - 1;
		lut[tmp]->data[lengths[tmp]] = i;
		lengths[tmp]++;
	}
	
	return (lut);
}
