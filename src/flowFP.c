/*
 * FILE: flowFP.c
 *
 * This module contains the functions called by the flowFP
 * R package to preform rectangular probability binning on
 * multi-dimensional FCS list mode data.
 *
 * 
 * Terms:
 *   parameter – Column in the data-set matrix. corresponds to 
 *               a detector in the flow cytometer.
 *   bin       – a subset of the events.
 *   bin model – is a set of parameters, and median values used
 *               to divide up a data-set. The model is organized
 *               as an array (num_levels long) of vectors. Each
 *               vector contains all of the split axis and values.
 *               The len of each vector is 2^level. (2,4,8,16...)
 *
 * Create a bin model - 
 *  Goal is to create a model for a data-set where we divide up
 *  the multi-dimensional space such that each region (bin) has
 *  an equal probability of containing events. In the simplest
 *  case if we had one dimensional data, we could calculate the
 *  median, and divide the set into roughly two even sections.
 *  With multi-dimensional data we first need to decide which 
 *  parameter (axis) to split on. We split on the parameter with
 *  the largest variance. Each time we split a parent bin into
 *  it's two child bins, we record the split value and parameter
 *  (axis), this gets stored in the model. We keep track of which
 *  events belong to which bin, by using the tags vector. The tags
 *  vector is as long as the number of events in the data-set, and
 *  contains the bin number for each event as we go through the
 *  binning process. 'tags' is initialized to all 1's the first
 *  time bin_level is called, and it gets changed to having 1's
 *  and 2's, etc.

 * Applying a model – Goal is to project a data-set onto an existing
 *  model. We use the model described above, to bin the data. We have
 *  to go level by level, first dividing based on the level one
 *  split_axis, and the split_val. Then we divide those two set into
 *  four, etc. We keep track of which events belong to each bin in the
 *  tags vector as described above.
 */
 
#include "flowFP.h"

/*-----------------------------------------------------------------
   interface to R with arguments:
     data :      matrix of fcs events 'provided by a call to exprs()'
                 read only.
     tags:       Results are returned in this vector of integers. The tags
                 vector keeps track of the assignment of an event to a bin
                 between calls to bin level.
     split_axis: Results are returned in this vector of integers, and records
                 which axis was split on.
     split_val:  Results are returned in this vector of numerics, and records
                 the split threshold values.
     level:      Current Level, counter from the caller's loop, (i.e. i)
     param_idx:  Read only vector of indexes, 1 based, parameters to consider
                 during the binning process.
------------------------------------------------------------------*/

SEXP bin_level(SEXP _data, SEXP _tags, SEXP _split_axis, SEXP _split_val, 
               SEXP _level, SEXP _param_idx)
{ 

	CMATRIX fcs;
	CIVECTOR tags, split_axis, param_idx;
	CVECTOR split_val;

	Rmatrix2C(_data, &fcs);
	Rvector2C(_split_val, &split_val);
	
	Rivector2C(_tags, &tags);
	Rivector2C(_split_axis, &split_axis);
	Rivector2C(_param_idx, &param_idx);

	int level = INTEGER_VALUE(_level);
	int num_splits = (int)(0x1) << (level - 1);
	
	CIVECTOR **event_idx_lut = create_idx_lut(&tags, num_splits);
	
	int i;
	int bin_number = 1;
	for(i=0; i < num_splits; i++) {
		if(split_axis.data[i] == 0) {
			split_axis.data[i] = find_split_axis(&fcs, &param_idx, event_idx_lut[i]);
			split_val.data[i] = find_split_val(&fcs, split_axis.data[i], event_idx_lut[i]);
		}
		split_bin(&fcs, &tags, event_idx_lut[i], split_axis.data[i], split_val.data[i], bin_number);
		bin_number += 2;
	}
		
	destroy_idx_lut(event_idx_lut, num_splits);
	return(NULL);
}

/*-----------------------------------------------------------------
   interface to R with arguments:
     fcs :       matrix of fcs events 'provided by a call to exprs()'
                 read only.
     level:      Current Level, counter from the caller's loop, (i.e. i)
     split_axis: Contains the 1 based axis numbers used to split this level.
     split_val:  Contains the threshold value as a double, used to divide
                 the event according to the model.
     tags:       Results are returned in this vector of integers. The tags
                 vector keeps track of the assignment of an event to a bin
                 between calls to bin level.
------------------------------------------------------------------*/

SEXP tag_events(SEXP _fcs, SEXP _level, SEXP _split_axis, SEXP _split_val, SEXP _tags) {

	CMATRIX fcs;
	int level = INTEGER_VALUE(_level);
	CIVECTOR split_axis, tags;
	CVECTOR split_val;
	
	Rmatrix2C(_fcs, &fcs);
	Rivector2C(_split_axis, &split_axis);
	Rvector2C(_split_val, &split_val);
	Rivector2C(_tags, &tags);
	int num_splits = (int)(0x1) << (level - 1);
	
	CIVECTOR **event_idx_lut = create_idx_lut(&tags, num_splits);
		
	int i;
	int bin_number = 1;
	for(i=0; i < num_splits; i++) {
		if(split_axis.data[i] == 0) {
			error("This model is broken\n");
		}
		split_bin(&fcs, &tags, event_idx_lut[i], split_axis.data[i], split_val.data[i], bin_number);
		bin_number += 2;
	}
		
	destroy_idx_lut(event_idx_lut, num_splits);
	return (NULL);
}

/*-----------------------------------------------------------------
   interface to R with arguments:
      cnts:      Results are returned in this vector of integers. The cnts
                 vector keeps track of number of an events in each bin. It
                 is a integer histogram.
                 between calls to bin level.
     tags:       This vector of integers, contains the assignment of each event
                 to a bin, 1 based.
------------------------------------------------------------------*/

SEXP count_events(SEXP _cnts, SEXP _tags) {

	CIVECTOR cnts, tags;
	Rivector2C(_cnts, &cnts);
	Rivector2C(_tags, &tags);
	int i;
	for(i=0; i < cnts.len; i++)
		cnts.data[i] = 0;
		
	for(i=0; i < tags.len; i++) {
		if (tags.data[i] < 1 || tags.data[i] > cnts.len)
			error("The tags data is out of range for this model, stopped at tag[%d] = %d (range 1 to %d)\n", i, tags.data[i], cnts.len);
			
		cnts.data[tags.data[i] - 1]++;
	}
	
	return(NULL);
}

/*
 * Function: Rmatrix2C
 *
 *  This function converts an rmatrix to a
 *  CMATRIX. CMATRIX is defined by this library
 *  as a pointer to a double and a dimensions.
 */
void Rmatrix2C(SEXP rmatrix, CMATRIX *cmatrix) {
	SEXP dim;
	/* check input argument _data */
	PROTECT(dim = getAttrib(rmatrix, R_DimSymbol));
	if(((!isReal(rmatrix)) & !isInteger(rmatrix)) | isNull(dim) | (LENGTH(dim)!=2))
		error("Invalid argument 'rmatrix': must be a matrix of reals");
	cmatrix->data = REAL(AS_NUMERIC(rmatrix));
	cmatrix->nrows = INTEGER(dim)[0];
	cmatrix->ncols = INTEGER(dim)[1];
	UNPROTECT(1); /* done with dimData */
	return;
}

/*
 * Function: Rivector2C
 *
 *  This function converts an r int vector to a
 *  CIVECTOR. CIVECTOR is defined by this library
 *  as a pointer to a int and a length.
 */
void Rivector2C(SEXP rivect, CIVECTOR *civect) {
		
	if (!isInteger(rivect))
		error("Invalid argument 'rivect': must be a vector of ints");

	civect->data = INTEGER_POINTER(rivect);
	civect->len = LENGTH(rivect);
	return;
}

/*
 * Function: Rvector2C
 *
 *  This function converts an r numeric vector to a
 *  CVECTOR. CVECTOR is defined by this library
 *  as a pointer to a double and a length.
 */
void Rvector2C(SEXP rvect, CVECTOR *cvect) {
		
	if (!isReal(rvect))
		error("Invalid argument 'rvect': must be a vector of Reals");
	if(LENGTH(rvect) <= 0)
		error("Invalid argument 'rvect': has zero length");
	cvect->data = NUMERIC_POINTER(rvect);
	cvect->len = LENGTH(rvect);
	return;
}

