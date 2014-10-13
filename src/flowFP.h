#include <R.h>
#include <Rmath.h>
#include <float.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Print.h>
#include <Rdefines.h>
#include <R_ext/Utils.h> 

#include <stdlib.h>
#include <stdarg.h>


typedef struct _CMATRIX {
	double *data;
	int nrows;
	int ncols;
} CMATRIX;



typedef struct _CVECTOR {
	double *data;
	int len;
} CVECTOR;

typedef struct _CIVECTOR {
	int *data;
	int len;
} CIVECTOR;


#define OFS(mat, row, col)	((mat)->nrows * col + row)
#define INRANGE(v, min, max)	((v) >= (min) && (v) <= (max))

/* proto types
*/
SEXP bin_level(SEXP _data, SEXP _tags, SEXP _split_axis, SEXP _split_val, SEXP _level, SEXP _param_idx);
SEXP tag_events(SEXP _fcs, SEXP _level, SEXP _split_axis, SEXP _split_val, SEXP _tags);
SEXP count_events(SEXP _cnts, SEXP _tags);

void Rmatrix2C(SEXP rmatrix, CMATRIX *cmatrix);
void Rivector2C(SEXP rivect, CIVECTOR *civect);
void Rvector2C(SEXP rvect, CVECTOR *cvect);
CIVECTOR **create_idx_lut(CIVECTOR *tags, int num_splits);
void destroy_idx_lut(CIVECTOR **lut, int num_splits);
double * extract_col(double *buf, CMATRIX *fcs, int col_idx, CIVECTOR *row_idx);
void split_bin(CMATRIX *fcs, CIVECTOR *tags, CIVECTOR *lut, int split_axis, double split_val, int low_bin_num);
double variance(double *data, int num_pts);
double median(double *data, int num_pts);
double find_split_val(CMATRIX *fcs, int split_axis, CIVECTOR *lut);
int find_split_axis(CMATRIX *fcs, CIVECTOR *param_idx, CIVECTOR *lut);
