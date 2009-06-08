#include <assert.h>
#include <setjmp.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "CuTest.h"
#include "../flowFP.h"

/*-------------------------------------------------------------------------*
 * Test
 *-------------------------------------------------------------------------*/

void Test_median(CuTest* tc)
{
	double pts_odd[7] = {1.0, 7.0, 2.0, 5.0, 4.0, 3.0, 6.0}; /* median 4.0 */
	double expect_med_odd = 4.0;
	double pts_even[8] = {8.0, 2.0, 7.0, 5.0, 6.0, 4.0, 3.0, 1.0}; /* median 4.5 */
	double expect_med_even = 4.5;	
	double med_odd = median(pts_odd, 7);

	double med_even = median(pts_even, 8);
	
	char buf[80];
	sprintf(buf, "Odd num els, med_odd = %.1f expected %.1f\n", med_odd, expect_med_odd);
	CuAssert(tc, buf, abs(med_odd - expect_med_odd) < DBL_EPSILON);
	sprintf(buf, "Odd num els, med_even = %.1f expected %.1f\n", med_even, expect_med_even);
	CuAssert(tc, buf, abs(med_even - expect_med_even) < DBL_EPSILON);
	
	double nan = median(NULL, 0);
	CuAssert(tc, "NULL test", isnan(nan));
	
}

