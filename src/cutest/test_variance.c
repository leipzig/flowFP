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

void Test_variance(CuTest* tc)
{
	double pts[5] = {1.0, 2.0, 3.0, 4.0, 5.0}; /* variance 2.5 */
	double expect_var = 2.5;
	double var = variance(pts, 5);

	char buf[80];
	sprintf(buf, "Variance, var = %.1f expected %.1f\n", var, expect_var);
	CuAssert(tc, buf, abs(var - expect_var) < DBL_EPSILON);
	
	double nan = variance(NULL, 0);
	CuAssert(tc, "NULL test", isnan(nan));
	
}
