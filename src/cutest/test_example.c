#include <assert.h>
#include <setjmp.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "CuTest.h"

/*-------------------------------------------------------------------------*
 * Test
 *-------------------------------------------------------------------------*/

void Test_example(CuTest* tc)
{
	char buf[80];
	
	sprintf(buf, "There is a big problem because %d does NOT equal %d\n", 1, 1);
	CuAssert(tc, buf, 1 == 1);
}

