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

void Test_create_idx_lut(CuTest* tc)
{
	int data[6] = {1, 2, 1, 1, 2, 2};
	CIVECTOR tags;
	tags.data = data;
	tags.len = 6;
	
	int expect_0[3] = {0, 2, 3};
	int expect_1[3] = {1, 4, 5};

	CIVECTOR **idx_lut = create_idx_lut(&tags, 2);
	CuAssert(tc, "NULL", idx_lut != NULL);
	CuAssert(tc, "Length[0]", idx_lut[0]->len == 3);
	CuAssert(tc, "Length[1]", idx_lut[1]->len == 3);
	int i;
	char buf[80];
	for(i=0; i < 3; i++) {
		sprintf (buf, "idx_lut[0]->data[%d] = %d, expected %d", i, idx_lut[0]->data[i], expect_0[i]);
		CuAssert(tc, buf, idx_lut[0]->data[i] == expect_0[i]);
		sprintf (buf, "idx_lut[1]->data[%d] = %d, expected %d", i, idx_lut[1]->data[i], expect_1[i]);
		CuAssert(tc, buf, idx_lut[0]->data[i] == expect_0[i]);
	}
	
	
	
}
