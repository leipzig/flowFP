# Simple fingerprint test.
#
#  create a model and make sure it works.
test.fingerprint1 <- function() {
	data(fs1)
	fp <- flowFP(fs1)
	checkTrue(!is.null(fp))
	checkTrue(is(fp)[1] == "flowFP")
}

# More advanced test.
#
# create a model with a name, and three parameters, check that the model got 
# created correctly.
test.FP2 <- function() {
	data(fs1)
	TEST_NAME = "bluepudding"
	fp <- flowFP(fs1, nRecursions=3, parameters=c(1,5,6), name=TEST_NAME)
	checkTrue(!is.null(fp))
	checkTrue(is(fp)[1] == "flowFP")
	checkTrue(nFeatures(fp) == 8)
	checkTrue(name(fp) == TEST_NAME)
	checkTrue(all(parameters(fp) == c("FS Lin", "FL3 Log", "FL4 Log")))
	names = c("FI05_000942_001.LMD", "FI05_000942_002.LMD", "FI05_000942_003.LMD",
	          "FI05_000942_004.LMD", "FI05_000942_005.LMD", "FI05_000942_006.LMD",
	          "FI05_000942_007.LMD")
	checkTrue(all(sampleNames(fp) == names))
	fp_counts_check =  rowSums(fp@counts)
	fs_lens = unlist(fsApply(fs1, nrow))
	checkTrue(all(fp_counts_check == fs_lens))
	checkTrue(sum(counts(fp)) == sum(fsApply(fs1, nrow)))
	
	
}

test.Reducing_nRecursions <- function() {
	data(fs1)
	TEST_NAME = "bluepudding"
	fp <- flowFP(fs1, nRecursions=10, parameters=c(1,5,6), name=TEST_NAME)
	checkTrue(!is.null(fp))
	checkTrue(is(fp)[1] == "flowFP")
	checkTrue(nRecursions(fp) == 10)
	
	checkTrue( all(rowSums(counts(fp)) == 30000) )
	nRecursions(fp) = 5
	
	checkTrue( all(rowSums(counts(fp)) == 30000) )
}

test.Reducing_model_nRecursions <- function () {
	data(fs1)
	TEST_NAME = "bluepudding"
	mod <- flowFPModel (fs1, nRecursions="auto", parameters=c(2,5), name=TEST_NAME)
	fp <- flowFP (fs1, mod)
	checkTrue (is.flowFP(fp))
	nRecursions(mod) <- 5
	fp <- flowFP (fs1, mod)
	checkTrue (is.flowFP(fp))	
}
