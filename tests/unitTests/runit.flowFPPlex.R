# Simple test.
#
#  create a fingerprint plex and make sure it works.
test.plex1 <- function() {
	data(fs1)
	fp <- flowFP(fs1)
	plex = flowFPPlex(c(fp, fp))
	checkTrue(!is.null(plex))
	checkTrue(is(plex) == "flowFPPlex")
}

# More advanced test.
#
#  create a fingerprint plex and make sure it works.
test.plex2 <- function() {
	data(fs1)
	fp <- flowFP(fs1)
	plex = flowFPPlex(c(fp, fp))
	checkTrue(!is.null(plex))
	checkTrue(is(plex) == "flowFPPlex")
}


test.plex_append <- function() {
	data(fs1)
	fp <- flowFP(fs1)
	plex = flowFPPlex()
	checkTrue(!is.null(plex))
	checkTrue(is(plex) == "flowFPPlex")
	plex = append(plex, fp)
	checkTrue(is.flowFPPlex(plex))
}


test.plex_testClasses <- function() {
	data(fs1)
	fp1 <- flowFP(fs1,sampleClasses=c(rep("A", 3), rep("B", 4)))
	plex = flowFPPlex(c(fp1, fp1))
	checkTrue(!is.null(plex))
	checkTrue(is.flowFPPlex(plex))
}

test.plex_append_plex <- function() {
	data(fs1)
	fp <- flowFP(fs1)
	plex = flowFPPlex()
	checkTrue(!is.null(plex))
	checkTrue(is(plex) == "flowFPPlex")
	plex = append(plex, fp)
	checkTrue(is.flowFPPlex(plex))
	plex2 = append(plex, c(fp, fp))
	checkTrue(is.flowFPPlex(plex2))
	plex3 = append(plex, plex2)
	checkTrue(is.flowFPPlex(plex3))
}
