TOP=../..
PKG=${shell cd ..;pwd}
SUITE=doRUnit.R
R=R
RBLIB=R_blib
ABS_RBLIB=${TOP}/${RBLIB}
LIB_DIR=/netusr/lib64/R/penn_library/hdb/

all: inst_local test_local

inst_local: ${ABS_RBLIB} #Install package in local blib directory
	cd ${TOP};\
	${R} CMD INSTALL ${PKG} -l ${RBLIB}

test_local: # test localling installed package
	export RCMDCHECK=FALSE;\
	export R_LIBS=../../${RBLIB}:$$R_LIBS;\
	${R} --vanilla --slave < ${SUITE}

inst: # Install package
	cd ${TOP};\
	${R} CMD INSTALL_HAH ${PKG} -l ${LIB_DIR}

test: # Run unit tests
	export RCMDCHECK=FALSE;\
	${R} --vanilla --slave < ${SUITE}
	
${ABS_RBLIB}:
	-@if [ ! -d "$@" ]; then\
		mkdir "$@"; \
		echo "Created Directory -->$@";	fi
