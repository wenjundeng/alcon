TARGET := alcon

FTN := mpif90
FFLAGS := -O2

.PHONY : build run clean

build : ${TARGET}

include ${SLEPC_DIR}/conf/slepc_common_variables

${TARGET} : ${TARGET}.o
	${FTN} ${FFLAGS} -o ${TARGET} ${TARGET}.o ${SLEPC_LIB}

${TARGET}.o : ${TARGET}.F90 ${TARGET}_dec.h90 ${TARGET}_init.h90 ${TARGET}_final.h90
	${FTN} ${FFLAGS} -c -o ${TARGET}.o ${TARGET}.F90 ${FC_INCLUDES}

run : ${TARGET}
	mpiexec -n 3 ./${TARGET} 2> ${TARGET}.err | tee ${TARGET}.out

clean :
	rm -f ${TARGET} ${TARGET}.o ${TARGET}.out ${TARGET}.err *.gp.dat fort.*

