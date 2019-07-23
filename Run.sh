#!/bin/bash

# ==============================================================
# This script use the SHaxi method to compute synthesis. Since
# the original program is easy to use, this script is dedicated
# to making ULVZ models.
#
# Shule Yu
# Aug 27 2015
# ==============================================================

# Export variables to all sub scripts.
set -a
CODEDIR=${PWD}
SRCDIR=${CODEDIR}/SRC
RunNumber=$$

#============================================
#            ! Test Files !
#============================================
if ! [ -e ${CODEDIR}/INFILE ]
then
    echo "INFILE not found ..."
    exit 1
fi

#============================================
#            ! Parameters !
#============================================

# DIRs.
WORKDIR=`grep "<WORKDIR>" ${CODEDIR}/INFILE | awk '{print $2}'`
ModelType=`grep "<ModelType>" ${CODEDIR}/INFILE | awk '{print $2}'`
WORKDIR="${WORKDIR}.${ModelType}"


trap "rm -f ${WORKDIR}/tmpfile*_$$; exit 1" SIGINT
mkdir -p ${WORKDIR}/LIST
mkdir -p ${WORKDIR}/INPUT
cp ${CODEDIR}/INFILE ${WORKDIR}/tmpfile_INFILE_$$
cp ${CODEDIR}/INFILE ${WORKDIR}/INPUT/INFILE_`date +%m%d_%H%M`
cp ${CODEDIR}/LIST.sh ${WORKDIR}/tmpfile_LIST_$$
cp ${CODEDIR}/LIST.sh ${WORKDIR}/LIST/LIST_`date +%m%d_%H%M`
chmod -x ${WORKDIR}/LIST/*
cd ${WORKDIR}

# Deal with single parameters.
grep -n "<" ${WORKDIR}/tmpfile_INFILE_$$     \
| grep ">" | grep -v "BEGIN" | grep -v "END" \
| awk 'BEGIN {FS="<"} {print $2}'            \
| awk 'BEGIN {FS=">"} {print $1,$2}' > tmpfile_$$
awk '{print $1}' tmpfile_$$ > tmpfile1_$$
awk '{$1="";print "\""$0"\""}' tmpfile_$$ > tmpfile2_$$
sed 's/\"[[:blank:]]/\"/' tmpfile2_$$ > tmpfile3_$$
paste -d= tmpfile1_$$ tmpfile3_$$ > tmpfile_$$
source ${WORKDIR}/tmpfile_$$

WORKDIR="${WORKDIR}.${ModelType}"

# Deal with multiple parameters.
# They are between <XXX_BEGIN> and <XXX_END>
# The list is put into ${WORKDIR}/tmpfile_XXX_${RunNumber}
grep -n "<" ${WORKDIR}/tmpfile_INFILE_$$ \
| grep ">" | grep "_BEGIN"               \
| awk 'BEGIN {FS=":<"} {print $2,$1}'    \
| awk 'BEGIN {FS="[> ]"} {print $1,$NF}' \
| sed 's/_BEGIN//g'                      \
| sort -g -k 2,2 > tmpfile1_$$

grep -n "<" ${WORKDIR}/tmpfile_INFILE_$$ \
| grep ">" | grep "_END"                 \
| awk 'BEGIN {FS=":<"} {print $2,$1}'    \
| awk 'BEGIN {FS="[> ]"} {print $1,$NF}' \
| sed 's/_END//g'                        \
| sort -g -k 2,2 > tmpfile2_$$

paste tmpfile1_$$ tmpfile2_$$ | awk '{print $1,$2,$4}' > tmpfile_parameters_$$

while read Name line1 line2
do
    Name=${Name%_*}
    awk -v N1=${line1} -v N2=${line2} '{ if ( $1!="" && N1<NR && NR<N2 ) print $0}' ${WORKDIR}/tmpfile_INFILE_$$ \
	| sed 's/^[[:blank:]]*//g' > ${WORKDIR}/tmpfile_${Name}_$$
done < tmpfile_parameters_$$

# Additional DIRs.
EXECDIR=${WORKDIR}/bin
PLOTDIR=${WORKDIR}/PLOTS

#============================================
#            ! Test Dependencies !
#============================================
CommandList="${CCOMP} ${FCOMP}"
for Command in ${CommandList}
do
    command -v ${Command} >/dev/null 2>&1 || { echo >&2 "Command ${Command} is not found. Exiting ... "; exit 1; }
done

#============================================
#            ! Compile !
#============================================
mkdir -p ${EXECDIR}
trap "rm -f ${EXECDIR}/*.o ${EXECDIR}/*.mod ${WORKDIR}/*_$$; exit 1" SIGINT

INCLUDEDIR="-I${SACDIR}/include -I${CCODEDIR} -I${CPPCODEDIR} -I/opt/gmt-5.4.3/include -I${NETCDFINCLUDEDIR}"
LIBRARYDIR="-L. -L${SACDIR}/lib -L${CCODEDIR} -L/opt/gmt-5.4.3/lib64 -L${NETCDFLIBDIR}"
LIBRARIES="-lASU_tools -lsac -lnetcdf -l${NETCDFLIB} -lgmt -lfftw3 -lsacio -lm"
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/opt/gmt-5.4.3/lib64

# ASU_tools Functions.
cd ${CCODEDIR}
make
cd ${EXECDIR}

# Fortran codes.
# ${FCOMP} -c ${SRCDIR}/mod_sac_io.f90 ${SRCDIR}/mod_wavelets.f90 ${SRCDIR}/convolve.f90
# ${FCOMP} -o convolve.out convolve.o mod_sac_io.o mod_wavelets.o


# C codes.
${CCOMP} -Wall -o ${EXECDIR}/MixThem.out ${SRCDIR}/MixThem.c ${INCLUDEDIR} ${LIBRARYDIR} -lASU_tools -lm
${CCOMP} -Wall -o ${EXECDIR}/MoreDist.out ${SRCDIR}/MoreDist.c ${INCLUDEDIR} ${LIBRARYDIR} -lASU_tools -lm
${CPPCOMP} ${CPPFLAG} -Wall -o ${EXECDIR}/MixModelInput.out ${SRCDIR}/MixModelInput.cpp ${INCLUDEDIR} ${LIBRARYDIR} -lASU_tools -lm

if [ ${RunOnOffice} -eq 1 ]
then
	for code in `ls ${SRCDIR}/*.c 2>/dev/null | grep -v fun.c`
	do
		name=`basename ${code}`
		name=${name%.c}

		${CCOMP} -Wall -o ${EXECDIR}/${name}.out ${code} ${INCLUDEDIR} ${LIBRARYDIR} ${LIBRARIES}

		if [ $? -ne 0 ]
		then
			echo "${name} C code is not compiled ..."
			rm -f ${EXECDIR}/*.o ${WORKDIR}/*_$$
			exit 1
		fi
	done

	for code in `ls ${SRCDIR}/*.cpp 2>/dev/null | grep -v fun.cpp`
	do
		name=`basename ${code}`
		name=${name%.cpp}

		${CPPCOMP} ${CPPFLAG} -Wall -o ${EXECDIR}/${name}.out ${code} ${INCLUDEDIR} ${LIBRARYDIR} ${LIBRARIES}

		if [ $? -ne 0 ]
		then
			echo "${name} C++ code is not compiled ..."
			rm -f ${EXECDIR}/*.o ${WORKDIR}/*_$$
			exit 1
		fi
	done
fi

rm -f ${EXECDIR}/*.o ${EXECDIR}/*.mod

# ==============================================
#           ! Work Begin !
# ==============================================

cat >> ${WORKDIR}/stdout << EOF

======================================
Run Date: `date`
EOF

bash ${WORKDIR}/tmpfile_LIST_$$ >> ${WORKDIR}/stdout 2>&1

cat >> ${WORKDIR}/stdout << EOF

End Date: `date`
======================================
EOF

# Clean up.
rm -f ${WORKDIR}/*_$$

exit 0
