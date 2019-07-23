#!/bin/bash

# ==============================================================
# This script run make SAC files for each models.
#
# Shule Yu
# Aug 28 2015
# ==============================================================

[ ${RunOnOffice} -ne 1 ] && echo "Not on office .." && exit 1
echo ""
echo "--> `basename $0` is running."

# ================================================
#             ! Work Begin !
# ================================================

if [ ${Integrate} -ne 0 ]
then
	Operate="int"
else
	Operate="mul 1"
fi

for dir in `ls -d ${WORKDIR}/Calculation/${ModelType}_*`
do
    DATADIR=${dir}/OUTPUT
    count=${dir##*_}
    EQname=`echo "201500000000 + ${count}" | bc `

    trap "rm -f ${dir}/OUTPUT/*T ${dir}/OUTPUT/tmpfile*; exit 1" SIGINT
    echo "    ==> Running SHaxi Postprocess on ${ModelType}_${count}.."

    cd ${DATADIR}
	cp Raw/*T .

    rm -f sac.macro
    for file in `ls ${ModelType}_${count}_*.T`
    do
		## proper normalization
		awk '{ if (NR == 31) print ("   1.000000","     ", $2, "     ", $3, "     ", $4, "     ", $5); else print ($0)}' ${file} > temp
		mv temp ${file}

		## convert outputs to SAC.
        outfile=${file##*_}
        outfile=${outfile%.T}.sac
        cat >> sac.macro << EOF
r alpha ${file}
w ${outfile}
EOF
	done

    # Execute SAC macro.
    sac >/dev/null << EOF
m sac.macro
q
EOF

	# next we need to convolve each record with a gaussian of the desired dominating period
	if [ `echo "${DominantPeriod}>0" | bc` -eq 1 ]
	then
		rm -f sac.macro
		for file in `ls ???.sac`
		do
			${EXECDIR}/convolve.out ${file} -d gauss -t ${DominantPeriod} > /dev/null
			cat >> sac.macro << EOF
r ${file}
cut 50 &1,E
r
${Operate}
w over
EOF
		done

		# Execute SAC macro.
		sac >/dev/null << EOF
m sac.macro
q
EOF
	fi

    # Add gcarc to the header.
    rm -f sac.macro
    echo "" > tmpfile_$$
    seq 1 `head -n 1 ${dir}/recs` >> tmpfile_$$
    paste tmpfile_$$ ${dir}/recs | awk 'NR>1 {printf "%.3d\t%.2lf\n",$1,$2}'> tmpfile1_$$

    while read code gcarc
    do
        cat >> sac.macro << EOF
r ${code}.sac
ch gcarc ${gcarc}
w ${EQname}.${ModelType}_${count}.${code}.THT.sac
EOF
    done < tmpfile1_$$

    # Execute SAC macro.
    cat >> sac.macro << EOF
q
EOF
    sac >/dev/null << EOF
m sac.macro
EOF
    RealD=`saclst evdp f 001.sac | awk '{print $2*1000}'`
    rm -f ???.sac

    # Set Omarker and Depth and other headers.
    rm -f sac.macro
    for file in `ls ${EQname}*.sac`
    do
        COMP=${file%%.sac}
        COMP=${COMP##*.}
        cat >> sac.macro << EOF
r ${file}
ch O 0 evdp ${RealD} LCALDA false KNETWK ${ModelType} KCMPNM ${COMP}
w over
EOF
    done

    # Execute SAC macro.
    sac >/dev/null << EOF
m sac.macro
q
EOF


    # Set time header.
    for file in `ls ${EQname}*.sac`
    do
        taup_setsac -mod prem -ph P-0,Pdiff-0,pP-1,S-2,Sdiff-2,sS-3,PP-4,SS-5,SKKS-6,PKP-7,SKS-8,ScS-9 ${file}
    done

    # Move SAC to OUTPUT dir.
    mkdir -p ${WORKDIR}/${EQname}
    rm -f ${WORKDIR}/${EQname}/*
    mv ${EQname}*.sac ${WORKDIR}/${EQname}

    # Creat false components.
    for file in `find ${WORKDIR}/${EQname} -iname "*T.sac"`
    do
        cp ${file} ${file%T.sac}Z.sac
        cp ${file} ${file%T.sac}R.sac
    done

    # Clean up.
    rm -f tmpfile* sac.macro *.T

done # Done Model loop.

cd ${CODEDIR}

exit 0
