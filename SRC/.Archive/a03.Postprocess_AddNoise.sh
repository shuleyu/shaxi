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

for count in `seq ${StartFrom} ${EndAt}`
do
	dir=${WORKDIR}/Calculation/${ModelType}_${count}
    DATADIR=${dir}/OUTPUT

	if ! [ -d ${DATADIR}/Raw ]
	then
		echo "    ==> Run SHaxi first on ${ModelType}_${count}.."
		continue
	fi

    EQname=`echo "201500000000 + ${count}" | bc `

    trap "rm -f ${WORKDIR}/*${RunNumber} ${dir}/OUTPUT/*sac ${dir}/OUTPUT/*T ${dir}/OUTPUT/tmpfile*; exit 1" SIGINT
    echo "    ==> Running SHaxi Postprocess on ${ModelType}_${count}.."

    cd ${DATADIR}
	cp Raw/*T .

    rm -f sac.macro
    for file in `ls ${ModelType}_${count}_*.T`
    do
		FixDelta=`head -n 1 ${file} | awk '{print $1-0.00001}'`

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

    RealD=`saclst evdp f 001.sac | awk '{print $2*1000}'`

	# Convolve SAC files with gaussian source of the desired dominating period.
	if [ `echo "${DominantPeriod}>0" | bc` -eq 1 ]
	then

		saclst O f `ls ???.sac` > tmpfile_$$
		${EXECDIR}/PostProcess_AddNoise.out 0 2 4 << EOF
tmpfile_$$
${SRCDIR}/noise.sac
${DominantPeriod}
1799.0
${FixDelta}
${NoiseLevel}
EOF
		if [ $? -ne 0 ]
		then
			echo "    !=> PostProcess C code failed ..."
			exit 1
		fi
	fi

    # Add headers.
    rm -f sac.macro
    echo "" > tmpfile_$$
    seq 1 `head -n 1 ${dir}/recs` >> tmpfile_$$
    paste tmpfile_$$ ${dir}/recs | awk 'NR>1 {printf "%.3d\t%.2lf\n",$1,$2}'> tmpfile1_$$

    while read code gcarc
    do
		if [ ${Integrate} -ne 0 ]
		then
			Velocity="IDISP"
		else
			Velocity="IVEL"
		fi

        cat >> sac.macro << EOF
r ${code}.sac
${Operate}
ch gcarc ${gcarc} KSTNM s${code} O 0 EVDP ${RealD} LCALDA false KNETWK ${ModelType}_${count} KCMPNM THT IDEP ${Velocity} CMPAZ 90 CMPINC 90 STLA 0.0 STLO ${gcarc} STEL 0.0 STDP 0.0 KEVNM ${EQname} EVLA 0.0 EVLO 0.0 AZ 90.0 BAZ 270.0
w ${EQname}.${ModelType}_${count}.s${code}.THT.sac
EOF
    done < tmpfile1_$$

    # Execute SAC macro.
    sac >/dev/null << EOF
m sac.macro
q
EOF

	# Clean up.
    rm -f ???.sac

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
	rm -f sac.macro
    for file in `find ${WORKDIR}/${EQname} -iname "*T.sac"`
    do
		cat >> sac.macro << EOF
cut b 0 1
r ${file}
ch KCMPNM THR
w ${file%T.sac}R.sac
ch KCMPNM THZ
w ${file%T.sac}Z.sac
EOF
    done

    sac >/dev/null << EOF
m sac.macro
q
EOF

    # Clean up.
    rm -f tmpfile* sac.macro *.T

done # Done Model loop.

cd ${CODEDIR}

exit 0
