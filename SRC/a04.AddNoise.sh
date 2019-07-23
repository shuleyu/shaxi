#!/bin/bash

# ==============================================================
# This script run make SAC files for each models.
#
# Shule Yu
# Aug 28 2015
# ==============================================================

[ ${RunOnOffice} -ne 1 ] && echo "Not on office .." && exit 1

mkdir -p ${WORKDIR}/tmpdir_$$
cd ${WORKDIR}/tmpdir_$$
trap "rm -rf ${WORKDIR}/tmpdir_$$; exit 1" SIGINT EXIT

echo ""
echo "--> `basename $0` is running."

# ================================================
#             ! Work Begin !
# ================================================

mysql -N ${DB} > tmpfile_data_gcarc_snr << EOF
select shift_gcarc,snr_scs,evlo,evla,stlo,stla,eq,stnm,hitlo,hitla from Master_a10 order by eq,gcarc;
EOF

ls ${SRCDIR}/Noises/*sac > tmpfile_noise_filename

for count in `seq ${StartFrom} ${EndAt}`
do
    ModelEQ=`echo "201500000000 + ${count}" | bc `

	rm -f tmpfile_$$
	for synfile in `ls ${WORKDIR}/${ModelEQ}/*THT.sac`
	do
		mysql -N ${SYNDB} >> tmpfile_$$ << EOF
select file,gcarc,round((S-Begin)/delta),round((Peak_ScS+ScS-Begin)/delta) from Master_a10 where file="${synfile}";
EOF
	done

	sort -g -k 2,2 tmpfile_$$ > tmpfile_syn_filename_gcarc_amp

	FirstTrace=`awk 'NR==1 {print $1}' tmpfile_syn_filename_gcarc_amp`
	SynDelta=`saclst delta f ${FirstTrace} | awk '{print $2}'`

	${EXECDIR}/AddNoise.out 1 4 3 << EOF
${AddNoise}
tmpfile_syn_filename_gcarc_amp
tmpfile_noise_filename
tmpfile_data_gcarc_snr
tmpfile_forheader
${SynDelta}
${F1}
${F2}
EOF
	echo "${ModelEQ}: `wc -l < tmpfile_forheader` / `wc -l < tmpfile_data_gcarc_snr`"

	# Make outputs
	Cnt=1
	mkdir -p ${OutputDirectory}/${ModelEQ}
	rm -f ${OutputDirectory}/${ModelEQ}/ScS_Hitting.txt
	while read originalFile shift_gcarc evlo evla stlo stla EQ stnm hitlo hitla
	do
		echo "${EQ}_${stnm} ${hitlo} ${hitla}" >> ${OutputDirectory}/${ModelEQ}/ScS_Hitting.txt
		mkdir -p ${OutputDirectory}/${ModelEQ}/${EQ}
		netwk=`saclst knetwk f ${originalFile} | awk '{print $2}'`
		evdp=`saclst evdp f ${originalFile} | awk '{print $2}'`
		[ ${Integrate} -ne 0 ] && Velocity="IDISP" || Velocity="IVEL"

		NewFile=${EQ}.${netwk}.${stnm}.THT.sac

		sac > /dev/null << EOF
r ${Cnt}.sac
ch gcarc ${shift_gcarc} KSTNM ${stnm} O 0 EVDP ${evdp} LCALDA false KNETWK ${netwk} KCMPNM THT IDEP ${Velocity} CMPAZ 90 CMPINC 90 STLA ${stla} STLO ${stlo} STEL 0.0 STDP 0.0 KEVNM ${EQ} EVLA ${evla} EVLO ${evlo} AZ 0.0 BAZ 0.0
w over
q
EOF

		sac > /dev/null << EOF
cut O 300 1790
r ${Cnt}.sac
w ${OutputDirectory}/${ModelEQ}/${EQ}/${NewFile}
q
EOF

# 		cp ${Cnt}.Noise.sac ${OutputDirectory}/${ModelEQ}/${EQ}/${NewFile}_Noise
		cp ${originalFile%T.sac}R.sac ${OutputDirectory}/${ModelEQ}/${EQ}/${NewFile%T.sac}R.sac
		cp ${originalFile%T.sac}Z.sac ${OutputDirectory}/${ModelEQ}/${EQ}/${NewFile%T.sac}Z.sac
        taup_setsac -mod prem -ph P-0,Pdiff-0,pP-1,S-2,Sdiff-2,sS-3,PP-4,SS-5,SKKS-6,PKP-7,SKS-8,ScS-9 ${OutputDirectory}/${ModelEQ}/${EQ}/${NewFile}

		Cnt=$((Cnt+1))
	done < tmpfile_forheader

done # Done Model loop.

cd ${WORKDIR}

exit 0
