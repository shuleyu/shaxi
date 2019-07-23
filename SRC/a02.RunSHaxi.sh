#!/bin/bash

# ==============================================================
# This script run SHaxi and make SAC files for each models.
#
# Shule Yu
# Aug 28 2015
# ==============================================================

echo ""
echo "--> `basename $0` is running."
cd ${WORKDIR}

# ================================================
#             ! Work Begin !
# ================================================

for count in `seq ${BeginRun} ${EndRun}`
do
	dir=${WORKDIR}/Calculation/${ModelType}_${count}
	mkdir -p ${dir}
    DATADIR=${dir}/OUTPUT

	if [ ${RunOnOffice} -eq 1 ]
	then
		CreateDir="mkdir -p ${dir}/OUTPUT"
		Command="${MPIRUN} -np ${RankNum} ./shaxi_${count}.out < Par_${ModelType}"
		MoveData=""
		DeleteDir=""
	else 
        read RunN < ${WORKDIR}/Calculation/${ModelType}_${count}/RunNumber.txt
		CreateDir="ssh shule@${HostName} 'mkdir -p /local/shule/${ModelType}_${RunN}_${count}/OUTPUT'"
		Command="${MPIRUN} --host ${HostName} -np ${RankNum} ./shaxi_${count}.out < Par_${ModelType}"
		MoveData="ssh shule@${HostName} 'mv /local/shule/${ModelType}_${RunN}_${count}/OUTPUT ${dir}/'"
		DeleteDir="ssh shule@${HostName} 'rm -r /local/shule/${ModelType}_${RunN}_${count}'"
	fi


	rm -rf ${dir}/OUTPUT
	cat > tmpfile_${count}_$$ << EOF
#!/bin/bash

echo "    ==> Running SHaxi on ${ModelType}_${count}.."

cd ${dir}

${CreateDir}
${Command}
${MoveData}
${DeleteDir}

mkdir -p ${DATADIR}/Raw
mv ${DATADIR}/*T ${DATADIR}/Raw

exit 0

EOF

# if [ \$? -ne 0 ]
# then
# 	echo "    !=> SHaxi Abort on ${ModelType}_${count} !"
# 	touch ${dir}/ERROR
# 	chmod +x ${dir}/ERROR
# else
# 	rm -f ${dir}/ERROR
# fi

	Running=`jobs|grep Running | wc -l`
	while [ $((Running*RankNum)) -eq ${TotalRank} ]
	do
		sleep 60
		Running=`jobs|grep Running | wc -l`
	done

	bash ./tmpfile_${count}_$$ &
	sleep 180

done # Done Model loop.

Running=`jobs|grep Running | wc -l`
while [ $((Running*RankNum)) -ne 0 ]
do
	sleep 60
	Running=`jobs|grep Running | wc -l`
done

rm ${WORKDIR}/tmpfile_*$$ 
cd ${CODEDIR}

exit 0
