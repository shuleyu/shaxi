#!/bin/bash

# ====================================================================
# This script plot the nc file calculated by SHaxi.
#
# Better executed while a02 is not running.
#
# Shule Yu
# Oct 19 2015
# ====================================================================

echo ""
echo "--> `basename $0` is running. "
mkdir -p ${PLOTDIR}/tmpdir_$$
cd ${PLOTDIR}/tmpdir_$$
trap "rm -rf ${PLOTDIR}/tmpdir_$$; exit 1" SIGINT EXIT

gmtset PAPER_MEDIA = letter
gmtset ANNOT_FONT_SIZE_PRIMARY = 8p
gmtset LABEL_FONT_SIZE = 9p
gmtset LABEL_OFFSET = 0.1c

# ==================================================
#              ! Work Begin !
# ==================================================

# Temtatively Run SHaxi for 120 seconds, let it create the *nc file.
for SHEXE in `find ${WORKDIR}/Calculation -iname "shaxi*.out"`
do
	dir=${SHEXE%/*}
	cd ${dir}

	ls OUTPUT/*nc > tmpfile_$$ 2>/dev/null

	if ! [ -s tmpfile_$$ ]
	then
		timeout 180 mpirun -np ${RankNum} ./shaxi_*.out < Par_*
		sleep 60
	fi

	rm -f tmpfile_$$

done

# Ploting.
cd ${PLOTDIR}/tmpdir_$$


count=0
echo "taup.distance.precision=3" > .taup
for receiver in `seq ${DISTMIN} 0.5 ${DISTMAX}`
do
	count=$((count+1))
	taup_path -mod prem -h ${EVDE} -ph ScS -deg ${receiver}
	Val=`grep 3480 taup_path.gmt | awk 'NR==1 {print $1}'`
	echo "${Val} 0" >> tmpfile_${EVDE}_symbol_$$
	if [ $((count%5)) -eq 1 ]
	then
		echo "${Val} -5 8 0 0 CT ${receiver}" >> tmpfile_${EVDE}_text_$$
		echo "${Val} -1.5 " >> tmpfile_${EVDE}_tick_$$
	fi

done

for receiver in `seq ${DISTMIN} 1 ${DISTMAX}`
do

	${EXECDIR}/PlotPath.out 0 1 3 << EOF
tmpfile_${EVDE}_path_${receiver}_$$
${EVDE}
${receiver}
`echo "${RMAX}" | awk '{print 2891-$1}'`
EOF

    if [ $? -ne 0 ]
    then
        echo "    !=> PlotPath C code failed on Gcarc: ${receiver} ..."
        exit 1
    fi

done

for count in `seq ${StartFrom_Plot} ${EndAt_Plot}`
do
	Model=${ModelType}_${count}
	echo "    ==> Plotting model ${Model}..."
	recfile=${WORKDIR}/Calculation/${Model}/recs
	OUTFILE=${count}_1.ps
	OUTFILE_Zoom=${count}_2.ps
    EQname=`echo "201500000000 + ${count}" | bc `
	EPulvz=`grep ${EQname} ${WORKDIR}/Calculation/index | awk '{print $3}'`
	${BASHCODEDIR}/Findfield.sh ${WORKDIR}/Calculation/index "<EQ> <Vs> <Rho>" > tmpfile_$$
	ColorScale_Max1=`grep ${EQname} tmpfile_$$ | awk '{print (1-$2)*100}'`
	ColorScale_Max2=`grep ${EQname} tmpfile_$$ | awk '{print (1-$3)*100}'`

	# Process the nc file.
	for NcFile in `find ${WORKDIR}/Calculation/${Model} -iname "*nc"`
	do

		ncdump ${NcFile} > all.ncdump
		${SRCDIR}/preprocess.sh

		Tsize=`awk 'BEGIN {FS=","} {if (NR==1) print $2-$1}' Theta.dat`
		Rsize=`awk 'BEGIN {FS=","} {if (NR==1) print $1-$2}' R.dat`
		Rank=${NcFile%_model.nc}
		Rank=${Rank#*rank}

	# C code.
	${EXECDIR}/CreateGrid.out 2 8 0 << EOF
`fgrep -o , Theta.dat  | wc -l`
`fgrep -o , R.dat  | wc -l`
Theta.dat
R.dat
Vs.dat
Rho.dat
Vs_${Rank}.grid
Rho_${Rank}.grid
Vs_${Rank}_zoom.grid
Rho_${Rank}_zoom.grid
EOF

	if [ $? -ne 0 ]
	then
		echo "    !=> CreateGrid C code failed on ${Model} / rank ${Rank}..."
		exit 1
	fi

	done

	NRank=`find ${WORKDIR}/Calculation/${Model} -iname "*nc" | wc -l`

	# Plotting Zoom.

	# Plot scale.

	Zoom_Min=`echo "${EPulvz} ${HalfLateralSize}"| awk '{print $1-$2}'`
	Zoom_Max=`echo "${EPulvz} ${HalfLateralSize}"| awk '{print $1+$2}'`

	REG_Z="-R${Zoom_Min}/${Zoom_Max}/0/${RMAX}"
	PROJ="-JX9i/2.5i"

	# some texts.
    pstext -JX11.0i/1i -R-1/1/-1/1 -X0i -Y7.5i -N -K > ${OUTFILE_Zoom} << EOF
0 0 20 0 0 CB Model: ${Model}. Reference: PREM
EOF
	psxy -J -R -X1i -Y0.4i -O -K >> ${OUTFILE_Zoom} << EOF
EOF

	for name in Vs Rho
	do
		psxy -J -R -Y-3.1i -O -K >> ${OUTFILE_Zoom} << EOF
EOF

		if [ ${name} = Rho ]
		then
			legend="@~\162@~"
		else
			legend="${name}"
		fi

		# Color scale.
		if [ "${ColorScale}" -ne 0 ]
		then
			scale=${ColorScale}
		else
			if [ ${name} = Rho ]
			then
				scale=${ColorScale_Max2}
			else
				scale=${ColorScale_Max1}
			fi
		fi
		if [ `echo "${scale}<1"|bc` -eq 1 ]
		then
			scale=1
		fi

		# Make cpt
		makecpt -Cpolar -I -T-${scale}/${scale}/0.5 -Z > tmp.cpt

		# Make grid file and plot.
		for Rank in `seq 0 $((NRank-1))`
		do
			num=`printf "%.3d" ${Rank}`
			REG="-R`echo "${ThetaMax}/${NRank}*${Rank}"|bc -l`/`echo "${ThetaMax}/${NRank}*(${Rank}+1)"|bc -l`/0/${RMAX}"
			xyz2grd ${name}_${num}_zoom.grid ${REG} -I${Tsize}/${Rsize} -G${name}.grd
			grdimage -Q ${name}.grd ${REG_Z} ${PROJ} -Ctmp.cpt -O -K >> ${OUTFILE_Zoom}
		done

		psbasemap ${REG_Z} ${PROJ} -Ba2f1:"@~\104@~(@~\260@~)":/a20f10:"Height above CMB (km)":WNe -O -K >> ${OUTFILE_Zoom}
		psxy -J -R -W1p -O -K >> ${OUTFILE_Zoom} << EOF
${Zoom_Min} 0
${Zoom_Max} 0
EOF
		awk -v Min=${Zoom_Min} -v Max=${Zoom_Max} '{if ($1>Min && $1<Max) print $0}' tmpfile_${EVDE}_tick_$$ | psxy -J -R -Sy0.07i -W1p,purple -N -O -K >> ${OUTFILE_Zoom}

		# Plot ray paths.
		for file in `ls tmpfile_*path*_$$`
		do
			psxy ${file} -J -R -Wthin -O -K >> ${OUTFILE_Zoom}
		done

		awk -v Min=${Zoom_Min} -v Max=${Zoom_Max} '{if ($1>Min && $1<Max) print $0}' tmpfile_${EVDE}_symbol_$$ | psxy -J -R -Sc0.05i -Ggreen -Wblack -N -O -K >> ${OUTFILE_Zoom}
		awk -v Min=${Zoom_Min} -v Max=${Zoom_Max} '{if ($1>Min && $1<Max) print $0}' tmpfile_${EVDE}_text_$$ | pstext -J -R -N -O -K >> ${OUTFILE_Zoom}

		## plot scale bar.
		B=`echo ${scale} | awk '{if ($1>2) print $1/5; else print "0.2"}'`
		psscale -Ctmp.cpt -D4.5i/0.22i/3.0i/0.13ih -B${B}/:"@~\144@~${legend} (%)": -Y-0.7i -N300 -O -K >> ${OUTFILE_Zoom}

	done

	psxy -J -R -O >> ${OUTFILE_Zoom} << EOF
EOF

	ps2pdf ${OUTFILE_Zoom} ${PLOTDIR}/${Model}_Zoom.pdf

	# Plotting Global.

	# Plot scale.
	size=`echo ${ThetaMax} | awk '{if ($1>=90) print "0.0017";else if ($1>=30) print 0.0017/sin($1/180*3.14159) ; else printf "0.0034"}'`
	PROJ="-Jpa${size}/`echo "${ThetaMax}/2" | bc -l`z"
	REG_G="-R0/${ThetaMax}/3480/6371"

	# some texts.
    pstext -JX8.5i/1i -R-1/1/-1/1 -X0i -Y9.5i -P -N -K > ${OUTFILE} << EOF
0 0 20 0 0 CB Model: ${Model}. Reference: PREM
EOF

	psxy -J -R -X1i -O -K >> ${OUTFILE} << EOF
EOF

	for name in Vs Rho
	do
		psxy -J -R -Y-3.5i -O -K >> ${OUTFILE} << EOF
EOF

		if [ ${name} = Rho ]
		then
			legend="@~\162@~"
		else
			legend="${name}"
		fi

		# Color scale.
		if [ "${ColorScale}" -ne 0 ]
		then
			scale=${ColorScale}
		else
			if [ ${name} = Rho ]
			then
				scale=${ColorScale_Max2}
			else
				scale=${ColorScale_Max1}
			fi
		fi
		if [ `echo "${scale}<1"|bc` -eq 1 ]
		then
			scale=1
		fi

		# Make cpt
		makecpt -Cpolar -I -T-${scale}/${scale}/0.5 -Z > tmp.cpt

		# Make grid file and plot.
		for Rank in `seq 0 $((NRank-1))`
		do
			num=`printf "%.3d" ${Rank}`
			REG="-R`echo "${ThetaMax}/${NRank}*${Rank}"|bc -l`/`echo "${ThetaMax}/${NRank}*(${Rank}+1)"|bc -l`/3480/6371"
			xyz2grd ${name}_${num}.grid ${REG} -I${Tsize}/${Rsize} -G${name}.grd
			grdimage -Q ${name}.grd ${REG_G} ${PROJ} -Ctmp.cpt -E400 -O -K >> ${OUTFILE}
		done


		psbasemap ${REG_G} ${PROJ} -Ba10f2/a1000f200NWe -O -K >> ${OUTFILE}

		## source.
		psxy ${REG_G} ${PROJ} -Sa0.2i -Gyellow -N -O -K >> ${OUTFILE} << EOF
0 `echo "6371-${EVDE}" | bc -l`
EOF

		## stations.
		awk 'NR>1 {print $1,"6371"}' ${recfile} > tmpfile_$$
		psxy tmpfile_$$ ${REG_G} ${PROJ} -Si0.05i -Gblue -N -O -K >> ${OUTFILE}

		## plot scale bar.
		B=`echo ${scale} | awk '{if ($1>2) print $1/5; else print "0.2"}'`
		psscale -Ctmp.cpt -D3.45i/0.2i/3.0i/0.13ih -B${B}/:"@~\144@~${legend} (%)": -Y-0.7i -N300 -O -K >> ${OUTFILE}

	done

	psxy -J -R -O >> ${OUTFILE} << EOF
EOF

	ps2pdf ${OUTFILE} ${PLOTDIR}/${Model}.pdf

done # Done model loop.

cat `ls *_2.ps | sort -n` > Zoom.ps
cat `ls *_1.ps | sort -n` > Global.ps
ps2pdf Zoom.ps ${PLOTDIR}/Zoom.pdf
ps2pdf Global.ps ${PLOTDIR}/Global.pdf

cd ${CODEDIR}

exit 0
