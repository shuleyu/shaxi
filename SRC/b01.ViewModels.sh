#!/bin/bash

# ====================================================================
# This script plot the nc file calculated by SHaxi.
#
# Better executed while a02 is not running.
#
# Shule Yu
# Oct 19 2015
# ====================================================================

[ ${RunOnOffice} -ne 1 ] && echo "Not on office .." && exit 1
echo ""
echo "--> `basename $0` is running. "
mkdir -p ${PLOTDIR}/tmpdir_$$
cd ${PLOTDIR}/tmpdir_$$
trap "rm -rf ${PLOTDIR}/tmpdir_$$; exit 1" SIGINT EXIT

gmt gmtset PS_MEDIA letter
gmt gmtset FONT_ANNOT_PRIMARY 8p
gmt gmtset FONT_LABEL 9p
gmt gmtset MAP_LABEL_OFFSET 0.1c
gmt gmtset MAP_FRAME_PEN 1.0p,black
gmt gmtset MAP_GRID_PEN_PRIMARY 0.2p,gray,-

# ==================================================
#              ! Work Begin !
# ==================================================

# If we don't have nc files for each model,
# temtatively run SHaxi for 120 seconds, let it create the *nc files.
for count in `seq ${StartFrom_Plot} ${EndAt_Plot}`
do
	dir=${WORKDIR}/Calculation/${ModelType}_${count}
	cd ${dir}

	ls OUTPUT/*nc > tmpfile_$$ 2>/dev/null

	if ! [ -s tmpfile_$$ ]
	then
		echo "    ==> Temptatively runing on model ${ModelType}_${count}..."
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
		echo "${Val} -5 ${receiver}" >> tmpfile_${EVDE}_text_$$
		echo "${Val} -1.5 " >> tmpfile_${EVDE}_tick_$$
	fi

done

for receiver in `seq ${DISTMIN} 1 ${DISTMAX}`
do

	${EXECDIR}/PlotPath.out 0 1 3 << EOF
tmpfile_${EVDE}_path_${receiver}_$$
${EVDE}
${receiver}
`echo "${MaxHeight}" | awk '{print 2891-$1-10}'`
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

	# 1. Plotting Zoom.

	# Plot scale.

	REG_Z="-R${Zoom_Min}/${Zoom_Max}/0/${MaxHeight}"
	PROJ="-JX9i/2.5i"

	# some texts.
    rm -f tmpfile_$$
    x=0.5
    while read line
    do
        echo "0 ${x} ${line}" >> tmpfile_$$
        x=`echo ${x} | awk '{print $1-0.3}'`
    done < ${WORKDIR}/indices/`echo 201500000000 + ${count} | bc`
# 0 0 Model: ${Model} ( @~\144@~Vs=-${ColorScale_Max1}%, @~\144@~@~\162@~=-${ColorScale_Max2}% )
    gmt pstext tmpfile_$$ -JX11.0i/1i -R-1/1/-1/1 -F+jCB+f12p -X0i -Y7.5i -N -K > ${OUTFILE_Zoom}

	gmt psxy -J -R -X1i -Y0.4i -O -K >> ${OUTFILE_Zoom} << EOF
EOF

	for name in Vs Q Rho
	do
		gmt psxy -J -R -Y-3.1i -O -K >> ${OUTFILE_Zoom} << EOF
EOF

		if [ ${name} = Rho ]
		then
			legend="@~\162@~"
			Mark="Rho"
            c=${ColorScale}
		elif [ ${name} = Vs ]
        then
			legend="${name}"
			Mark="Vs1"
            c=${ColorScale}
        else
			legend="${name}"
			Mark="Q"
            c=500
		fi

		# Plot grid.
		gmt makecpt -Cpolar -D -I -T-${c}/${c}/0.5 -Z > tmp.cpt


		find ${WORKDIR}/Calculation/${Model} -iname "*nc" > tmpfile_FileNames_$$
		${EXECDIR}/AddGrid.out 0 7 2 << EOF
${PROJ}
${REG_Z}
${Mark}
tmpfile_FileNames_$$
tmp.cpt
${OUTFILE_Zoom}
tmpfile_Ref_${name}_$$
${Zoom_Min}
${Zoom_Max}
EOF



		gmt psbasemap ${REG_Z} ${PROJ} -Ba2f1:"@~\104@~(@~\260@~)":/a20f10:"Height above CMB (km)":WNe -O -K >> ${OUTFILE_Zoom}
		gmt psxy -J -R -W1p -O -K >> ${OUTFILE_Zoom} << EOF
${Zoom_Min} 0
${Zoom_Max} 0
EOF
		awk -v Min=${Zoom_Min} -v Max=${Zoom_Max} '{if ($1>Min && $1<Max) print $0}' tmpfile_${EVDE}_tick_$$ | gmt psxy -J -R -Sy0.07i -W1p,purple -N -O -K >> ${OUTFILE_Zoom}

		# Plot ray paths. (This is using the EVDE specified in the INFILE, if EVDE is changing, these ray paths are wrong: comment out next lines.)
# 		for file in `ls tmpfile_*path*_$$`
# 		do
# 			gmt psxy ${file} -J -R -Wthin -O -K >> ${OUTFILE_Zoom}
# 		done

		awk -v Min=${Zoom_Min} -v Max=${Zoom_Max} '{if ($1>Min && $1<Max) print $0}' tmpfile_${EVDE}_symbol_$$ | gmt psxy -J -R -Sc0.05i -Ggreen -Wblack -N -O -K >> ${OUTFILE_Zoom}
		awk -v Min=${Zoom_Min} -v Max=${Zoom_Max} '{if ($1>Min && $1<Max) print $0}' tmpfile_${EVDE}_text_$$ > tmpfile_$$
		gmt pstext tmpfile_$$ -F+f8+jCT -J -R -N -O -K >> ${OUTFILE_Zoom}

		## plot scale bar.
		B=`echo ${c} | awk '{if ($1>2) print $1/5; else print "0.2"}'`
		gmt psscale -Ctmp.cpt -D4.5i/0.18i/3.0i/0.13ih -B${B}/:"@~\144@~${legend} (%)": -Y-0.7i -N300 -O -K >> ${OUTFILE_Zoom}

	done # Done rho and Vs loop.


	gmt psxy -J -R -O >> ${OUTFILE_Zoom} << EOF
EOF

	ps2pdf ${OUTFILE_Zoom} ${PLOTDIR}/${Model}_Zoom.pdf

	# A reference plot comparing with original prem. (the value is get at theta=0)
	OUTFILE_Ref=prem.ps
	cat > tmpfile_$$  << EOF
0 0 Model: ${Model} (1D refernce model,green)
EOF
    gmt pstext tmpfile_$$ -JX8.5i/1i -R-1/1/-1/1 -F+jCB+f20p -X0i -Y10i -N -P -K > ${OUTFILE_Ref}
	gmt psbasemap -JX6.5i/-9i -R0/13/0/2891 -Ba5f1g1:"Density(purple), Vs(red)":/a500f100g100:"Depth(km)":WSne -Xc -Yc -O -K >> ${OUTFILE_Ref}
	awk '{print $2,$1}' tmpfile_Ref_Vs_$$ | gmt psxy -J -R -W2p,red -O -K >> ${OUTFILE_Ref}
	awk '{print $3,$1}' tmpfile_Ref_Vs_$$ | gmt psxy -J -R -W1p,green,- -O -K >> ${OUTFILE_Ref}
	awk '{print $2,$1}' tmpfile_Ref_Rho_$$ | gmt psxy -J -R -W2p,purple -O -K >> ${OUTFILE_Ref}
	awk '{print $3,$1}' tmpfile_Ref_Rho_$$ | gmt psxy -J -R -W1p,green,- -O >> ${OUTFILE_Ref}

	ps2pdf ${OUTFILE_Ref} ${PLOTDIR}/${Model}_Ref.pdf

	continue

	# 2. Plotting Global.

	# Plot scale.
	size=`echo ${ThetaMax} | awk '{if ($1>=90) print "0.0017";else if ($1>=30) print 0.0017/sin($1/180*3.14159) ; else printf "0.0034"}'`
	PROJ="-Jpa${size}/`echo "${ThetaMax}/2" | bc -l`z"
	REG_G="-R0/${ThetaMax}/3480/6371"

	# some texts.

    rm -f tmpfile_$$
    while read line
    do
        echo "0 0 ${line}" >> tmpfile_$$
    done < ${WORKDIR}/indices/`echo 201500000000 + ${count} | bc`
# 0 0 Model: ${Model} ( @~\144@~Vs=-${ColorScale_Max1}% )
    gmt pstext tmpfile_$$ -JX8.5i/1i -R-1/1/-1/1 -F+jCB+f20p -X0i -Y9.5i -P -N -K > ${OUTFILE}

	gmt psxy -J -R -X1i -O -K >> ${OUTFILE} << EOF
EOF

	for name in Vs Q Rho
	do
		gmt psxy -J -R -Y-3.5i -O -K >> ${OUTFILE} << EOF
EOF

		if [ ${name} = Rho ]
		then
			legend="@~\162@~"
			Mark="Rho"
            c=${ColorScale}
		elif [ ${name} = Vs ]
        then
			legend="${name}"
			Mark="Vs1"
            c=${ColorScale}
        else
			legend="${name}"
			Mark="Q"
            c=500
		fi

		# Plot grid.
		gmt makecpt -Cpolar -D -I -T-${c}/${c}/0.5 -Z > tmp.cpt
		find ${WORKDIR}/Calculation/${Model} -iname "*nc" > tmpfile_FileNames_$$
		${EXECDIR}/AddGrid.out 0 7 2 << EOF
${PROJ}
${REG_G}
${Mark}
tmpfile_FileNames_$$
tmp.cpt
${OUTFILE}
tmpfile_$$
0
180
EOF
		gmt psbasemap ${REG_G} ${PROJ} -Ba10f2/a1000f200NWe -O -K >> ${OUTFILE}

		## source.
		gmt psxy ${REG_G} ${PROJ} -Sa0.2i -Gyellow -N -O -K >> ${OUTFILE} << EOF
0 `echo "6371-${EVDE}" | bc -l`
EOF

		## stations.
		awk 'NR>1 {print $1,"6371"}' ${recfile} > tmpfile_$$
		gmt psxy tmpfile_$$ ${REG_G} ${PROJ} -Si0.05i -Gblue -N -O -K >> ${OUTFILE}

		## plot scale bar.
		B=`echo ${c} | awk '{if ($1>2) print $1/5; else print "0.2"}'`
		gmt psscale -Ctmp.cpt -D3.45i/0.18i/3.0i/0.13ih -B${B}/:"@~\144@~${legend} (%)": -Y-0.7i -N300 -O -K >> ${OUTFILE}

	done

	gmt psxy -J -R -O >> ${OUTFILE} << EOF
EOF

	ps2pdf ${OUTFILE} ${PLOTDIR}/${Model}.pdf

done # Done model loop.

# Combine PDFs.
cat `ls *_2.ps | sort -n` > Zoom.ps
# cat `ls *_1.ps | sort -n` > Global.ps
ps2pdf Zoom.ps ${PLOTDIR}/Zoom.pdf
convert -flatten Zoom.ps ${PLOTDIR}/Zoom.png
# ps2pdf Global.ps ${PLOTDIR}/Global.pdf

cd ${WORKDIR}

exit 0
