#!/bin/bash

rm -f *.ps *.psx
c++ Post.cpp -o Post.out -lm
MinSnap=197
MaxSnap=307
RankXN=16
RankYN=1
DT="2"

gmt gmtset PS_MEDIA letter

cat > tmp.cpt << EOF
-1	0/0/255	0	255/255/255
0	255/255/255	1	255/0/0
B	0/0/255
F	255/0/0
N	0/0/0
EOF

# ULVZ structure.
DegMin=`echo 29.189 300 60.737 | awk '{print $1-$2/$3}'`
DegMax=`echo 29.189 300 60.737 | awk '{print $1+$2/$3}'`
DegMin1=`echo 29.189 200 60.737 | awk '{print $1-$2/$3}'`
DegMax1=`echo 29.189 200 60.737 | awk '{print $1+$2/$3}'`
Delta=1.64644286

echo "${DegMin} 3480" | awk '{print -$1,$2}' > ulvz.tmp
for Deg in `seq ${DegMin} 0.1 ${DegMin1}`
do
	echo "${Deg} ${DegMin1} ${Delta}" | awk '{print -$1,3480+20*exp(1-1/(1-(($1-$2)/$3)**2))}' >> ulvz.tmp
done
echo "${DegMin1}" | awk '{print -$1,3500}' >> ulvz.tmp
for Deg in `seq ${DegMin1} 0.1 ${DegMax1}`
do
	echo "${Deg} 3500" | awk '{print -$1,$2}' >> ulvz.tmp
done
echo "${DegMax1}" | awk '{print -$1,3500}' >> ulvz.tmp
for Deg in `seq ${DegMax1} 0.1 ${DegMax}`
do
	echo "${Deg} ${DegMax1} ${Delta}" | awk '{print -$1,3480+20*exp(1-1/(1-(($1-$2)/$3)**2))}' >> ulvz.tmp
done
echo "${DegMax}" | awk '{print -$1,3480}' >> ulvz.tmp
for Deg in `seq ${DegMax} -0.1 ${DegMin}`
do
	echo "${Deg} 3480" | awk '{print -$1,$2}' >> ulvz.tmp
done

# Boundary
echo "0 6371" > boundary.tmp
for Deg in `seq 0 0.1 110`
do
	echo "-${Deg} 6371" >> boundary.tmp
done
echo "-110 3480" >> boundary.tmp
for Deg in `seq 110 -0.1 0`
do
	echo "-${Deg} 3480" >> boundary.tmp
done
echo "0 6371" >> boundary.tmp

# Stations.
rm -f stations.tmp
for Deg in `seq 40 5 85`
do
	echo "-${Deg} 6371" >> stations.tmp
done

# Ray paths.
echo "taup.distance.precision=3" > .taup
taup_path -mod prem -h 500 -ph ScS -deg 50 -o stdout | awk 'NR>1 {print -$1,$2}' > path.tmp
echo ">" >> path.tmp
taup_path -mod prem -h 500 -ph ScS -deg 60 -o stdout | awk 'NR>1 {print -$1,$2}' >> path.tmp
echo ">" >> path.tmp
taup_path -mod prem -h 500 -ph ScS -deg 70 -o stdout | awk 'NR>1 {print -$1,$2}' >> path.tmp


# Plot !
for ns in `seq ${MinSnap} ${MaxSnap}`
do
	echo "Plotting Snapshot ${ns} ..."
	ns=`echo ${ns} | awk '{printf "%03d",$1}'`

	OUTFILE=${ns}.ps
	PROJ="-JP9i/215"
	REG="-R-110/0/3480/6371"


	OUTFILE2=${ns}.psx
	PROJ2="-JX-9i/2.5i"
	REG2="-R-39.18/-19.18/3480/3821.58"


	# Plot basemap & text.

	gmt psbasemap ${PROJ} ${REG} -Bnews -X1i -Y2i -K > ${OUTFILE}
	gmt psbasemap ${PROJ2} ${REG2} -Ba2f1/a20f10WSne -X1i -Y2i -K > ${OUTFILE2}

	cat > text.tmp << EOF
-55 1700 `echo ${DT} ${ns} | awk '{print $1*$2}'` (sec)
EOF
	gmt pstext text.tmp ${PROJ} ${REG} -F+jCB+f15p,Times-Bold -N -O -K >> ${OUTFILE}

	cat > text.tmp << EOF
-29.18 3900 `echo ${DT} ${ns} | awk '{print $1*$2}'` (sec)
EOF
	gmt pstext text.tmp ${PROJ2} ${REG2} -F+jCB+f15p,Times-Bold -N -O -K >> ${OUTFILE2}

	# generate data.
	rm -f List.tmp
	for nx in `seq 0 $((RankXN-1))`
	do
		nx=`echo ${nx} | awk '{printf "%02d",$1}'`
		for ny in `seq 0 $((RankYN-1))`
		do
			ny=`echo ${ny} | awk '{printf "%02d",$1}'`

			DataFile=`ls *rx${nx}rz${ny}snap${ns}`
			echo "${DataFile} ${nx} ${ny}" >> List.tmp
		done
	done
	./Post.out List.tmp
	[ $? -ne 0 ] && exit 1

	# Plot snap shot.
	read xmin xmax xinc ymin ymax yinc < tmp.desc
	gmt xyz2grd tmp.dat -Gtmp.grd -R${xmin}/${xmax}/${ymin}/${ymax} -I${xinc}/${yinc}
	gmt grdimage -Q ${PROJ} ${REG} tmp.grd -Ctmp.cpt -O -K >> ${OUTFILE}
	gmt grdimage -Q ${PROJ2} ${REG2} tmp.grd -Ctmp.cpt -O -K >> ${OUTFILE2}

	# Add ray paths.
	gmt psxy path.tmp -W0.1p,gray,- ${PROJ} ${REG} -O -K >> ${OUTFILE}
	gmt psxy path.tmp -W0.1p,gray,- ${PROJ2} ${REG2} -O -K >> ${OUTFILE2}


	# Add ULVZ.
	gmt psxy ulvz.tmp -Ggreen ${PROJ} ${REG} -O -K >> ${OUTFILE}
	gmt psxy ulvz.tmp -Wgreen,- ${PROJ2} ${REG2} -O -K >> ${OUTFILE2}


	# plot boundary.
	gmt psxy boundary.tmp ${PROJ} ${REG} -W1p,black -O -K >> ${OUTFILE}

	# Add stations.
	gmt psxy stations.tmp -Sc0.05i -Ggreen ${PROJ} ${REG} -Wblack -N -O -K >> ${OUTFILE}

	# Seal it.
	gmt psxy ${PROJ} ${REG} -O >> ${OUTFILE} <<EOF
EOF
	gmt psxy ${PROJ2} ${REG2} -O >> ${OUTFILE2} <<EOF
EOF

done # done snapshot loop.

cat `ls -rt *ps` > tmp.ps
ps2pdf tmp.ps global.pdf

cat `ls -rt *psx` > tmp.ps
ps2pdf tmp.ps zoom.pdf

rm -f *.dat *.desc *.grd tmp.ps
exit 0
