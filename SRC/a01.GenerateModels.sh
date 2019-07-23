#!/bin/bash

# =======================================
# This script generate ULVZ models !
#
# Shule Yu
# Aug 27 2015
# =======================================

echo ""
echo "--> `basename $0` is running."
mkdir -p ${WORKDIR}/Calculation
cd ${WORKDIR}
trap "rm -rf ${WORKDIR}/tmpfile*${RunNumber}; exit 1" SIGINT

CMB_DEG2KM="60.737458"

# Extract from tar file.
rm -rf ${WORKDIR}/shaxi
tar xf ${SRCDIR}/${OriginFile} --directory=${WORKDIR}

# Make ScS center postion table.

# if ! [ ${ModelType} = "PREM" ]
# then
# 	echo "taup.distance.precision=3" > .taup
# 	rm -f ${WORKDIR}/Calculation/tmpfile_${EVDE}_table_$$
# 	for receiver in `seq ${DISTMIN} 0.2 ${DISTMAX}`
# 	do
# 		taup_path -mod prem -h ${EVDE} -ph ScS -deg ${receiver}
# 		echo "${receiver} `grep 3480 taup_path.gmt | awk 'NR==1 {print $1}'`" >> ${WORKDIR}/Calculation/tmpfile_${EVDE}_table_$$
# 	done
# 
# 	${EXECDIR}/MoreDist.out 1 2 0 << EOF
# `wc -l < ${WORKDIR}/Calculation/tmpfile_${EVDE}_table_$$`
# ${WORKDIR}/Calculation/tmpfile_${EVDE}_table_$$
# ${WORKDIR}/Calculation/tmpfile_$$
# EOF
# 
# 	mv ${WORKDIR}/Calculation/tmpfile_$$ ${WORKDIR}/Calculation/tmpfile_${EVDE}_table_$$
# 	rm -f taup_path.gmt
# fi

# Modify fortran code (fd2_model.f90, mod_global.f90):
cd ${WORKDIR}/shaxi/Source

# 1. Insert "use ulvz"
Line=`grep -n "USE circle" fd2_model.f90 | awk 'BEGIN {FS=":"} {print $1}'`
ed -s fd2_model.f90 >/dev/null << EOF
${Line}
a
  USE ulvz
.
wq
EOF

## Add default no atten to the model.

cat >> mod_attenuation.f90 << EOF
SUBROUTINE Qprem2(Q)
USE global, ONLY: nx,nz
IMPLICIT NONE
REAL, DIMENSION(nx,nz), INTENT(OUT) :: Q
INTEGER :: I

    DO i=1,nz
      Q(:,i) = 1E37
    ENDDO

END SUBROUTINE Qprem2
EOF

##   Insert ulvz model to model_type=300
Line=`grep -n "CALL tomodriver" fd2_model.f90 | awk 'BEGIN {FS=":"} {print $1}'`
ed -s fd2_model.f90 >/dev/null << EOF
${Line}
a

  ELSEIF (model_type == 300) THEN  ! ULVZ
    write(420,*)' Initializing ULVZ model'
    IF (attenuate == 1) THEN
        IF ( ${GlobalAttenuate} == 1) THEN
            CALL Qprem(Q1,1)
            CALL Qprem(Q2,2)
        ELSE
            CALL Qprem2(Q1)
            CALL Qprem2(Q2)
        ENDIF
    ENDIF
    CALL ulvzdriver(rho,mu1,mu2,${StructureAttenuate},${StructureQ},Q1,Q2)
.
wq
EOF

## Crust,220,400,670 km discontinuity removal.
if [ ${RemoveCrust} -eq 1 ] || [ ${Remove220} -eq 1 ] || [ ${Remove400} -eq 1 ] || [ ${Remove670} -eq 1 ] || [ ${SplineTo400} -eq 1 ] || [ ${Straight600To771} -eq 1 ] || [ ${Degree2From600} -eq 1 ]
then
	Line=`grep -n "INTEGER :: k,i" fd2_model.f90 | awk 'BEGIN {FS=":"} {print $1}'`
	ed -s fd2_model.f90 >/dev/null << EOF
${Line}
a
  REAL :: R11, VS11, RHO11, R12, VS12, RHO12
.
wq
EOF
	Line=`grep -n "End modification of crust removal" fd2_model.f90 | awk 'BEGIN {FS=":"} {print $1}'`
	ed -s fd2_model.f90 >/dev/null << EOF
${Line}
a

  ! Marker For More.
.
wq
EOF
fi

if [ ${RemoveCrust} -eq 1 ]
then
	Line=`grep -n "Marker For More." fd2_model.f90 | awk 'BEGIN {FS=":"} {print $1}'`
	ed -s fd2_model.f90 >/dev/null << EOF
${Line}c
      IF (dist >= 6346.6) THEN
          muu = 4.49100712
          ro = 3.38074821
      ENDIF

  ! Marker For More.
.
wq
EOF
fi

if [ ${Remove220} -eq 1 ]
then
	Line=`grep -n "Marker For More." fd2_model.f90 | awk 'BEGIN {FS=":"} {print $1}'`
	ed -s fd2_model.f90 >/dev/null << EOF
${Line}c
      IF (dist.lt.6181.0 .AND. dist.ge.6121.0) THEN
           R11 = 6181.0
           VS11 = 4.42997347
           RHO11 = 3.36275081
           R12 = 6121.0
           VS12 = 4.66490000
           RHO12 = 3.45368975
           m = (VS12-VS11)/(R12-R11)
           muu = m*(dist-R11)+VS11
           m = (RHO12-RHO11)/(R12-R11)
           ro = m*(dist-R11)+RHO11
      ENDIF

  ! Marker For More.
.
wq
EOF
fi

if [ ${Remove400} -eq 1 ]
then
	Line=`grep -n "Marker For More." fd2_model.f90 | awk 'BEGIN {FS=":"} {print $1}'`
	ed -s fd2_model.f90 >/dev/null << EOF
${Line}c
      IF (dist.lt.6001.0 .AND. dist.ge.5941.0) THEN
           R11 = 6001.0
           VS11 = 4.74890000
           RHO11 = 3.52534883
           R12 = 5941.0
           VS12 = 5.02000402
           RHO12 = 3.76155793
           m = (VS12-VS11)/(R12-R11)
           muu = m*(dist-R11)+VS11
           m = (RHO12-RHO11)/(R12-R11)
           ro = m*(dist-R11)+RHO11
      ENDIF

  ! Marker For More.
.
wq
EOF
fi

if [ ${Remove670} -eq 1 ]
then
	Line=`grep -n "Marker For More." fd2_model.f90 | awk 'BEGIN {FS=":"} {print $1}'`
	ed -s fd2_model.f90 >/dev/null << EOF
${Line}c
      IF (dist.lt.5731.0 .AND. dist.ge.5671.0) THEN
           R11 = 5731.0
           VS11 = 5.54698517
           RHO11 = 3.98513532
           R12 = 5671.0
           VS12 = 6.03284432
           RHO12 = 4.39943638
           m = (VS12-VS11)/(R12-R11)
           muu = m*(dist-R11)+VS11
           m = (RHO12-RHO11)/(R12-R11)
           ro = m*(dist-R11)+RHO11
      ENDIF

  ! Marker For More.
.
wq
EOF
fi

if [ ${SplineTo400} -eq 1 ]
then
	Line=`grep -n "Marker For More." fd2_model.f90 | awk 'BEGIN {FS=":"} {print $1}'`
	ed -s fd2_model.f90 >/dev/null << EOF
${Line}c
      IF (dist.ge.5971.0) THEN
           muu = -1147.226569441673*x**3+3481.648723871561*x*x &
                 -3521.617739418105*x+1191.686592108216
           ro = 734.780041069559*x**3-2071.193615789281*x*x &
                +1938.047108369886*x-598.252785440164
      ENDIF

  ! Marker For More.
.
wq
EOF

fi

if [ ${Straight600To771} -eq 1 ]
then
	Line=`grep -n "Marker For More." fd2_model.f90 | awk 'BEGIN {FS=":"} {print $1}'`
	ed -s fd2_model.f90 >/dev/null << EOF
${Line}c

      IF (dist.lt.5771.0 .AND. dist.ge.5600) THEN
           R11=5771
           R12=5600
           VS11 = 5.51593119
           VS12 = 6.24053584
           RHO11 = 3.97582037
           RHO12 = 4.44320419
           m = (VS12-VS11)/(R12-R11)
           muu = m*(dist-R11)+VS11
           m = (RHO12-RHO11)/(R12-R11)
           ro = m*(dist-R11)+RHO11
      ENDIF

  ! Marker For More.
.
wq
EOF

fi

if [ ${Degree2From600} -eq 1 ]
then
	Line=`grep -n "Marker For More." fd2_model.f90 | awk 'BEGIN {FS=":"} {print $1}'`
	ed -s fd2_model.f90 >/dev/null << EOF
${Line}c

      IF (dist.lt.5771.0 .AND. dist.ge.5215.326) THEN
           muu = -84.46554942*x*x+134.4361189*x-46.95411628
      ENDIF

      IF (dist.lt.5771.0 .AND. dist.ge.4644.677) THEN
           ro =  -13.29902926*x*x+16.06334012*x+0.3373366175
      ENDIF

  ! Marker For More.
.
wq
EOF

fi


if [ ${RemoveDepthPhase} -eq 1 ]
then

	Line=`grep -n "initialize 2D damping function window" fd2_init.f90 | awk 'BEGIN {FS=":"} {print $1}'`
	ed -s fd2_init.f90 >/dev/null << EOF
${Line}
a
  ! modified by Shule.
  ! Remove depth phase sS by adding an absorbing boundary
  ! between 0 ~ 30 deg at the surface.

  if((mpi_xrank+1)*(nx-2*BSIZE)*dtheta<30) then
    DO i=1,nabs
      window(:,i)=taper(i)
    ENDDO
  endif

  if (mpi_xrank*(nx-2*BSIZE)*dtheta<30 &
      .AND. (mpi_xrank+1)*(nx-2*BSIZE)*dtheta>30) then
    DO k=1,nx
      if ((mpi_xrank*(nx-2*BSIZE)-BSIZE+k)*dtheta<30) then
        DO i=1,nabs
          window(k,i)=taper(i)
        ENDDO
      endif
    ENDDO
  endif

.
wq
EOF
fi

# 2. Upgrade the maximum of stations.
Line=`grep -n "maxnr=" mod_global.f90 | awk 'BEGIN {FS=":"} {print $1}'`
ed -s mod_global.f90 >/dev/null << EOF
${Line}c
  INTEGER, PARAMETER :: maxnr=500      ! modified by Shule.
.
wq
EOF

##   Change grid size.
Line=`grep -n "^\ \ integer" ./mod_global.f90 | awk 'BEGIN {FS=":"} {print $1}'`
ed -s mod_global.f90 >/dev/null << EOF
${Line}c
  integer, parameter :: XSIZE=${XSIZE}/${RankNum},  ZSIZE=${ZSIZE}   ! modified by Shule.
.
wq
EOF

# Make tmp receivers.
echo `seq ${DISTMIN} ${DISTINC} ${DISTMAX}| wc -l` > ${WORKDIR}/tmpfile.recs_${RunNumber}
for receiver in `seq ${DISTMIN} ${DISTINC} ${DISTMAX}`
do
    printf "%11.6lf\n" ${receiver} >> ${WORKDIR}/tmpfile.recs_${RunNumber}
done

# Make tmp Par file.
Attenuate=0
[ ${GlobalAttenuate} -ne 0 ] && Attenuate=1
[ ${StructureAttenuate} -ne 0 ] && Attenuate=1
cp ${SRCDIR}/Par.txt ${WORKDIR}/tmpfile.Par_${RunNumber}
ed -s ${WORKDIR}/tmpfile.Par_${RunNumber} >/dev/null << EOF
5c
time          =${Length}
.
9c
attenuate     =${Attenuate}
.
12c
model_type    =300
.
16c
Tmax          =${ThetaMax}
.
19,20c
sdepth        =${EVDE}
aa            =${aa}
.
wq
EOF


if [ ${NETCDFLIB} = "stub" ]
then
	NETCDFLIB=""
	FileList="mod_global.f90 mod_cart2D.f90 mod_sacoutput.f90 mod_ncoutput_stub.f90 mod_attenuation.f90 mod_ulvz.f90 mod_simpleslab.f90 mod_simplerift.f90 mod_circle.f90 mod_modifycrust.f90 mod_tomogrand.f90 fd2_output.f90 fd2_model.f90 fd2_main.f90 fd2_input.f90 fd2_init.f90 fd2_allocate.f90 fd2_evolution.f90 fd2_check.f90"
else
	NETCDFLIB="-l${NETCDFLIB}"
	FileList="mod_global.f90 mod_cart2D.f90 mod_sacoutput.f90 mod_ncoutput.f90 mod_attenuation.f90 mod_ulvz.f90 mod_simpleslab.f90 mod_simplerift.f90 mod_circle.f90 mod_modifycrust.f90 mod_tomogrand.f90 fd2_output.f90 fd2_model.f90 fd2_main.f90 fd2_input.f90 fd2_init.f90 fd2_allocate.f90 fd2_evolution.f90 fd2_check.f90"
fi

# ================================================
#      ! Make Models and re-compile SHaxi !
# ================================================

if ! [ ${ModelType} = "PREM" ]
then

	${EXECDIR}/MixModelInput.out 0 2 0 << EOF
${WORKDIR}/tmpfile_DefineShapes_${RunNumber}
${WORKDIR}/Calculation/tmpfile_models$$_
EOF

	if [ $? -ne 0 ]
	then
		echo "    !=> Mix C++ code failed ..."
		exit 1;
	fi
else
    echo "PREM" > ${WORKDIR}/Calculation/tmpfile_models$$_1
fi

ModelCnt=`ls ${WORKDIR}/Calculation/tmpfile_models$$_* | wc -l`

# Model loop.
for count in `seq ${BeginIndex} $((BeginIndex+ModelCnt-1))`
do
    ModelFile=${WORKDIR}/Calculation/tmpfile_models$$_$((count-BeginIndex+1))

    mkdir -p ${WORKDIR}/Calculation/${ModelType}_${count}/OUTPUT
    echo $$ > ${WORKDIR}/Calculation/${ModelType}_${count}/RunNumber.txt
	rm -rf ${WORKDIR}/Calculation/${ModelType}_${count}/OUTPUT/*
    cd ${WORKDIR}/Calculation/${ModelType}_${count}

    # Make receivers.
    cp ${WORKDIR}/tmpfile.recs_${RunNumber} ${WORKDIR}/Calculation/${ModelType}_${count}/recs

    # Make ${ModelType}_${count}.par
    cp ${WORKDIR}/tmpfile.Par_${RunNumber} Par_${ModelType}
    
    # OUTDIR has a length limit: 50
	[ ${RunOnOffice} -eq 0 ] && OUTDIR="/local/shule/${ModelType}_$$_${count}/OUTPUT/${ModelType}_${count}" || OUTDIR="OUTPUT/${ModelType}_${count}"

    ed -s Par_${ModelType} >/dev/null << EOF
2c
seisfile      =${OUTDIR}
.
32c
kevnm         =${ModelType}_${count}
.
wq
EOF

    if [ ${RunOnOffice} -eq 2 ]
    then

cat > ${WORKDIR##*/}_${ModelType}_${count}.sh << EOF
#!/bin/bash
 
#SBATCH -N 1                        # number of nodes 
#SBATCH -n ${RankNum}               # number of cores
#SBATCH -t ${WallTime}              # wall time (D-HH:MM)
##SBATCH -A shuleyu                 # Account hours will be pulled from (commented out with double # in front)
#SBATCH -o slurm.%j.out             # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.err             # STDERR (%j = JobId)
#SBATCH --mail-type=FAIL            # Send a notification when the job fails
#SBATCH --mail-user=shuleyu@asu.edu # send-to address

module purge

srun --mpi=pmi2 ./shaxi_${count}.out < ./Par_ULVZ
EOF

    fi

    # Make mod_ulvz_${count}.f90
    cp ${SRCDIR}/ULVZ.f90 mod_${ModelType}_${count}.f90
    Line=`grep -n "MarkerHere" mod_${ModelType}_${count}.f90 | awk 'BEGIN {FS=":"} {print $1}'`

    if [ "${ModelType}" = ULVZ ]
    then

        rm -f tmpfile_$$
        while read line
        do
            Shape=`echo ${line} | awk '{print $1}'`

            Rulvz=`    echo ${line} | awk '{print 6371-$2}'`
            Thickness=`echo ${line} | awk '{print $3}'`
            # NOTE: ulvz properties in %, for example dvs=-20.0
            Vs=`       echo ${line} | awk '{printf "%.2lf",($4-1)*100}'`
            Rho=`      echo ${line} | awk '{printf "%.2lf",($5-1)*100}'`
            EPulvz=`   echo ${line} | awk '{print $6}'`

            case "${Shape}" in

                Flat ) # 1. Box-car

                    LSize=`    echo ${line} | awk '{print $7}'`

                    cat >> tmpfile_$$ << EOF
Rulvz      = ${Rulvz}
Thickness  = ${Thickness}
EPulvz     = ${EPulvz}
LSize      = ${LSize}

IF (r <= Rulvz+Thickness &
.AND. r >= Rulvz &
.AND. theta >= EPulvz-LSize/2/CMB_deg2km &
.AND. theta <= EPulvz+LSize/2/CMB_deg2km) THEN
flag=1
dvs=${Vs}
drho=${Rho}
ENDIF
EOF
                    ;;

                Tilting )  # 2. Tilting top Layer ULVZ.

                    LSize=`    echo ${line} | awk '{print $7}'`
                    TiltAngle=`echo ${line} | awk '{print $8}'`

                    cat >> tmpfile_$$ << EOF
Rulvz      = ${Rulvz}
Thickness  = ${Thickness}
EPulvz     = ${EPulvz}
LSize      = ${LSize}
TiltAngle  = ${TiltAngle}*deg2rad

IF (r <= Rulvz+Thickness+(theta-EPulvz)*CMB_deg2km*TAN(TiltAngle) &
.AND. r >= Rulvz &
.AND. theta >= EPulvz-LSize/2/CMB_deg2km &
.AND. theta <= EPulvz+LSize/2/CMB_deg2km ) THEN
flag=1
dvs=${Vs}
drho=${Rho}
ENDIF
EOF
                    ;;

                Trapezoid )  # 3. Trapezoid shape ULVZ.

                    UpperSize=`echo ${line} | awk '{print $7}'`
                    LowerSize=`echo ${line} | awk '{print $8}'`

                    cat >> tmpfile_$$ << EOF
Rulvz      = ${Rulvz}
Thickness  = ${Thickness}
EPulvz     = ${EPulvz}
ESize      = (${LowerSize}-${UpperSize})/2/CMB_deg2km
LSize      = ${LowerSize}

IF (r <= Rulvz+Thickness .AND. r >= Rulvz &
.AND. ABS(theta-EPulvz) <= LSize/2/CMB_deg2km-ESize) THEN
flag=1
dvs=${Vs}
drho=${Rho}
ELSE IF (theta >= EPulvz-LSize/2/CMB_deg2km &
.AND.   theta <= EPulvz+LSize/2/CMB_deg2km &
.AND. r >= Rulvz &
.AND. r <= Rulvz+Thickness-ABS(ABS(theta-EPulvz)+ESize-LSize/2/CMB_deg2km)*Thickness/ESize) THEN
flag=1
dvs=${Vs}
drho=${Rho}
ENDIF
EOF

                    ;;

                Hill ) # 4. Asymmetric triangle shape ULVZ.

                    SourceSize=`  echo ${line} | awk '{print $7}'`
                    ReceiverSize=`echo ${line} | awk '{print $8}'`

                    cat >> tmpfile_$$ << EOF
Rulvz      = ${Rulvz}
Thickness  = ${Thickness}
EPulvz     = ${EPulvz}
LSize      = ${SourceSize}
RSize      = ${ReceiverSize}

IF ( theta >= EPulvz-LSize/CMB_deg2km .AND. theta < EPulvz .AND. r >= Rulvz &
.AND. r <= Rulvz + ( theta-EPulvz+LSize/CMB_deg2km )*CMB_deg2km/LSize*Thickness ) THEN
flag=1
dvs=${Vs}
drho=${Rho}
ELSE IF ( theta >= EPulvz .AND. theta < EPulvz+RSize/CMB_deg2km .AND. r >= Rulvz &
.AND. r <= Rulvz + ( EPulvz+RSize/CMB_deg2km-theta )*CMB_deg2km/RSize*Thickness ) THEN
flag=1
dvs=${Vs}
drho=${Rho}
ENDIF
EOF
                    ;;

                Doom ) # 5. Ellipse shape ULVZ.

                    LSize=`    echo ${line} | awk '{print $7}'`

                    cat >> tmpfile_$$ << EOF
Rulvz      = ${Rulvz}
Thickness  = ${Thickness}
EPulvz     = ${EPulvz}
LSize      = ${LSize}

IF (r <= Rulvz+Thickness .AND. r >= Rulvz &
.AND. theta >= EPulvz-LSize/2/CMB_deg2km &
.AND. theta <= EPulvz+LSize/2/CMB_deg2km &
.AND. r <= Rulvz+Thickness*SQRT(1-(2*(theta-EPulvz)*CMB_deg2km/LSize)**2) ) THEN
flag=1
dvs=${Vs}
drho=${Rho}
ENDIF
EOF
                    ;;

                Mollifier ) # 6. Mollifier at the boundary of flat layer ULVZ.

                    UpperSize=`echo ${line} | awk '{print $7}'`
                    LowerSize=`echo ${line} | awk '{print $8}'`

                    cat >> tmpfile_$$ << EOF
Rulvz      = ${Rulvz}
Thickness  = ${Thickness}
EPulvz     = ${EPulvz}
LSize      = ${UpperSize}
Delta      = (${LowerSize}-${UpperSize})/2/CMB_deg2km

IF (r <= Rulvz+Thickness .AND. r >= Rulvz &
.AND. theta >= EPulvz-LSize/2/CMB_deg2km &
.AND. theta <= EPulvz+LSize/2/CMB_deg2km) THEN
flag=1
dvs=${Vs}
drho=${Rho}
ENDIF

tmp        = EPulvz-LSize/2/CMB_deg2km

IF (EPulvz-LSize/2/CMB_deg2km-Delta <= theta &
.AND. theta <= EPulvz-LSize/2/CMB_deg2km &
.AND. r >= Rulvz &
.AND. r<=Rulvz+Thickness*EXP(1.0)*EXP(-1.0/(1-((theta-tmp)/Delta)**2))) THEN
flag=1
dvs=${Vs}
drho=${Rho}
ENDIF

tmp        = EPulvz+LSize/2/CMB_deg2km

IF (EPulvz+LSize/2/CMB_deg2km <= theta &
.AND. theta <= EPulvz+LSize/2/CMB_deg2km+Delta &
.AND. r >= Rulvz &
.AND. r<=Rulvz+Thickness*EXP(1.0)*EXP(-1.0/(1-((theta-tmp)/Delta)**2))) THEN
flag=1
dvs=${Vs}
drho=${Rho}
ENDIF
EOF
                    ;;

                Gaussian ) # 7. Gaussian shape.

                    UpperSize=`          echo ${line} | awk '{print $7}'`
                    LowerSize=`          echo ${line} | awk '{print $8}'`
                    HalfHeightHalfWidth=`echo ${line} | awk '{print $9}'`

                    cat >> tmpfile_$$ << EOF
Rulvz      = ${Rulvz}
Thickness  = ${Thickness}
EPulvz     = ${EPulvz}
LSize      = ${UpperSize}
RSize      = ${LowerSize}
Delta      = ${HalfHeightHalfWidth}/CMB_deg2km
Delta      = sqrt(Delta*Delta/2.0/log(2.0))

IF (  theta >= EPulvz-LSize/2/CMB_deg2km &
.AND. theta <= EPulvz+LSize/2/CMB_deg2km &
.AND. r >= Rulvz .AND. r <= Rulvz+Thickness ) THEN
flag=1
dvs=${Vs}
drho=${Rho}
ENDIF

IF (  theta >= EPulvz-RSize/2/CMB_deg2km &
.AND. theta <= EPulvz-LSize/2/CMB_deg2km &
.AND. r >= Rulvz & 
.AND. r <= Rulvz+Thickness*exp(-0.5*((theta-EPulvz+LSize/2/CMB_deg2km)/Delta)**2)) THEN
flag=1
dvs=${Vs}
drho=${Rho}
ENDIF

IF (  theta >= EPulvz+LSize/2/CMB_deg2km &
.AND. theta <= EPulvz+RSize/2/CMB_deg2km &
.AND. r >= Rulvz & 
.AND. r <= Rulvz+Thickness*exp(-0.5*((theta-EPulvz-LSize/2/CMB_deg2km)/Delta)**2)) THEN
flag=1
dvs=${Vs}
drho=${Rho}
ENDIF
EOF

                    ;;


                *)
                    echo "    !=> Structure shape code wrong, exiting ..."
                    rm -f ${WORKDIR}/Calculation/tmpfile*$$ ${WORKDIR}/Calculation/tmpfile_models$$_*
                    exit 1
                    ;;
            esac


        done < ${ModelFile}

        awk -v L=${Line} 'NR<L {print $0}' mod_${ModelType}_${count}.f90 > x
        cat tmpfile_$$ >> x
        awk -v L=${Line} 'NR>L {print $0}' mod_${ModelType}_${count}.f90 >> x
        mv x mod_${ModelType}_${count}.f90
        rm -f tmpfile_$$
    fi


    cp mod_${ModelType}_${count}.f90 ${WORKDIR}/shaxi/Source/mod_ulvz.f90

    # Re-compile SHaxi for this Model.
    echo "    ==> Compiling ${ModelType}_${count} ..."
    cd ${WORKDIR}/shaxi/Source
	if [ ${RunOnOffice} -eq 1 ]
	then
		${FCOMP} ${FFLAGS} -o ${WORKDIR}/Calculation/${ModelType}_${count}/shaxi_${count}.out ${FileList} ${NETCDFLIB}
	else
		${FCOMP} ${FFLAGS} -o ${WORKDIR}/Calculation/${ModelType}_${count}/shaxi_${count}.out ${FileList} ${NETCDFLIB} -I${NETCDFINCLUDEDIR} -L${NETCDFLIBDIR}
	fi

    if [ $? -ne 0 ]
    then
        echo "    !=> SHaxi didn't compile successfully on ${ModelType}_${count} !"
        exit 1;
    fi


done # Done Model loop.

# Make index files.
mkdir -p ${WORKDIR}/indices
for count in `seq 1 ${ModelCnt}`
do
    EQname=`echo "201500000000 + ${count} + ${BeginIndex} - 1" | bc `
    mv ${WORKDIR}/Calculation/tmpfile_models$$_${count} ${WORKDIR}/indices/${EQname}
done

rm -f ${WORKDIR}/Calculation/tmpfile*$$

cd ${WORKDIR}

exit 0
