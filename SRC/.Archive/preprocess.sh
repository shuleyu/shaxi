#!/bin/bash

a=`grep -n ";" all.ncdump | awk 'BEGIN {FS=":"} {print $1}'`
e=`echo ${a} |awk '{print $1}'`

# velocity data
Vs1=`grep -n "Vs1 =" all.ncdump | awk 'BEGIN {FS=":"} {print $1}'`
Vs2=`grep -n "Vs2 =" all.ncdump | awk 'BEGIN {FS=":"} {print $1}'`

count=${Vs2}
while [ "${count}" -gt "${e}" ]
do
    echo $a | grep "[[:blank:]]${count}[[:blank:]]" > /dev/null
    if [ $? -eq 0 ]
    then
        Vs2=${count}
        break
    fi
    count=$((count-1))
done

sed -n $((Vs1+1)),${Vs2}p all.ncdump > Vs.dat
ed -s Vs.dat << EOF
%s/;/,/g
wq
EOF

# R data.
R1=`grep -n "r =" all.ncdump | awk 'BEGIN {FS=":"} NR==2 {print $1}'`
R2=`grep -n "Vs1 =" all.ncdump | awk 'BEGIN {FS=":"} {print $1}'`

count=${R2}
while [ "${count}" -gt "${e}" ]
do
    echo $a | grep "[[:blank:]]${count}[[:blank:]]" > /dev/null
    if [ $? -eq 0 ]
    then
        R2=${count}
        break
    fi
    count=$((count-1))
done

sed -n ${R1},${R2}p all.ncdump > R.dat
ed -s R.dat << EOF
%s/;/,/g
%s/r =//g
wq
EOF

# Theta data.
theta1=`grep -n "theta =" all.ncdump | awk 'BEGIN {FS=":"} NR==2 {print $1}'`
theta2=`grep -n "r =" all.ncdump | awk 'BEGIN {FS=":"} NR==2 {print $1}'`

count=${theta2}
while [ "${count}" -gt "${e}" ]
do
    echo $a | grep "[[:blank:]]${count}[[:blank:]]" > /dev/null
    if [ $? -eq 0 ]
    then
        theta2=${count}
        break
    fi
    count=$((count-1))
done

sed -n ${theta1},${theta2}p all.ncdump > Theta.dat
ed -s Theta.dat << EOF
%s/;/,/g
%s/theta =//g
wq
EOF

# Density data
Rho1=`grep -n "Rho =" all.ncdump | awk 'BEGIN {FS=":"} {print $1}'`
Rho2=`grep -n "Q =" all.ncdump | awk 'BEGIN {FS=":"} {print $1}'`

count=${Rho2}
while [ "${count}" -gt "${e}" ]
do
    echo $a | grep "[[:blank:]]${count}[[:blank:]]" > /dev/null
    if [ $? -eq 0 ]
    then
        Rho2=${count}
        break
    fi
    count=$((count-1))
done

sed -n $((Rho1+1)),${Rho2}p all.ncdump > Rho.dat
ed -s Rho.dat << EOF
%s/;/,/g
wq
EOF

exit 0
