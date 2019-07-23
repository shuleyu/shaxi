#!/bin/bash

rm -f *sac

count=1
while read file
do
	sac << EOF
r ${file}
lh delta gcarc
evaluate to begin &1,t0-510
evaluate to end &1,t0-60
hp co 0.02 n 2 p 2
w junk.sac
cut %begin %end
r junk.sac
rtr
decimate 4
w noise${count}.sac
q
EOF
	count=$((count+1))
done << EOF
/home/shule/soduser/200601022213/200601022213.TA.Y12C.BHT.sac
/home/shule/soduser/200601022213/200601022213.TA.Y22C.BHT.sac
/home/shule/soduser/200601022213/200601022213.TA.M02C.BHT.sac
/home/shule/soduser/200601022213/200601022213.TA.O03C.BHT.sac
/home/shule/soduser/200601022213/200601022213.TA.O04C.BHT.sac
/home/shule/soduser/200601022213/200601022213.TA.O05C.BHT.sac
/home/shule/soduser/200601022213/200601022213.TA.P01C.BHT.sac
/home/shule/soduser/200601022213/200601022213.TA.P05C.BHT.sac
/home/shule/soduser/200601022213/200601022213.TA.Q03C.BHT.sac
/home/shule/soduser/200601022213/200601022213.TA.R04C.BHT.sac
/home/shule/soduser/200601022213/200601022213.TA.R06C.BHT.sac
/home/shule/soduser/200601022213/200601022213.TA.S04C.BHT.sac
/home/shule/soduser/200601022213/200601022213.TA.S06C.BHT.sac
/home/shule/soduser/200601022213/200601022213.TA.T05C.BHT.sac
/home/shule/soduser/200601022213/200601022213.TA.U04C.BHT.sac
/home/shule/soduser/200601022213/200601022213.TA.U05C.BHT.sac
/home/shule/soduser/200601022213/200601022213.TA.V03C.BHT.sac
/home/shule/soduser/200707211327/200707211327.BK.BDM.BHT.sac   
/home/shule/soduser/200707211327/200707211327.BK.BKS.BHT.sac   
/home/shule/soduser/200707211327/200707211327.BK.BRIB.BHT.sac  
/home/shule/soduser/200707211327/200707211327.BK.CMB.BHT.sac   
/home/shule/soduser/200707211327/200707211327.BK.CVS.BHT.sac   
/home/shule/soduser/200707211327/200707211327.BK.ELFS.BHT.sac  
/home/shule/soduser/200707211327/200707211327.BK.GASB.BHT.sac  
/home/shule/soduser/200707211327/200707211327.BK.HAST.BHT.sac  
/home/shule/soduser/200707211327/200707211327.BK.HATC.BHT.sac  
/home/shule/soduser/200707211327/200707211327.BK.HELL.BHT.sac  
/home/shule/soduser/200707211327/200707211327.BK.HUMO.BHT.sac  
/home/shule/soduser/200707211327/200707211327.BK.JCC.BHT.sac   
/home/shule/soduser/200707211327/200707211327.BK.KCC.BHT.sac   
EOF


rm -r junk.sac

exit 0

