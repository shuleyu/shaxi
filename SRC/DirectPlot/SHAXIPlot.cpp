#include<iostream>
#include<fstream>
#include<sstream>
#include<cstdio>
#include<cstdlib>
#include<vector>
#include<string>
#include<cmath>
#include<algorithm>
#include<map>
#include<sys/types.h>
#include<unistd.h>

#include<GMTPlotSignal.hpp>
#include<SACSignals.hpp>
#include<ShellExecVec.hpp>
#include<ShellExec.hpp>
#include<PREM.hpp>

extern "C" {
#include<netcdf.h>
}

using namespace std;

string d2s(const double &val, const size_t &N) {
    stringstream ss;
    ss << fixed << setprecision(N) << val;
    return ss.str();
}

int main(){



    // const string premDIR="/home/shule/PROJ/t052.MultipleGaussian.ULVZ/201500000022";
    const string premDIR="/home/shule/PROJ/t052.PREMX.PREM/201500000000";

    // const string chosenDIR="t052.LateralSizeAndPosition.ULVZ_Flat";
    // const string chosenDIR="t052.Shapes";
    // const string chosenDIR="t052.MultipleGaussian.ULVZ";
    const string chosenDIR="t052.QTest.ULVZ";

    const vector<string> Phases{"S","ScS"};
    vector<size_t> Models{12,13,14,15,16};

    const string dataDIR="/home/shule/PROJ/"+chosenDIR;
    const string ncDIR="/home/shule/PROJ/"+chosenDIR+"/Calculation/ULVZ_";
    const string indexDIR="/home/shule/PROJ/"+chosenDIR+"/indices";

    const long EQ=201500000000;
    const double f1=0.033,f2=0.3;
    const double dt=0.05,cutWindowT1=-50,cutWindowT2=50,receiverDistInterval=5;
    const vector<double> frsTickX=CreateGrid(cutWindowT1,cutWindowT2,receiverDistInterval,1);
    const double plotThetaMin=15,plotThetaMax=50,plotHeightMin=0,plotHeightMax=50;

    /****************************************************************

                              Job begin.

    ****************************************************************/


    for (string phase: Phases){

        // PREM waveforms.
        SACSignals premSACs(ShellExecVec("ls "+premDIR+"/*T.sac"));
        premSACs.Interpolate(dt);
        premSACs.FindPeakAround(premSACs.GetTravelTimes("S"));
        premSACs.NormalizeToPeak();
        premSACs.FlipPeakUp();

        premSACs.ShiftTime(premSACs.GetTravelTimes(phase));
        premSACs.CheckAndCutToWindow(cutWindowT1-100,cutWindowT2+100);
        premSACs.HannTaper();
        premSACs.Butterworth(f1,f2,2,2);
        premSACs.CheckAndCutToWindow(cutWindowT1,cutWindowT2);

        premSACs.FindPeakAround(0,5,false);

        auto premDistances=premSACs.GetDistances();
        if (phase=="ScS") premSACs*=3;

        for (size_t Cnt: Models) {

            string eq=to_string(EQ+Cnt);

            // Structure waveforms.
            SACSignals structureSACs(ShellExecVec("ls "+dataDIR+"/"+eq+"/*T.sac"));
            structureSACs.Interpolate(dt);
            structureSACs.FindPeakAround(structureSACs.GetTravelTimes("S"));
            structureSACs.NormalizeToPeak();
            structureSACs.FlipPeakUp();

            structureSACs.ShiftTime(structureSACs.GetTravelTimes(phase));
            structureSACs.CheckAndCutToWindow(cutWindowT1-100,cutWindowT2+100);
            structureSACs.HannTaper();
            structureSACs.Butterworth(f1,f2,2,2);
            structureSACs.CheckAndCutToWindow(cutWindowT1,cutWindowT2);

            structureSACs.FindPeakAround(0,5,false);

            auto structureDistances=structureSACs.GetDistances();
            if (phase=="ScS") structureSACs*=3;

            // Plotting.
            string outfile=to_string(Cnt)+"_"+phase+".ps";
            GMT::set("PS_MEDIA 20ix10i");
            GMT::set("FONT_ANNOT 8p");
            GMT::set("FONT_LABEL 12p");
            GMT::set("MAP_ANNOT_OFFSET_PRIMARY 6p");
            GMT::set("MAP_FRAME_PEN 0.4p,black");
            GMT::set("MAP_TICK_PEN_PRIMARY 0.4p,black");
            GMT::BeginPlot(outfile,"-P",ShellExec("pwd",true)+"/"+string(__FILE__));
            GMT::MoveReferencePoint(outfile,"-Xf1i -Yf5i");
            GMT::makecpt("-Cpolar -T-30/30/0.5 -I -Z > tmp.cpt");

            // Plot grid.
            auto gridFileNames=ShellExecVec("ls "+ncDIR+to_string(Cnt)+"/OUTPUT/*nc");

            int ncid,varid,thetaNum=20000/gridFileNames.size()+8,radiusNum=3808;
            float *t= new float [thetaNum];
            float *r= new float [radiusNum];
            float *val= new float [radiusNum*thetaNum];
            for (auto ncfile:gridFileNames) {

                // Read in NC files.
                nc_open(ncfile.c_str(),NC_NOWRITE,&ncid);

                nc_inq_varid(ncid,"theta",&varid);
                nc_get_var_float(ncid,varid,t);

                // this grid file is not within our plot range.
                if (t[thetaNum-1]<plotThetaMin || t[0]>plotThetaMax){
                    nc_close(ncid);
                    continue;
                }

                nc_inq_varid(ncid,"r",&varid);
                nc_get_var_float(ncid,varid,r);

                nc_inq_varid(ncid,"Vs1",&varid);
                nc_get_var_float(ncid,varid,val);

                nc_close(ncid);

                // Make grid.
                vector<vector<double>> plotGrid;
                for (size_t i=0;i<radiusNum*thetaNum;++i)
                    plotGrid.push_back({t[i/radiusNum],r[i%radiusNum]-3480,(val[i]/RvsX(r[i%radiusNum],false)*100-100)});
                GMT::grdimage(outfile,plotGrid,t[1]-t[0],r[0]-r[1],"-JX18i/4i -R"+to_string(plotThetaMin)+"/"+to_string(plotThetaMax)+"/"+to_string(plotHeightMin)+"/"+to_string(plotHeightMax)+" -Ctmp.cpt -O -K");
            }
            delete t;
            delete r;
            delete val;

            GMT::psbasemap(outfile,"-J -R -Bxa5f1 -Bya10f5 -BWSne -O -K");


            // Add ScS paths, ScS bounce point.
            map<double,double> receiverDistMapToBounceDist;
            vector<double> bounceTickX;
            vector<GMT::Text> plotTexts;
            for (double receiverDist=40;receiverDist<=85;receiverDist+=1) {

                auto res=ShellExecVec("taup_path -mod prem -h 500 -ph ScS -deg "+to_string(receiverDist)+" -o stdout | grep -v \">\"");
                vector<double> rayPathTheta,rayPathHeight;
                double bounceDist=0.0,maxRadius=6371;
                for (const string &line: res) {
                    stringstream ss(line);
                    double tt,rr;
                    while (ss >> tt >> rr) {
                        if (rr>3480+plotHeightMax+50) continue;
                        rayPathTheta.push_back(tt);
                        rayPathHeight.push_back(rr-3480);
                        if (rr<maxRadius){
                            maxRadius=rr;
                            bounceDist=tt;
                        }
                    }
                }
                receiverDistMapToBounceDist[receiverDist]=bounceDist;
                if (((int)receiverDist)%5==0) {
                    plotTexts.push_back(GMT::Text(bounceDist,0,d2s(receiverDist,0),16,"CT"));
                    bounceTickX.push_back(bounceDist);
                }
                GMT::psxy(outfile,rayPathTheta,rayPathHeight,"-J -R -W0.5p,black -O -K");
            }

            // Add receivers.
            GMT::psxy(outfile,bounceTickX,vector<double> (bounceTickX.size(),0),"-J -R -Sc0.1i -Gdarkgreen -W0.5p -N -O -K");
            GMT::MoveReferencePoint(outfile,"-Y-0.2i");
            GMT::pstext(outfile,plotTexts,"-J -R -N -O -K");


            // Make FRS.

            auto premFrsSACs=premSACs;
            premFrsSACs.FlipReverseSum(premSACs.PeakTime());
            premFrsSACs.CheckAndCutToWindow(0,15-dt*0.8);


            auto structureFrsSACs=structureSACs;
            structureFrsSACs.FlipReverseSum(structureSACs.PeakTime());
            structureFrsSACs.CheckAndCutToWindow(0,15-dt*0.8);

            // normalize to ScS amplitude.
            premFrsSACs.AmplitudeDivision(premSACs.PeakAmp());
            structureFrsSACs.AmplitudeDivision(structureSACs.PeakAmp());

            // Plot waveforms and distance.
            double XBegin=2.6,YBegin=3.5,XINC=1.27,YINC=0.7;
            auto premPeakAmps=premSACs.PeakAmp();
            auto premPeakTimes=premSACs.PeakTime();
            auto structurePeakAmps=structureSACs.PeakAmp();
            auto structurePeakTimes=structureSACs.PeakTime();

            for (double receiverDist=40;receiverDist<=85;receiverDist+=1) {

                size_t index1=distance(premDistances.begin(),find(premDistances.begin(),premDistances.end(),receiverDist));
                size_t index2=distance(structureDistances.begin(),find(structureDistances.begin(),structureDistances.end(),receiverDist));

                double XP=XBegin+((int)receiverDist-40)/5*XINC,YP=YBegin-((int)receiverDist-40)%5*YINC;
                GMT::MoveReferencePoint(outfile,"-Xf"+to_string(XP)+"i -Yf"+to_string(YP)+"i");

                GMT::psxy(outfile,vector<double> {-15,15,15,-15,-15}, vector<double> {-0.5,-0.5,1,1,-0.5},"-JX1.1i/0.5i -R"+to_string(cutWindowT1)+"/"+to_string(cutWindowT2)+"/-1/1 -L -Glightgray -O -K");
                GMT::psxy(outfile,vector<double> {cutWindowT1,cutWindowT2},vector<double> {0,0},"-J -R -W0.5p,black,. -O -K");
                GMT::psxy(outfile,frsTickX,vector<double> (frsTickX.size(),0),"-J -R -Sy0.02i -W0.5p,black -O -K");

                // Add the peak position.
                GMT::psxy(outfile,premSACs,index1,"-J -R -W1p,black -O -K");
                GMT::psxy(outfile,vector<double> {premPeakTimes[index1]},vector<double> {premPeakAmps[index1]},"-J -R -Sy0.1i -W0.1p,cyan -O -K");
                GMT::psxy(outfile,structureSACs,index2,"-J -R -W1p,red -O -K");
                GMT::psxy(outfile,vector<double> {structurePeakTimes[index2]},vector<double> {structurePeakAmps[index2]},"-J -R -Sy0.1i -W0.1p,green -O -K");

                plotTexts.clear();
                plotTexts.push_back(GMT::Text(cutWindowT1,1,d2s(receiverDist,0),8,"LT"));
                GMT::pstext(outfile,plotTexts,"-J -R -N -O -K");

                // Plot FRS.
                GMT::psxy(outfile,premFrsSACs,index1,"-JX0.55i/0.5i -R0/"+to_string(cutWindowT2)+"/-1/1 -W1p,black -O -K -X0.55i -Y-0.2i");
                GMT::psxy(outfile,structureFrsSACs,index2,"-J -R -W1p,red -O -K");

                GMT::psxy(outfile,structureFrsSACs.GetData()[index2]-premFrsSACs.GetData()[index1],"-J -R -W1p,purple -O -K -X-0.385i");
            }

            // plot frs time axis once and for all.
            double XP=XBegin+(86-40)/5*XINC,YP=YBegin-(86-40)%5*YINC;
            GMT::MoveReferencePoint(outfile,"-Xf"+to_string(XP)+"i -Yf"+to_string(YP)+"i");
            GMT::set("FONT_ANNOT 4p");
            GMT::set("MAP_TICK_PEN_PRIMARY 0.4p,black");
            GMT::psbasemap(outfile,"-J -R -Bxa10 -Bya0.5 -BWS -O -K");


            // texts on Model properties.
            XP=XBegin+(90-40)/5*XINC,YP=YBegin-(90-40)%5*YINC;
            GMT::MoveReferencePoint(outfile,"-Xf"+to_string(XP)+"i -Yf"+to_string(YP)+"i");

            string line;
            double textYPosition=1;
            plotTexts.clear();
            ifstream fpin(indexDIR+"/"+eq);
            while (getline(fpin,line))
                plotTexts.push_back(GMT::Text(0,textYPosition-=0.5,line+"  M_"+to_string(Cnt),8,"LB"));
            fpin.close();

            GMT::pstext(outfile,plotTexts,"-J -R -N -O -K");
            GMT::SealPlot(outfile);
        }

        ShellExec("cat `ls -rt *_"+phase+".ps` > "+phase+".ps");
        ShellExec("ps2pdf "+phase+".ps");
        ShellExec("rm *_"+phase+".ps "+phase+".ps");
    }

    ShellExec("rm tmp.cpt");

    return 0;
}

