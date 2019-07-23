#include<iostream>
#include<fstream>
#include<sstream>
#include<cstdio>
#include<cstdlib>
#include<vector>
#include<string>
#include<cstring>
#include<cmath>
#include<random>
#include<chrono>
#include<unistd.h>
#include<algorithm>
extern "C"{
#include<sac.h>
#include<sacio.h>
#include<ASU_tools.h>
}

using namespace std;

inline int Window(double L,double delta){
	return (int)ceil(L/delta);
}

int main(int argc, char **argv){

    enum PIenum{AddNoise,FLAG1};
    enum PSenum{SynList,NoiseList,TargetList,OriginalSyn,FLAG2};
    enum Penum{synDelta,F1,F2,FLAG3};

    /****************************************************************

				Deal with inputs. (Store them in PI,PS,P)

    ****************************************************************/

	if (argc!=4){
		cerr << "In C++: Argument Error!" << endl;
		return 1;
	}

    int int_num,string_num,double_num;

    vector<int> PI;
    vector<string> PS;
    vector<double> P;

    int_num=atoi(argv[1]);
    string_num=atoi(argv[2]);
    double_num=atoi(argv[3]);

	if (FLAG1!=int_num){
		cerr << "In C++: Ints Naming Error !" << endl;
	}
	if (FLAG2!=string_num){
		cerr << "In C++: Strings Naming Error !" << endl;
	}
	if (FLAG3!=double_num){
		cerr << "In C++: Doubles Naming Error !" << endl;
	}

	string tmpstr;
	int tmpint,Cnt;
	double tmpval;

	Cnt=0;
	while (getline(cin,tmpstr)){
		++Cnt;
		stringstream ss(tmpstr);
		if (Cnt<=int_num){
			if (ss >> tmpint && ss.eof()){
				PI.push_back(tmpint);
			}
			else{
				cerr << "In C++: Ints reading Error !" << endl;
				return 1;
			}
		}
		else if (Cnt<=int_num+string_num){
			PS.push_back(tmpstr);
		}
		else if (Cnt<=int_num+string_num+double_num){
			if (ss >> tmpval && ss.eof()){
				P.push_back(tmpval);
			}
			else{
				cerr << "In C++: Doubles reading Error !" << endl;
				return 1;
			}
		}
		else{
			cerr << "In C++: Redundant inputs !" << endl;
			return 1;
		}
	}
	if (Cnt!=int_num+string_num+double_num){
		cerr << "In C++: Not enough inputs !" << endl;
		return 1;
	}

    /****************************************************************

                              Job begin.

    ****************************************************************/


	// Read in noise file names.
	vector<string> NoiseFileNames;
	ifstream fpin;
	fpin.open(PS[NoiseList]);
	while (fpin >> tmpstr) NoiseFileNames.push_back(tmpstr);
	fpin.close();


	// Chose a random noise for each trace.
    auto seed=chrono::system_clock::now().time_since_epoch().count();
    default_random_engine e(seed);
    uniform_int_distribution<int> NoiseIndex(0,NoiseFileNames.size()-1);


	// Read in noises.
	vector<vector<double>> Noises;

	int maxl=2000000,rawnpts,nerr,MaxNoiseLength=0;
	float rawbeg,rawdel;
	float *amp=(float *)malloc(maxl*sizeof(float));
	char tmpfilename[200];

	for (auto noisefile: NoiseFileNames){

		// Read in noise SAC.
		strcpy(tmpfilename,noisefile.c_str());
		rsac1(tmpfilename,amp,&rawnpts,&rawbeg,&rawdel,&maxl,&nerr,noisefile.size());
		normalize(amp,rawnpts);


		// Taper this noise so that it's reasonable to periodically when wrapping around.
		taper(amp,rawnpts,10.0/rawnpts);

		// Interpolate it to synthetic sampling rate.
		float *time1=(float *)malloc(rawnpts*sizeof(float));
		for (int i=0;i<rawnpts;++i) time1[i]=rawbeg+rawdel*i;

		int new_npts=(int)floor((rawnpts-1)*rawdel/P[synDelta]);
		double *amp2=(double *)malloc(new_npts*sizeof(double));
		double *time2=(double *)malloc(new_npts*sizeof(double));
		for (int i=0;i<new_npts;++i) time2[i]=rawbeg+P[synDelta]*i;

		wiginterp_f(time1,amp,rawnpts,time2,amp2,new_npts,0);
		free(time1);
		free(time2);

		// Push it into the noise set.
		vector<double> NoiseTrace(new_npts);
		for(int i=0;i<new_npts;++i) NoiseTrace[i]=amp2[i];
		free(amp2);

		Noises.push_back(NoiseTrace);
		MaxNoiseLength=max(MaxNoiseLength,rawnpts);
	}

	// For each signal, assign a random start ponit.
	uniform_int_distribution<int> NoiseBegin(0,MaxNoiseLength-1);

	// Main stuff.
	//
	// For each data, select the similar distance synthetics,
	// chose a similar noise level,
	// periodically create the noise trace,
	// add it to synthetics,
	// output both noise and noise added synthetics.

	ifstream synlist(PS[SynList]);
	vector<string> synFile;
	vector<vector<double>> Syn;
	vector<double> synGcarc;
	vector<int> synS,synPeakScS;
	int tmpint2;
	double synBegin;

	while (synlist >> tmpstr >> tmpval >> tmpint >> tmpint2){
		synFile.push_back(tmpstr);
		synGcarc.push_back(tmpval);
		synS.push_back(tmpint);
		synPeakScS.push_back(tmpint2);

		// Read in synthetic traces.
		strcpy(tmpfilename,tmpstr.c_str());
		rsac1(tmpfilename,amp,&rawnpts,&rawbeg,&rawdel,&maxl,&nerr,tmpstr.size());
		synBegin=rawbeg;
		vector<double> SynTrace(rawnpts);
		for(int index=0;index<rawnpts;index++) SynTrace[index]=amp[index];
		Syn.push_back(SynTrace);
	}
	synlist.close();

	ifstream targetlist(PS[TargetList]);
	ofstream fpout(PS[OriginalSyn]);
	string EQ,stnm;
	double SNR_ScS,Shift_Gcarc,stlo,stla,evlo,evla,hitlo,hitla;
	Cnt=1;
	while (targetlist >> Shift_Gcarc >> SNR_ScS
					  >> evlo >> evla >> stlo >> stla >> EQ >> stnm
					  >> hitlo >> hitla){

		// Select the correct gcarc synthetics.
		auto it=lower_bound(synGcarc.begin(),synGcarc.end(),Shift_Gcarc);
		if (it==synGcarc.end()) continue;
		if (it==synGcarc.begin() && fabs(*it-Shift_Gcarc)>0.1) continue;
		if (fabs(*prev(it)-Shift_Gcarc)>0.1 && fabs(*next(it)-Shift_Gcarc)>0.1) continue;
		int Index=distance(synGcarc.begin(),it);
		if (fabs(synGcarc[Index-1]-Shift_Gcarc)<fabs(synGcarc[Index]-Shift_Gcarc)) --Index;

		// Wrap around the random selected noise trace.
		int Index2=NoiseIndex(e);
		vector<double> NoiseTrace=Noises[Index2];
		while (NoiseTrace.size()<Syn[Index].size())
			NoiseTrace.insert(NoiseTrace.end(),Noises[Index2].begin(),Noises[Index2].end());
		NoiseTrace.resize(Syn[Index].size());

		int RandomBegin=NoiseBegin(e);

		// Add the proper noise level.
		double l=0,r=1,mid=0,NewSNR=-10;
		double *NoiseAdded=(double *)malloc(Syn[Index].size()*sizeof(double));
		int Cnt2=0;

		if (PI[AddNoise]!=0){
			while (fabs(NewSNR-SNR_ScS)>0.05) {

				mid=(l+r)/2;
				for (size_t i=0;i<Syn[Index].size();++i)
					NoiseAdded[i]=Syn[Index][i]+mid*NoiseTrace[(i+RandomBegin)%NoiseTrace.size()];

				butterworth_bp(&NoiseAdded,1,Syn[Index].size(),P[synDelta],2,2,P[F1],P[F2],&NoiseAdded);
				// Find the new peak on ScS trace.
				int NewPeak=synPeakScS[Index];
				findpeak(NoiseAdded,Syn[Index].size(),&NewPeak,Window(-5,P[synDelta]),Window(10,P[synDelta]));

				NewSNR=snr_envelope(NoiseAdded,Syn[Index].size(),
									synS[Index]-Window(180,P[synDelta]),Window(120,P[synDelta]),
									NewPeak-Window(5,P[synDelta]),Window(11,P[synDelta]));
				if (NewSNR<SNR_ScS) r=mid;
				else l=mid;
				if (Cnt2++>50) break;
			}
		}

		for (size_t i=0;i<Syn[Index].size();++i)
			NoiseAdded[i]=Syn[Index][i]+mid*NoiseTrace[(i+RandomBegin)%NoiseTrace.size()];

		string outfile;

		// Output noise.
// 		outfile=to_string(Cnt)+".Noise.sac";
// 		strcpy(tmpfilename,outfile.c_str());
// 		rawnpts=NoiseTrace.size();
// 		rawdel=P[synDelta];
// 		rawbeg=synBegin;
// 		for (int i=0;i<rawnpts;++i) amp[i]=NoiseTrace[i];
// 		wsac1(tmpfilename,amp,&rawnpts,&rawbeg,&rawdel,&nerr,outfile.size());

		// Output noise added synthetics.
		outfile=to_string(Cnt)+".sac";
		strcpy(tmpfilename,outfile.c_str());
		rawnpts=Syn[Index].size();
		rawdel=P[synDelta];
		rawbeg=synBegin;
		for (int i=0;i<rawnpts;++i) amp[i]=NoiseAdded[i];
		wsac1(tmpfilename,amp,&rawnpts,&rawbeg,&rawdel,&nerr,outfile.size());

		// Output info.
		//
		waypoint(stlo,stla,evlo,evla,Shift_Gcarc,&evlo,&evla);

		fpout << synFile[Index] << " " << Shift_Gcarc << " "
		      << evlo << " " << evla << " " << stlo << " " << stla
			  << " " << EQ << " " << stnm << " " << hitlo << " " << hitla << endl;

		++Cnt;
		free(NoiseAdded);
	}
	targetlist.close();
	fpout.close();

	free(amp);

    return 0;
}
