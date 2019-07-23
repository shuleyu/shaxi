#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<sac.h>
#include<sacio.h>
#include<ASU_tools.h>

void free_parameters(int *PI, char **PS, double *P, int string_num){

	int Cnt;

	for (Cnt=0;Cnt<string_num;++Cnt){
		free(PS[Cnt]);
	}
	free(P);
	free(PI);
	free(PS);

	return;
}

int main(int argc, char **argv){
	if (argc!=4){
		printf("In C : Argument Error!\n");
		return 1;
	}

    // Deal with inputs.
    int    int_num,string_num,double_num,Cnt;
    int    *PI;
    char   **PS;
    double *P;

//     enum PIenum {};
    enum PSenum {infile};
    enum Penum  {Period,Length,delta};

    int_num=atoi(argv[1]);
    string_num=atoi(argv[2]);
    double_num=atoi(argv[3]);

    PI=(int *)malloc(int_num*sizeof(int));
    PS=(char **)malloc(string_num*sizeof(char *));
    P=(double *)malloc(double_num*sizeof(double));
    for (Cnt=0;Cnt<string_num;Cnt++){
        PS[Cnt]=(char *)malloc(200*sizeof(char));
    }

    for (Cnt=0;Cnt<int_num;Cnt++){
        if (scanf("%d",PI+Cnt)!=1){
            printf("In C : Int parameter reading Error !\n");
			free_parameters(PI,PS,P,string_num);
            return 1;
        }
    }

    for (Cnt=0;Cnt<string_num;Cnt++){
        if (scanf("%s",PS[Cnt])!=1){
            printf("In C : String parameter reading Error !\n");
			free_parameters(PI,PS,P,string_num);
            return 1;
        }
    }

    for (Cnt=0;Cnt<double_num;Cnt++){
        if (scanf("%lf",P+Cnt)!=1){
            printf("In C : Double parameter reading Error !\n");
			free_parameters(PI,PS,P,string_num);
            return 1;
        }
    }


//     for (Cnt=0;Cnt<int_num;Cnt++){
//             printf("%d\n",PI[Cnt]);
//     }
//     for (Cnt=0;Cnt<string_num;Cnt++){
//             printf("%s\n",PS[Cnt]);
//     }
//     for (Cnt=0;Cnt<double_num;Cnt++){
//             printf("%lf\n",P[Cnt]);
//     }
//     return 1;

    /****************************************************************

                              Job begin.

    ****************************************************************/

	float  *tmpdata,beg,del_float;
	double **data,*ptime,del_double;
	int    Cnt2,fileN,nptsy,*bad,nptsx,nerr;
	char   **in_list;
	FILE   *fpin;

	nptsx=filenr(PS[infile]);

	beg=0.0;

	nptsy=(int)floor(P[Length]/P[delta]);

	bad=(int *)malloc(nptsx*sizeof(int));
	ptime=(double *)malloc(nptsx*sizeof(double));
	tmpdata=(float *)malloc(nptsy*sizeof(float));

	data=(double **)malloc(nptsx*sizeof(double *));
	in_list=(char **)malloc(nptsx*sizeof(char *));
	for (Cnt=0;Cnt<nptsx;++Cnt){
		data[Cnt]=(double *)malloc(nptsy*sizeof(double));
		in_list[Cnt]=(char *)malloc(200*sizeof(char));
	}

	// Read in SAC.
	fpin=fopen(PS[infile],"r");
	for (Cnt=0;Cnt<nptsx;++Cnt){
		fscanf(fpin,"%s%lf",in_list[Cnt],&ptime[Cnt]);
	}
	fclose(fpin);

	del_double=P[delta];
	fileN=read_sac_fixdel(data,nptsx,&nptsy,ptime,0.0,&del_double,0.0,0.0,0,0,0,0,0.0,in_list,bad);
	del_float=del_double;

	if (fileN!=nptsx){
		for (Cnt=0;Cnt<nptsx;++Cnt){
			if (bad[Cnt]==1){
				printf("SAC file %s read error ...",in_list[Cnt]);
				exit(1);
			}
		}
	}

	// Convolve.
	gaussblur_1d(data,nptsx,nptsy,del_double,P[Period]/2/M_PI,50.0,data);

	// Write out SAC.
	for (Cnt=0;Cnt<nptsx;++Cnt){
		for (Cnt2=0;Cnt2<nptsy;++Cnt2){
			tmpdata[Cnt2]=data[Cnt][Cnt2];
		}
		wsac1(in_list[Cnt],tmpdata,&nptsy,&beg,&del_float,&nerr,strlen(in_list[Cnt]));
	}

    // Free spaces.
	free(bad);
	free(ptime);
	free(tmpdata);

	for (Cnt=0;Cnt<nptsx;++Cnt){
		free(data[Cnt]);
		free(in_list[Cnt]);
	}
	free(data);
	free(in_list);

	free_parameters(PI,PS,P,string_num);

    return 0;
}
