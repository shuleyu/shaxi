#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
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
    enum PSenum {outfile};
    enum Penum  {EVDE,gcarc,depth};

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


    /****************************************************************

                              Job begin.

    ****************************************************************/

	int    npts,res;
    double **ans;
    char   Phase[10];
	FILE   *fpout;

    strcpy(Phase,"ScS");

	npts=15000;
	ans=(double **)malloc(npts*sizeof(double *));
	for (Cnt=0;Cnt<npts;Cnt++){
		ans[Cnt]=(double *)malloc(3*sizeof(double));
	}

    // Use function.
	res=waypoint_sectionpath(Phase,0.0,0.0,P[EVDE],0.0,P[gcarc],P[depth],npts,ans);

	if (res<0){
		exit(1);
	}

	fpout=fopen(PS[outfile],"w");
	for (Cnt=0;Cnt<res;Cnt++){
		fprintf(fpout,"%11.3lf%11.2lf\n",ans[Cnt][1],2891-ans[Cnt][2]);
	}
	fclose(fpout);

    // Free spaces.
	for (Cnt=0;Cnt<npts;Cnt++){
		free(ans[Cnt]);
	}
	free(ans);
	free_parameters(PI,PS,P,string_num);
    
    return 0;    
}
