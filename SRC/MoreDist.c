#include<stdio.h>
#include<stdlib.h>
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

    enum PIenum {NPTS_Old};
    enum PSenum {infile,outfile};
//     enum Penum  {};

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

	int    NPTS_New;
	FILE   *fpin,*fpout;
	double *dist,*center,dt,*dist_new,*center_new;

	dist=(double *)malloc(PI[NPTS_Old]*sizeof(double));
	center=(double *)malloc(PI[NPTS_Old]*sizeof(double));

	fpin=fopen(PS[infile],"r");
	for (Cnt=0;Cnt<PI[NPTS_Old];++Cnt){
		fscanf(fpin,"%lf%lf",&dist[Cnt],&center[Cnt]);
	}
	fclose(fpin);

	dt=0.01;
	NPTS_New=1+(int)floor(dist[PI[NPTS_Old]-1]-dist[0])/dt;

	dist_new=(double *)malloc(NPTS_New*sizeof(double));
	center_new=(double *)malloc(NPTS_New*sizeof(double));

	for (Cnt=0;Cnt<NPTS_New;++Cnt){
		dist_new[Cnt]=dist[0]+dt*Cnt;
	}

	wiginterpd(dist,center,PI[NPTS_Old],dist_new,center_new,NPTS_New,0);

	fpout=fopen(PS[outfile],"w");
	for (Cnt=0;Cnt<NPTS_New;++Cnt){
		fprintf(fpout,"%.3lf\t%.3lf\n",center_new[Cnt],dist_new[Cnt]);
	}
	fclose(fpout);

    // Free spaces.
	free(dist);
	free(center);
	free(dist_new);
	free(center_new);
	free_parameters(PI,PS,P,string_num);

    return 0;
}

