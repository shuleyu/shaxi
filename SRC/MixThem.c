#include<stdio.h>
#include<stdlib.h>
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

//     enum PIenum {NPTS};
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

	FILE   *fpin,*fpout;
	int    Cnt2,NPTS,Size;
	double *p,tmpvar,**res;

	// Read input.
	fpin=fopen(PS[infile],"r");
	NPTS=0;
	while (fscanf(fpin,"%lf",&tmpvar)==1){
		++NPTS;
	}
	p=(double *)malloc(NPTS*sizeof(double));
	fclose(fpin);

	fpin=fopen(PS[infile],"r");
	NPTS=0;
	while (fscanf(fpin,"%lf",&p[NPTS])==1){
		++NPTS;
	}
	fclose(fpin);

	// Count how many models.
	Size=mixsize(p,NPTS);

	res=(double **)malloc(Size*sizeof(double *));
	for (Cnt=0;Cnt<Size;++Cnt){
		res[Cnt]=(double *)malloc(NPTS/3*sizeof(double));
	}

	// Mix them.
	mixthem(p,NPTS,res);

	// Output.
	fpout=fopen(PS[outfile],"w");
	for (Cnt=0;Cnt<Size;++Cnt){
		for (Cnt2=0;Cnt2<NPTS/3;++Cnt2){
			fprintf(fpout,"%lf\t",res[Cnt][Cnt2]);
		}
		fprintf(fpout,"\n");
	}

	fclose(fpout);

    // Free spaces.
	free_parameters(PI,PS,P,string_num);
	for (Cnt=0;Cnt<Size;++Cnt){
		free(res[Cnt]);
	}
	free(res);
	free(p);
	

    return 0;
}

