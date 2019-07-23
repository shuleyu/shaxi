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

    enum PIenum {num_theta,num_r};
    enum PSenum {infile_theta,infile_r,infile_vs,infile_rho,outfile_vs,outfile_rho,outfile_vs_zoom,outfile_rho_zoom};
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

	FILE   *fpin, *fpout, *fpin2, *fpout2,*fpout3,*fpout4;
	int    index1,index2;
	double *theta,*r,vs,rho,dvs,drho;

	theta=(double *)malloc(PI[num_theta]*sizeof(double));
	r=(double *)malloc(PI[num_r]*sizeof(double));

	// Read in data.
	fpin=fopen(PS[infile_theta],"r");
	for (Cnt=0;Cnt<PI[num_theta];++Cnt){
		fscanf(fpin,"%lf,",&theta[Cnt]);
	}
	fclose(fpin);

	fpin=fopen(PS[infile_r],"r");
	for (Cnt=0;Cnt<PI[num_r];++Cnt){
		fscanf(fpin,"%lf,",&r[Cnt]);
	}
	fclose(fpin);
	

	// Output grid file.
	fpin=fopen(PS[infile_vs],"r");
	fpin2=fopen(PS[infile_rho],"r");
	fpout=fopen(PS[outfile_vs],"w");
	fpout2=fopen(PS[outfile_rho],"w");
	fpout3=fopen(PS[outfile_vs_zoom],"w");
	fpout4=fopen(PS[outfile_rho_zoom],"w");
	for (Cnt=0;Cnt<PI[num_r]*PI[num_theta];++Cnt){
		fscanf(fpin,"%lf,",&vs);
		fscanf(fpin2,"%lf,",&rho);

		index2=Cnt%PI[num_r];
		index1=Cnt/PI[num_r];

		if (r[index2]<=6371.0){
			if (r_vs(r[index2])>0){
				dvs=(vs/r_vs(r[index2])-1)*100;
			}
			else{
				dvs=0;
			}
			drho=(rho/r_rho(r[index2])-1)*100;
		}
		else{
			dvs=0;
			drho=0;
		}

		fprintf(fpout,"%.4lf\t%.4lf\t%.2lf\n",theta[index1],r[index2],dvs);
		fprintf(fpout2,"%.4lf\t%.4lf\t%.2lf\n",theta[index1],r[index2],drho);

		fprintf(fpout3,"%.4lf\t%.4lf\t%.2lf\n",theta[index1],r[index2]-3480.0,dvs);
		fprintf(fpout4,"%.4lf\t%.4lf\t%.2lf\n",theta[index1],r[index2]-3480.0,drho);

	}
	fclose(fpin);
	fclose(fpin2);
	fclose(fpout);
	fclose(fpout2);
	fclose(fpout3);
	fclose(fpout4);


    // Free spaces.
	free_parameters(PI,PS,P,string_num);
	free(theta);
	free(r);

    return 0;
}

