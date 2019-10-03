#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<netcdf.h>
#include<gmt/gmt.h>
#include<unistd.h>
#include<ASU_tools.h>

#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); return 2;}

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
    enum PSenum {PROJ,REG,Mark,INFILEs,cptfile,OUTFILE,RefFile};
    enum Penum  {MIN,MAX};

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

	int    ncid,varid,Theta_Size,R_Size,gmtid,Count;
	float  *Val,*Val_gmt,*theta,*r;
	double *PREM_Val,wesn[4],inc[2];
	void   *API;
	struct GMT_GRID *Grd;
	char   grdfile[100],command[200],filename[200];
	FILE   *fpin,*fpout;

	// Allocate memories.
	//Theta_Size=32512/filenr(PS[INFILEs])+8;
	Theta_Size=20000/filenr(PS[INFILEs])+8;
	R_Size=3808;

	PREM_Val=(double *)malloc(R_Size*sizeof(double));
	theta=(float *)malloc(Theta_Size*sizeof(float));
	r=(float *)malloc(R_Size*sizeof(float));
	Val=(float *)malloc(R_Size*Theta_Size*sizeof(float));
	Val_gmt=(float *)malloc(R_Size*Theta_Size*sizeof(float));

	// GMT_API and plot first line.
	API=GMT_Create_Session("My session",0,0,NULL);

	fpin=fopen(PS[INFILEs],"r");
	Count=0;
	while (fscanf(fpin,"%s",filename)==1){

		// Read in NC files.
		nc_open(filename,NC_NOWRITE,&ncid);

		nc_inq_varid(ncid,"theta",&varid);
		nc_get_var_float(ncid,varid,theta);

		if (theta[Theta_Size-1]<P[MIN] || theta[0]>P[MAX]){
			nc_close(ncid);
			continue;
		}

		nc_inq_varid(ncid,"r",&varid);
		nc_get_var_float(ncid,varid,r);

		nc_inq_varid(ncid,PS[Mark],&varid);
		nc_get_var_float(ncid,varid,Val);

		nc_close(ncid);

		// Get PREM value for this nc file.
		if (strcmp(PS[Mark],"Vs1")==0){
			for (Cnt=0;Cnt<R_Size;Cnt++){
				PREM_Val[Cnt]=r_vs(r[Cnt]);
			}
		}
		else if (strcmp(PS[Mark],"Rho")==0){
			for (Cnt=0;Cnt<R_Size;Cnt++){
				PREM_Val[Cnt]=r_rho(r[Cnt]);
			}
		}
        else {
			for (Cnt=0;Cnt<R_Size;Cnt++){
				PREM_Val[Cnt]=0;
			}
        }

		// Output reference model.
		if (Count==0) {
			fpout=fopen(PS[RefFile],"w");
			for (Cnt=0;Cnt<R_Size;Cnt++){
				fprintf(fpout,"%.2lf\t%.2lf\t%.2lf\n",6371.0-r[Cnt],PREM_Val[Cnt],Val[Cnt]);
			}
			fclose(fpout);
		}

		// Modify absolute val to relative val to PREM.
		// Pad for eps.
		for (Cnt=0;Cnt<R_Size*Theta_Size;Cnt++){
			if (PREM_Val[Cnt%R_Size]!=0){
				Val[Cnt]=(Val[Cnt]-PREM_Val[Cnt%R_Size])*100.0/PREM_Val[Cnt%R_Size]+0.001;
			}
		}

		// Modify data array for GMT API (which behaves very weird).
		redirect_f(Val,Theta_Size,R_Size,Val_gmt);
		if (P[MAX]<180){
			for (Cnt=0;Cnt<R_Size;Cnt++){
				r[Cnt]-=3480.0;
			}
		}
		reverse_array_f(r,R_Size);

		// GMT_plot.
		wesn[0]=theta[0];
		wesn[1]=theta[Theta_Size-1];
		wesn[2]=r[0];
		wesn[3]=r[R_Size-1];
		inc[0]=(theta[Theta_Size-1]-theta[0])/(Theta_Size-1);
		inc[1]=(r[R_Size-1]-r[0])/(R_Size-1);

		Grd=GMT_Create_Data(API,GMT_IS_GRID,GMT_IS_SURFACE,GMT_GRID_HEADER_ONLY,NULL,wesn,inc,0,0,NULL);

		Grd->data=Val_gmt;
		gmtid=GMT_Register_IO(API,GMT_IS_GRID,GMT_IS_REFERENCE,GMT_IS_SURFACE,GMT_IN,NULL,Grd);
		GMT_Encode_ID(API,grdfile,gmtid);
		sprintf(command,"%s %s -<%s -C%s -O -K ->>%s",PS[PROJ],PS[REG],grdfile,PS[cptfile],PS[OUTFILE]);
		GMT_Call_Module(API,"grdimage",GMT_MODULE_CMD,command);
		Grd->data=NULL;

		GMT_Destroy_Data(API,&Grd);
		Count++;
	}
	fclose(fpin);

	// GMT_plot, last line.
	GMT_Destroy_Session(API);
	unlink("gmt.history");

    // Free spaces.
	free(PREM_Val);
	free(theta);
	free(r);
	free(Val);
	free(Val_gmt);

	free_parameters(PI,PS,P,string_num);

	return 0;
}
