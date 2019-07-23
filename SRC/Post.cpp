#include<stdio.h>
#include<iostream>
#include<fstream>
#include<vector>
#include<string>
#include<cmath>
#include<algorithm>

using namespace std;

int main(int argc, char **argv){

	if (argc!=2) return 1;
    string FileList(argv[1]);

	// read infiles
	vector<string> filenames,fileRX,fileRZ;
	ifstream fpin(FileList);
	string tmpstr,n1,n2;
	while (fpin >> tmpstr >> n1 >> n2) {
		filenames.push_back(tmpstr);
		fileRX.push_back(n1);
		fileRZ.push_back(n2);
	}
	fpin.close();


	// for each file.
	int iout,it,nx,ny;
	double dt,xmin,xmax,ymin,ymax,dx,dy;
	vector<int> NX,NY;
	vector<double> XMIN,XMAX,YMIN,YMAX,DX,DY;
	int n=filenames.size();

	for (int i=0;i<n;++i){
		fpin.open(filenames[i]);

		// read meta data.
		fpin >> iout >> dt >> it >> ny >> ymin >> ymax 
	     >> nx >> xmin >> xmax;
		fpin.close();

		xmin*=180/M_PI;
		xmax*=180/M_PI;
		ymin/=1000;
		ymax/=1000;
		dx=(xmax-xmin)/(nx-1);
		dy=(ymax-ymin)/(ny-1);

		// adjust gird according to iout.
		// output for plotting.
		nx=1+(int)floor((xmax-xmin)/(iout*dx));
		ny=1+(int)floor((ymax-ymin)/(iout*dy));
		dx*=iout;
		dy*=iout;
		xmax=xmin+(nx-1)*dx;
		ymax=ymin+(ny-1)*dy;

		// Store meta data for further use.
		NX.push_back(nx);
		NY.push_back(ny);
		DX.push_back(dx);
		DY.push_back(dy);
		XMIN.push_back(xmin);
		XMAX.push_back(xmax);
		YMIN.push_back(ymin);
		YMAX.push_back(ymax);
	}

	// read values.
	// Find global max value.

	vector<vector<double>> data(n);
	double Max=-1;
	for (int i=0;i<n;++i){
		fpin.open(filenames[i]);

		// read meta data.
		fpin >> iout >> dt >> it >> ny >> ymin >> ymax 
	     >> nx >> xmin >> xmax;

		// read grid value.
		double val;
		while (fpin >> val) {
			data[i].push_back(val);
			Max=max(Max,fabs(val));
		}

		fpin.close();
	}


	// Rescale.
	for (int i=0;i<n;++i){
		for (size_t j=0;j<data[i].size();++j){
			data[i][j]=(data[i][j]>0?1:-1)*pow(fabs(data[i][j]/Max),0.33);
// 			data[i][j]=(data[i][j]>0?1:-1)*pow(fabs(data[i][j]/Max),0.1);
		}
	}


	// Resample with linear interpolation.
	// To join grids from every rank into a big grid.
	double xmin1=*min_element(XMIN.begin(),XMIN.end());
	double xmax1=*max_element(XMAX.begin(),XMAX.end());
	double ymin1=*min_element(YMIN.begin(),YMIN.end());
	double ymax1=*max_element(YMAX.begin(),YMAX.end());

	// adjust this big grid.
	int nx1=1+(int)floor((xmax1-xmin1)/dx);
	int ny1=1+(int)floor((ymax1-ymin1)/dy);
	xmax1=xmin1+(nx1-1)*dx;
	ymax1=ymin1+(ny1-1)*dy;

	FILE *fpout=fopen("tmp.desc","w");
	fprintf(fpout,"%.8lf\t%.8lf\t%.8lf\t%.8lf\t%.8lf\t%.8lf\n"
                 ,-xmax1,-xmin1,dx,ymin1,ymax1,dy);
	fclose(fpout);

	// interpolate and output.
	fpout=fopen("tmp.dat","w");
	for (int i=0;i<nx1;++i){
		for (int j=0;j<ny1;++j){
			double x=xmin1+i*dx,y=ymax1-j*dy;

			// find the rank containing this grid piont.
			bool flag=true;
			size_t I;
			for (I=0;I<XMIN.size();++I){
				if (XMIN[I]<=x && x<=XMAX[I] && YMIN[I]<=y && y<=YMAX[I]){
					flag=false;
					break;
				}
			}

			if (flag) {
				cerr << "Can't find the rank containing current grid point !"
				     << endl;
				return 1;
			}

			// Do interpolation.
			double v1,v2,v3,v4,vc1,vc2,vc;
			int Ix=(x-XMIN[I])/dx,Iy=(YMAX[I]-y)/dy;
			v1=data[I][Ix*NY[I]+Iy];
			v2=data[I][Ix*NY[I]+Iy+1];
			v3=data[I][Ix*NY[I]+Iy+NY[I]];
			v4=data[I][Ix*NY[I]+Iy+NY[I]+1];

			vc1=v1+(v2-v1)*(YMAX[I]-Iy*dy-y)/dy;
			vc2=v3+(v4-v3)*(YMAX[I]-Iy*dy-y)/dy;

			vc=vc1+(vc2-vc1)*(x-(XMIN[I]+Ix*dx))/dx;

			fprintf(fpout,"%.8e\t%.8e\t%.8e\n",-x,y,vc);
		}
	}
	fclose(fpout);

	return 0;
}
