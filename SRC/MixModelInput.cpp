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

#include<ReadParameters.hpp>
#include<MeshGrid.hpp>
#include<Enumerate.hpp>

using namespace std;

int main(int argc, char **argv){

    enum PI{FLAG1};
    enum PS{ModelInput,ModelOutputPrefix,FLAG2};
    enum PF{FLAG3};
    auto P=ReadParameters<PI,PS,PF>(argc,argv,cin,FLAG1,FLAG2,FLAG3);

    /****************************************************************

                              Job begin.

    ****************************************************************/

    ifstream fpin(P[ModelInput]);
    vector<string> types;
    vector<size_t> Cnt;
    vector<vector<vector<double>>> grid;
    string one_line;
    while (getline(fpin,one_line)){
        stringstream ss(one_line);
        double left,right,inc;
        vector<vector<double>> data;
        string type;
        ss >> type;
        while (ss >> left){
            ss >> right >> inc;
            data.push_back({left,right,inc});
        }
        grid.push_back(MeshGrid(data,1));
        Cnt.push_back(grid.back().size());
        types.push_back(type);
    }
    fpin.close();

    auto res=Enumerate(Cnt);
    for (size_t i=0;i<res.size();++i) {
        ofstream fpout(P[ModelOutputPrefix]+to_string(1+i));
        for (size_t j=0;j<types.size();++j){
            fpout << types[j] << " ";
            for (const auto &item: grid[j][res[i][j]])
                fpout << item << " " ;
            fpout << '\n';
        }
        fpout.close();
    }

    return 0;
}
