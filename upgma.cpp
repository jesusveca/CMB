#include <iostream>
#include <map>
#include <fstream>
#include <math.h>
#include <sstream>
#include<set>
#include <algorithm>
#include <vector>
#include <map>

using namespace std;

unsigned _initial_sub, _initial;
float _minimo_general=5, acumulador;
int _pos_u,_pos_v,_indicator=0,max_iteraciones, _iter;
double _time_concurrent=0;
vector<vector<float>> numeros,matrix_numeros_final,copia;
vector<string> headers, _header_inds;
string tipo, _ph,str;

void _from_txt_to_vect(){
    ifstream _file;
    _file.open(("genes.csv"));
    while (!_file.eof()) {
        getline(_file,_ph);
        if(_indicator!=0){
            istringstream _istream(_ph);
            int _header_ind=0;
            string _key;
            vector<float> _tmp_vec;
            while(getline(_istream,str,',')){
                if(_header_ind==0)
                    _key=str;
                else
                    _tmp_vec.push_back(stod(str));
                _header_ind++;
            }
            numeros.push_back(_tmp_vec);
            _header_inds.push_back(_key);
        }
        _indicator++;
    }
    _file.close();
}


void _build_data_to_distance_matrix(){
    vector< vector<float> > matrix_numeros(_header_inds.size()-1,vector<float>(_header_inds.size()-1));
    for (int i=0;i<matrix_numeros.size();i++){ // to fill differents values from zero into diagonal
        for (int j=0;j<matrix_numeros.size();j++){    if(j==i)    matrix_numeros[i][j]=9999;} }

    for(int i=1; i< _header_inds.size()-1; i++){
            for(int j=0; j<i; j++){
                for(int k=0; k<numeros[i].size();k++){
                    if(j!=i)
                        acumulador+=((numeros[j][k]-numeros[i][k])*(numeros[j][k]-numeros[i][k]));
                }
                if (j!=i){
                    float _calculus=sqrt(acumulador);   acumulador=0;
                    matrix_numeros[i][j]=matrix_numeros[j][i]=_calculus;
                    if(_calculus<_minimo_general){
                        _minimo_general = _calculus;
                        _pos_v =i;  _pos_u =j;
                    }
                }
            }
    }
    matrix_numeros_final=matrix_numeros;

    // for (int i=0;i<matrix_numeros_final.size();i++){
    //     cout<<_header_inds[i]<<",";
    //     for(int j=0;j<matrix_numeros_final.size();j++){
    //         cout<<matrix_numeros_final[i][j]<<",";
    //     }
    //     cout<<endl;
    //     // cout<<endl;
    // }
    // cout<<endl; cout<<"\nminimo: "<<_minimo_general<<endl;  cout<<"coords   x:"<<_pos_v<<"     y:"<<_pos_u<<endl<<endl;
}

float _f_minimo(float a, float b){
    return min(a,b);
}
float _f_maximo(float a, float b){
    return max(a,b);
}

float _f_promedio(float a, float b){
    return (a+b)/2.0;
}

void clustering(){
    int _pos_a = _pos_v;
    int _pos_b = _pos_u;
    _header_inds[_pos_b]=_header_inds[_pos_a]=string(_header_inds[_pos_a])+","+string(_header_inds[_pos_b]);

    for (int i=0;i<_header_inds.size()-1;i++){
        if (tipo=="minimo") matrix_numeros_final[_pos_a][i]=matrix_numeros_final[i][_pos_a]=matrix_numeros_final[i][_pos_b]=matrix_numeros_final[_pos_b][i]=_f_minimo(matrix_numeros_final[_pos_b][i],matrix_numeros_final[_pos_a][i]);
        if (tipo=="maximo") matrix_numeros_final[_pos_a][i]=matrix_numeros_final[i][_pos_a]=matrix_numeros_final[i][_pos_b]=matrix_numeros_final[_pos_b][i]=_f_maximo(matrix_numeros_final[_pos_b][i],matrix_numeros_final[_pos_a][i]);
        if (tipo=="promedio") matrix_numeros_final[_pos_a][i]=matrix_numeros_final[i][_pos_a]=matrix_numeros_final[i][_pos_b]=matrix_numeros_final[_pos_b][i]=_f_promedio(matrix_numeros_final[_pos_b][i],matrix_numeros_final[_pos_a][i]);
    }

    //copiar a otro vector solo los que son necesarios;
    for (int i=0;i<matrix_numeros_final.size();i++){
        vector<float> _tmp_vec;
        if (i!=_pos_b){
            for(int j=0;j<matrix_numeros_final.size();j++){
                if(j!=_pos_a)
                    _tmp_vec.push_back(matrix_numeros_final[i][j]);
            }
            headers.push_back(_header_inds[i]);
            copia.push_back(_tmp_vec);
    }}
    _header_inds=headers;
    matrix_numeros_final=copia;
    headers.clear();
    copia.clear();

    _minimo_general=100000;
	for(int i=0; i<_header_inds.size(); i++){
		for(int j=0; j<i; j++){
            matrix_numeros_final[i][j]=matrix_numeros_final[j][i];
            if (j!=i){
                if (matrix_numeros_final[i][j]<_minimo_general){
                    _minimo_general = matrix_numeros_final[i][j];
                    _pos_v =i;
                    _pos_u =j;
                }
            }
        }
    }

    // cout<<endl;
    // cout<<endl;
    // for (int i=0;i<matrix_numeros_final.size();i++){
    //     cout<<_header_inds[i]<<"  ";
    //     for(int j=0;j<matrix_numeros_final.size();j++){
    //         if (i==j)  cout<<0<<"   ";
    //         else    cout<<matrix_numeros_final[i][j]<<"   ";
    //     }
    //     cout<<endl;
    // }
    // cout<<endl;
    // cout<<endl;
    // cout<<"\nminimo: "<<_minimo_general<<endl;  cout<<"coords   x:"<<_pos_v<<"     y:"<<_pos_u<<endl;   cout<<" \n\n";
    _initial = clock(); double _res_time = (double(_initial-_initial_sub))/100000000;  _time_concurrent =_time_concurrent+_res_time;
}

int num_clusteres=0;

void _clusterization(){
    int n=0;
    while(n<10000){

      // cout<<endl;
      // // cout<<endl;
      // for (int i=0;i<matrix_numeros_final.size();i++){
      //     if (_header_inds[i].size()==7){
      //         cout<<"Cluster "<<i+1<<".- "<<_header_inds[i]<<" \n";
      //     }
      //     else {
      //         cout<<"Cluster "<<i+1<<".- "<<_header_inds[i]<<endl;
      //     }
      // }
      // cout<<endl;

        clustering();
        // cout<<"numero de clusteres: "<<_header_inds.size()<<endl;
        n=_header_inds.size();

        if(n==num_clusteres){
          cout<<endl << "\ntime :   " << _time_concurrent <<" secs"<< endl;
          break;
        }

        // n++;
    }
    cout<<endl;
    // cout<<endl;
    for (int i=0;i<matrix_numeros_final.size();i++){
        if (_header_inds[i].size()==7){
            cout<<"Cluster "<<i+1<<".- "<<_header_inds[i]<<" \n";
        }
        else {
            cout<<"Cluster "<<i+1<<".- "<<_header_inds[i]<<endl;
        }
    }
    cout<<endl;

}



int main(){
    tipo ="promedio";

    num_clusteres = 50;
	_from_txt_to_vect();
    _build_data_to_distance_matrix();
    _clusterization();
	return 0;
}
