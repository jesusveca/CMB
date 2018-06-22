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

#define ANSI_COLOR_RED     "\x1b[31m"
float acumulador;
int _pos_u, _pos_v, _indicator=0, max_iteraciones, _iter;
double _time_concurrent=0;
vector<vector<float>> numeros, matrix_dist_inicial, matrix_distancias, copia;
vector<string> _genes_iniciales, _genes_finales, headers;
string _ph, str;

vector<vector<string> > _genes_clusteres;


void _from_txt_to_vect(){
    ifstream _file;
    _file.open(("genes.csv"));
    while (!_file.eof()) {
        getline(_file,_ph);
        if(_indicator!=0){
            istringstream _istream(_ph);
            int _header_ind=0;  string _key;    vector<float> _tmp_vec;
            while(getline(_istream,str,',')){
                if(_header_ind==0)  _key=str;
                else    _tmp_vec.push_back(stod(str));
                _header_ind++;
            }
            numeros.push_back(_tmp_vec);
            _genes_iniciales.push_back(_key);
            vector<string> temp;
            temp.push_back(_key);
            _genes_clusteres.push_back(temp);

        }
        _indicator++;
    }
    _file.close();
}

void print_matrix_distances(string mensaje){
  // printf ("\033[0m%s\n-------------------------------------------\n",mensaje.c_str());
  for(int i=0;i<_genes_clusteres.size();i++){
    cout<<"Cluster "<<i+1<<"    ";
      for(int j=0;j<_genes_clusteres[i].size();j++)
        cout<<_genes_clusteres[i][j]<<" ";
      cout<<endl;
  }
  // for (int i=0;i<matrix_dist_inicial.size();i++){
  //   for(int j=0;j<matrix_dist_inicial.size();j++){
  //     printf ("\033[0m%5.2f ", matrix_dist_inicial[i][j]);
  //   }
  //   cout<<endl;
  // }
  // printf ("\033[0m    \n-------------------------------------------\n");

}

int pos_u=0;
int pos_v=0;
float menor_ ;

void _get_minimum(){
    float menor=400;
    for(int i=0;i<matrix_dist_inicial.size();i++){
        for  (int j=0;j<matrix_dist_inicial[i].size();j++){
            if(i!=j){
                if(matrix_dist_inicial[i][j]<menor){
                    menor=matrix_dist_inicial[i][j];
                    pos_u=i;
                    pos_v=j;
                    menor_=menor;
                }
            }
        }
    }
}

void _build_data_to_distance_matrix(){
    vector< vector<float> > matrix_numeros(_genes_iniciales.size()-1,vector<float>(_genes_iniciales.size()-1));
    for(int i=1; i< _genes_iniciales.size()-1; i++){
        for(int j=0; j<i; j++){
            for(int k=0; k<numeros[i].size();k++){
                if(j!=i)
                  acumulador+=((numeros[j][k]-numeros[i][k])*(numeros[j][k]-numeros[i][k]));
            }
            if (j!=i){
                float _calculus=sqrt(acumulador);   acumulador=0;
                matrix_numeros[i][j]=matrix_numeros[j][i]=_calculus;
            }
        }
    }
    matrix_dist_inicial=matrix_numeros;
    // print_matrix_distances("DISTANCIAS");
}

float _f_promedio(float a, float b){
    return (a+b)/2.0;
}

vector<vector<string> > clusteres;
void _clustering_action_uno(){
    _get_minimum();
    int _pos_a = pos_u;     int _pos_b = pos_v;
    if(_pos_b>_pos_a)   swap(_pos_b,_pos_a);

    vector<string > tmp;
    tmp.push_back(_genes_clusteres[_pos_a][0]);     tmp.push_back(_genes_clusteres[_pos_b][0]);

    _genes_clusteres[_pos_a]=tmp;       _genes_clusteres[_pos_b]=tmp;
    _genes_clusteres.erase(_genes_clusteres.cbegin()+_pos_b);


    for (int i=0;i<_genes_clusteres.size();i++)
        matrix_dist_inicial[_pos_a][i]=matrix_dist_inicial[i][_pos_a]=matrix_dist_inicial[i][_pos_b]=matrix_dist_inicial[_pos_b][i]=_f_promedio(matrix_dist_inicial[_pos_b][i],matrix_dist_inicial[_pos_a][i]);

    for (int i=0;i<matrix_dist_inicial.size();i++){
        vector<float> _tmp_vec;
        if(i!=_pos_a){
          for(int j=0;j<matrix_dist_inicial[i].size();j++){
              if(j!=_pos_b)   _tmp_vec.push_back(matrix_dist_inicial[i][j]);
          }
          copia.push_back(_tmp_vec);
        }
    }

    matrix_dist_inicial = copia;
    copia.clear();

    for (int i=0;i<matrix_dist_inicial.size();i++){ // para que el triangular superior sea igual al triangular inferior
        for(int j=0;j< matrix_dist_inicial[i].size();j++){
            if(i!=j)    matrix_dist_inicial[j][i]=matrix_dist_inicial[i][j];
            else        matrix_dist_inicial[i][i]=0;
        }
    }

}



void _clustering_action_dos(){
    _get_minimum();
    int _pos_a = pos_u;     int _pos_b = pos_v;
    if(_pos_b>_pos_a)   swap(_pos_b,_pos_a);

    // cout<<"minimos son : "<<menor_<<"   "<<pos_u<<"   "<<pos_v<<endl;

    vector<string > tmp;
    for(int i=0;i<_genes_clusteres[_pos_a].size();i++){
        tmp.push_back(_genes_clusteres[_pos_a][i]);
    }
    for(int j=0;j<_genes_clusteres[_pos_b].size();j++){
      tmp.push_back(_genes_clusteres[_pos_b][j]);
    }

    _genes_clusteres[_pos_a]=tmp;
    _genes_clusteres[_pos_b]=tmp;
    _genes_clusteres.erase(_genes_clusteres.cbegin()+_genes_clusteres.size()-1);

    float sum=0;
    float prom=0;
    for(int i=0;i<_genes_clusteres.size();i++){
        for(int j=0;j<_genes_clusteres[i].size();j++){
            int enum_=0;
            if((i!=_pos_b and i!=_pos_a)){
              for(int p=0;p<_genes_clusteres[_pos_a].size();p++){
                int pos_nido_uno = find(_genes_iniciales.begin(), _genes_iniciales.end(), _genes_clusteres[i][j]) - _genes_iniciales.begin();
                int pos_nido_dos = find(_genes_iniciales.begin(), _genes_iniciales.end(), _genes_clusteres[_pos_a][p]) - _genes_iniciales.begin();
                sum += matrix_distancias[pos_nido_dos][pos_nido_uno];
                enum_++;
              }
              prom = sum/enum_;
              matrix_dist_inicial[i][_pos_b] = prom;
              matrix_dist_inicial[i][_pos_a] = prom;
              matrix_dist_inicial[_pos_a][i] = prom;
              matrix_dist_inicial[_pos_b][i] = prom;
              sum=0;
              prom=0;
            }
        }
    }


    for (int i=0;i<matrix_dist_inicial.size();i++){
        vector<float> _tmp_vec;
        if(i!=_pos_a){
          for(int j=0;j<matrix_dist_inicial[i].size();j++){
              if(j!=_pos_a)   _tmp_vec.push_back(matrix_dist_inicial[i][j]);
          }
          copia.push_back(_tmp_vec);
        }
    }

    matrix_dist_inicial = copia;
    copia.clear();
    //
    for (int i=0;i<matrix_dist_inicial.size();i++){ // para que el triangular superior sea igual al triangular inferior
        for(int j=0;j< matrix_dist_inicial[i].size();j++){
             if(i!=j)   matrix_dist_inicial[j][i]=matrix_dist_inicial[i][j];
            else        matrix_dist_inicial[i][i]=0;
        }
    }
    _genes_clusteres.erase(_genes_clusteres.cbegin()+_pos_a);
}




void _clustering_action_dos_punto_cero(){
    _get_minimum();
    int _pos_a = pos_u;     int _pos_b = pos_v;
    if(_pos_b>_pos_a)   swap(_pos_b,_pos_a);

    // cout<<"minimos son : "<<menor_<<"   "<<pos_u<<"   "<<pos_v<<endl;

    vector<string > tmp;
    for(int i=0;i<_genes_clusteres[_pos_a].size();i++){
        tmp.push_back(_genes_clusteres[_pos_a][i]);
    }

    _genes_clusteres[_pos_a]=tmp;
    _genes_clusteres[_pos_b]=tmp;
    _genes_clusteres.erase(_genes_clusteres.cbegin()+_genes_clusteres.size()-1);

    float sum=0;
    float prom=0;
    for(int i=0;i<_genes_clusteres.size();i++){
        for(int j=0;j<_genes_clusteres[i].size();j++){
            int enum_=0;
            if((i!=_pos_b and i!=_pos_a)){
              for(int p=0;p<_genes_clusteres[_pos_a].size();p++){
                int pos_nido_uno = find(_genes_iniciales.begin(), _genes_iniciales.end(), _genes_clusteres[i][j]) - _genes_iniciales.begin();
                int pos_nido_dos = find(_genes_iniciales.begin(), _genes_iniciales.end(), _genes_clusteres[_pos_a][p]) - _genes_iniciales.begin();
                sum += matrix_distancias[pos_nido_dos][pos_nido_uno];
                enum_++;
              }
              prom = sum/enum_;
              matrix_dist_inicial[i][_pos_b] = prom;
              matrix_dist_inicial[i][_pos_a] = prom;
              matrix_dist_inicial[_pos_a][i] = prom;
              matrix_dist_inicial[_pos_b][i] = prom;
              sum=0;
              prom=0;
            }
        }
    }


    for (int i=0;i<matrix_dist_inicial.size();i++){
        vector<float> _tmp_vec;
        if(i!=_pos_a){
          for(int j=0;j<matrix_dist_inicial[i].size();j++){
              if(j!=_pos_a)   _tmp_vec.push_back(matrix_dist_inicial[i][j]);
          }
          copia.push_back(_tmp_vec);
        }
    }

    matrix_dist_inicial = copia;
    copia.clear();
    //
    for (int i=0;i<matrix_dist_inicial.size();i++){ // para que el triangular superior sea igual al triangular inferior
        for(int j=0;j< matrix_dist_inicial[i].size();j++){
             if(i!=j)   matrix_dist_inicial[j][i]=matrix_dist_inicial[i][j];
            else        matrix_dist_inicial[i][i]=0;
        }
    }
    _genes_clusteres.erase(_genes_clusteres.cbegin()+_pos_a);
}



void _cluster_decidir_tamano(){
    _get_minimum();
    // cout<<"minimos son : "<<menor_<<"   "<<pos_u<<"   "<<pos_v<<endl;
    //
    // cout<<"size of x:"<<_genes_clusteres[pos_u].size()<<endl;
    // cout<<"size of y:"<<_genes_clusteres[pos_v].size()<<endl;

    if(_genes_clusteres[pos_u].size()==1 or _genes_clusteres[pos_v].size()==1){
      // if(_genes_clusteres[pos_v].size()==1){
        _clustering_action_uno();
      // }
    }
    if(_genes_clusteres[pos_u].size()>1 or _genes_clusteres[pos_v].size()>1){
      // if(_genes_clusteres[pos_v].size()==1){
        // cout<<"1er if aqui\n";
        // _clustering_action_dos_punto_cero();
      // }
      // else {
        // cout<<"2do if aqui\n";
        _clustering_action_dos();
      // }

    }


}


void _clusterization(){
    _get_minimum();
    // cout<<"minimos son : "<<menor_<<"   "<<pos_u<<"   "<<pos_v<<endl;
    _clustering_action_uno(); // for the first action, when you only compare exactly 2 genes
    // print_matrix_distances("h");

    int contador =0;
    int jk=0;
    // contador--;
    jk = _genes_clusteres.size();
    cout<<jk<<endl;
    // while(contador>0){
    while(contador<100000){
      _cluster_decidir_tamano();
      // print_matrix_distances("h");
      contador++;
      jk = _genes_clusteres.size();
      cout<<jk<<endl;
      if(jk==670)
      break;
    }
    print_matrix_distances("h");

}




int main(){

    _from_txt_to_vect();
    _build_data_to_distance_matrix();
    matrix_distancias = matrix_dist_inicial;
    _genes_finales  = _genes_iniciales;
    _clusterization();
    return 0;
}
