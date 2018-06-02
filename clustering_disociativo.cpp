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

unsigned _initial_sub, _initial;
float _minimo_general=5, acumulador;
int _pos_u,_pos_v,_indicator=0,max_iteraciones, _iter;
double _time_concurrent=0;
vector<vector<float>> numeros,matrix_numeros_final,copia;
vector<string> headers, _genes_iniciales;
string tipo, _ph,str;

vector<float> _vector_fill;
vector<vector<string> > _genes_clusteres;
vector<vector<float> > diferencias;

map <string, vector<float> > map_diferencias;
int num_clusteres_finales;

int num_clusters_here=0;


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
            numeros.push_back(_tmp_vec);    _genes_iniciales.push_back(_key);
        }
        _indicator++;
    }
    _file.close();
}

void print_genes_clusteres(){
    int cont=0;
    for(int i=0;i<_genes_clusteres.size();i++){
        cout<<"Cluster "<<i+1<<": ";
        for(int j=0;j<_genes_clusteres[i].size();j++){
            cout<<_genes_clusteres[i][j]<<", ";
        }
        cont++;
        cout<<endl;
    }
    cont++;
    cout<<"Cluster "<<cont<<": ";
    for(auto it2 = map_diferencias.cbegin(); it2 != map_diferencias.cend(); ++it2){
        cout<<it2->first<<", ";
    }
    cout<<endl;
}

void print_matrix_distances(string mensaje){
    printf ("\033[0m%s\n-------------------------------------------\n",mensaje.c_str());  
    printf ("\033[0m     ");  
    for(int i=0;i<_genes_iniciales.size();i++)
        printf (ANSI_COLOR_RED "%5.2s ", _genes_iniciales[i].c_str()); 
    printf ("\033[0m \n");  
    for (int i=0;i<matrix_numeros_final.size();i++){
        printf (ANSI_COLOR_RED "%5.2s ", _genes_iniciales[i].c_str()); 
        for(int j=0;j<matrix_numeros_final.size();j++){
            if(i!=j)    printf ("\033[0m%5.2f ", matrix_numeros_final[i][j]);  
            else        printf (ANSI_COLOR_RED "%5.2f ", matrix_numeros_final[i][j]);                
        }
        printf ("\033[0m \n");  
    }
    printf ("\033[0m    \n-------------------------------------------\n");  
}

void print_map_diferencias(string mensaje){
    cout<<endl;
    cout<<"DIFERENCIAS "<<mensaje<<"\n\n";
    for(auto it2 = map_diferencias.cbegin(); it2 != map_diferencias.cend(); ++it2){
        printf (ANSI_COLOR_RED "%5.2s ", it2->first.c_str()); 
        for(int i=0;i<it2->second.size();i++){
          printf ("\033[0m%5.2f ", it2->second[i]); 
        }
        cout<<endl;
      }
    cout<<endl;
    cout<<"-------------------------------------------\n";
}


void _build_data_to_distance_matrix(){
    vector< vector<float> > matrix_numeros(_genes_iniciales.size()-1,vector<float>(_genes_iniciales.size()-1));
    for (int i=0;i<matrix_numeros.size();i++){ // to fill differents values from zero into diagonal
        for (int j=0;j<matrix_numeros.size();j++){    if(j==i)    matrix_numeros[i][j]=0;} }

    for(int i=1; i< _genes_iniciales.size()-1; i++){
        for(int j=0; j<i; j++){
            for(int k=0; k<numeros[i].size();k++){
                if(j!=i)    acumulador+=((numeros[j][k]-numeros[i][k])*(numeros[j][k]-numeros[i][k])); 	    }
            if (j!=i){
                float _calculus=sqrt(acumulador);   acumulador=0;
                matrix_numeros[i][j]=matrix_numeros[j][i]=_calculus;
            }
        }
    }
    matrix_numeros_final=matrix_numeros;
    // print_matrix_distances("DISTANCIAS");
}

void _get_minimum(){
    for(int i=0;i<matrix_numeros_final.size();i++){
        float menor=40;
        for  (int j=0;j<matrix_numeros_final.size();j++){
            if(i!=j){
                if(matrix_numeros_final[i][j]<menor)    menor=matrix_numeros_final[i][j];
            }}
        _vector_fill.push_back(menor);}
}
void _get_maximum(){
    for(int i=0;i<matrix_numeros_final.size();i++){
        float mayor=-40;
        for  (int j=0;j<matrix_numeros_final.size();j++){
            if(i!=j){
                if(matrix_numeros_final[i][j]>mayor)    mayor=matrix_numeros_final[i][j];}}
        _vector_fill.push_back(mayor);}
}



void _from_get_minimun_to_diferencias(){
    if(tipo=="minimo")  _get_minimum();
    if(tipo=="maximo")  _get_maximum();

  for(int i=0;i<_vector_fill.size();i++){
    vector <float> aux;
    aux.push_back(_vector_fill[i]); aux.push_back(0);   aux.push_back(0);
    map_diferencias[headers[i]] = (aux);
  }
}

string _get_element_maximum_from_diferencias_vect(){
  float _inicio=-19;
  string gen;
  for(auto it2 = map_diferencias.cbegin(); it2 != map_diferencias.cend(); ++it2){
    if(it2->second[0]>=_inicio){
      _inicio=it2->second[0];       gen = it2->first;
    }
  }
  return gen;
}

void _armado_de_matrix_primer_paso(){
   string gen = _get_element_maximum_from_diferencias_vect();
    map_diferencias.erase (gen);
    int pos_mayor = find(headers.begin(), headers.end(), gen) - headers.begin();

    vector<string> vect_string;
    vect_string.push_back(headers[pos_mayor]);

    // print_map_diferencias("INICIO, NO GENES, inicial");
    for(auto it3 = map_diferencias.cbegin(); it3 != map_diferencias.cend(); ++it3){
      string gen_uno = it3->first;
      int pos_uno = find(headers.begin(), headers.end(), gen_uno) - headers.begin();
      map_diferencias[gen_uno][1]=matrix_numeros_final[pos_uno][pos_mayor];
      map_diferencias[gen_uno][2]=map_diferencias[gen_uno][0]-map_diferencias[gen_uno][1];
    }
    // print_map_diferencias("INICIO, NO GENES, final");

    vector<string> cluster_chico;
    int flag=0;
    for(auto it3 = map_diferencias.cbegin(); it3 != map_diferencias.cend(); ++it3){
      if(it3->second[2]>=0){
        cluster_chico.push_back(it3->first);
        flag=1;
      }
    }

    if(flag==0){ // quiere decir que ningun elemento se agrega al cluster actual;
      _genes_clusteres.push_back(vect_string);
    }
    if(flag==1){//quiere decir que hay positivos y se agregaran al cluster actual
      for(int i=0;i<cluster_chico.size();i++){
        vect_string.push_back(cluster_chico[i]);
        map_diferencias.erase(cluster_chico[i]);
      }
      _genes_clusteres.push_back(vect_string);
    }
}

void _loop_for_genes_clusteres(){
    string gen = _get_element_maximum_from_diferencias_vect();
    map_diferencias.erase (gen);
    int pos_mayor = find(headers.begin(), headers.end(), gen) - headers.begin();

    // print_map_diferencias("UN GEN, inicial");
    vector<string> vect_string;
    vect_string.push_back(headers[pos_mayor]);
    for(auto it3 = map_diferencias.cbegin(); it3 != map_diferencias.cend(); ++it3){
        string gen_uno = it3->first;
        int pos_uno = find(headers.begin(), headers.end(), gen_uno) - headers.begin();
        map_diferencias[gen_uno][1]=matrix_numeros_final[pos_uno][pos_mayor];
        map_diferencias[gen_uno][2]=it3->second[0]-it3->second[1];
    }
    // print_map_diferencias("UN GEN, final");

    vector<string> cluster_chico;
    int flag=0;
    for(auto it3 = map_diferencias.cbegin(); it3 != map_diferencias.cend(); ++it3){
        if(it3->second[2]>=0){
            cluster_chico.push_back(it3->first);
            flag=1;
        }
    }

    if(flag==0){ // quiere decir que ningun elemento se agrega al cluster actual;
        _genes_clusteres.push_back(vect_string);
    }
    if(flag==1){//quiere decir que hay positivos y se agregaran al cluster actual
        for(int i=0;i<cluster_chico.size();i++){
            vect_string.push_back(cluster_chico[i]);
            map_diferencias.erase(cluster_chico[i]);
        }
        _genes_clusteres.push_back(vect_string);
    }
}

void _loop_genes_clusteres_mas_de_dos(){

    vector<string> copia;
    for(int i=0;i<_genes_clusteres[_genes_clusteres.size()-1].size();i++){
        copia.push_back(_genes_clusteres[_genes_clusteres.size()-1][i]);
    }
    // print_map_diferencias("DOS GENES, inicial");
    for(auto it3 = map_diferencias.cbegin(); it3 != map_diferencias.cend(); ++it3){
        string gen_uno = it3->first;
        int pos_uno = find(headers.begin(), headers.end(), gen_uno) - headers.begin();
        float numero;
        if(tipo=="minimo"){
            float menor=40;
            for(int i=0;i<copia.size();i++){
                int pos_nido = find(headers.begin(), headers.end(), copia[i]) - headers.begin();
                if(matrix_numeros_final[pos_nido][pos_uno]<menor){
                    menor=matrix_numeros_final[pos_nido][pos_uno];
                }
            }
            numero = menor;
        }
        map_diferencias[gen_uno][1]=numero;
        map_diferencias[gen_uno][2]=map_diferencias[gen_uno][0]-numero;
    }
    // print_map_diferencias("DOS GENES, final");
    vector<string> vect_string;
    vector<string> cluster_chico;
    int flag=0;
    for(auto it3 = map_diferencias.cbegin(); it3 != map_diferencias.cend(); ++it3){
        if(it3->second[2]>=0){
            cluster_chico.push_back(it3->first);
            flag=1;
        }
    }
    
    if(flag==1){//quiere decir que hay positivos y se agregaran al cluster actual
        for(int i=0;i<cluster_chico.size();i++){
            vect_string.push_back(cluster_chico[i]);
            map_diferencias.erase(cluster_chico[i]);
        }
        _genes_clusteres.push_back(vect_string);
    }
    if (flag==0){
        _loop_for_genes_clusteres();
    }

}


void _clustering_action(){
    if(_genes_clusteres.size()>0 and _genes_clusteres[_genes_clusteres.size()-1].size()==1){
            _loop_for_genes_clusteres();
            // print_genes_clusteres();
    }
    else if(_genes_clusteres.size()>0 and _genes_clusteres[_genes_clusteres.size()-1].size()>1){
            _loop_genes_clusteres_mas_de_dos();
            // print_genes_clusteres();
    }

    if(_genes_clusteres.size()<1){
        _armado_de_matrix_primer_paso();
        // print_genes_clusteres();
    }
    // print_genes_clusteres();
    _initial = clock(); double _res_time = (double(_initial-_initial_sub))/100000000;  _time_concurrent =_time_concurrent+_res_time;
}



void _clusterization(){
    headers=_genes_iniciales;
    _from_get_minimun_to_diferencias();

    int num_clusters_here=0;
    do{
        _clustering_action();
        num_clusters_here=_genes_clusteres.size()+1;
        // cout<<"clusteres contados "<<num_clusters_here<<endl;
    }
    // while(num_clusters_here<num_clusteres_finales and map_diferencias.size()>2);
    while(num_clusters_here<num_clusteres_finales );

    cout<<"\nClusteres contados "<<num_clusters_here<<endl;
    cout<<"\nTiempo: " << _time_concurrent <<" secs \n\n"<< endl;
    print_genes_clusteres();
    cout<<endl;
    
}



int main(){
    tipo ="minimo";    

    num_clusteres_finales=200;

	_from_txt_to_vect();
    _build_data_to_distance_matrix();
    _clusterization();
	return 0;
}
