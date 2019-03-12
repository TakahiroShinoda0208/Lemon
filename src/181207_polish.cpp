#include <iostream>//標準データ入出力
#include <fstream>//ファイルの入出力
#include <vector>//vector(動的配列クラス)
#include <string>//string(文字列クラス)
#include <sstream>//文字列の入出力
#include <unordered_map>//unordered_map関数の導入
#include <map>//map関数の導入
#include <set>//set関数の導入
#include <stdlib.h>//絶対値を用いるための関数
#include <algorithm>
#include <iomanip>//少数点以下表示
#include <time.h>
#include <bitset>
#include <boost/dynamic_bitset.hpp>
using namespace std;

//関数プロトタイプ宣言
vector<string> split(const string &str, char sep);
string ItoS(int number);

int main(int argc,char**argv)
{
      ifstream fin,fin2; //fin1=gff result,fin2=weight.txt
      ofstream fout("phase5.gff");
      string line;
      for(int i=0;i<argc;i++){
            string ss=argv[i];
            if(ss=="-f"){fin.open(argv[i+1]);}
      }

//入力ファイルが存在しなかった時の出力
      if(argc<=2){            
            cout << endl;
            cout <<"\t"<<"gff polish tool"<<"\n\n\n";
            cout <<"version 1.0" <<"\n";
            cout <<"updated 2018/12/07"<<"\t"<<"TAKAHIRO SHINODA"<<"\n\n\n";

            cout <<"---How to use---"<<"\n";
            cout << "[Polish mode]"<<"\t\t\t"<< "0. Polish -f rna.gff" <<"\n";
            cout <<"------------"<<"\n\n\n";
            return 1;
      }
      while(getline(fin,line)){
            vector<string> tmp;
            tmp=split(line,'\t');
            if(tmp[7]=="1"){tmp[7]="2";}
            else if(tmp[7]=="2"){tmp[7]="1";}
            fout << tmp[0]<<"\t"<<tmp[1]<<"\t"<<tmp[2]<<"\t"<<tmp[3]<<"\t"<<tmp[4]<<"\t"<<tmp[5]<<"\t"<<tmp[6]<<"\t"<<tmp[7]<<"\t"<<tmp[8]<<"\n";
      }
      

}

//split関数
vector<string> split(const string &str, char sep)
{
      vector<string> v;
      stringstream ss(str);
      string buffer;
      while( getline(ss, buffer, sep) ) {
            v.push_back(buffer);     
      }
      return v;      
}

//########intをstringに変える関数
string ItoS(int number)
{
      stringstream ss;
      ss << number;      
      return ss.str();     
}

