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
using namespace std;


//struct宣言
struct group
{      
      int st,ed;
      string chr,tool;
      vector<string> gff;
      
      bool operator<(const group &another) const
      {     
            return chr == another.chr ? st < another.st : chr < another.chr;
      };
};

//関数プロトタイプ宣言
vector<string> split(const string &str, char sep);
string ItoS(int number);
void Grouping(ifstream &ifs);
void Output(ofstream &fout,vector<group> &Box,int &gro_num);


//main関数
int main(int argc,char**argv)
{
      ifstream fin; //mapping result
      int num = 1;
      int flag=0;
      for(int i=0;i<argc;i++){
            string ss=argv[i];
            if(ss=="-f"){fin.open(argv[i+1]);}
      }

//入力ファイルが存在しなかった時の出力
      if(argc<=2){            
            cout << endl;
            cout <<"\t"<<"annotation_grouping tool"<<"\n\n\n";
            cout <<"version 1.0" <<"\n";
            cout <<"updated 2018/10/04"<<"\t"<<"TAKAHIRO SHINODA"<<"\n\n\n";

            cout <<"---How to use---"<<"\n";
            cout << "[Group mode]"<<"\t\t\t"<< "0. Group -f rna.gff" <<"\n";
            cout <<"------------"<<"\n\n\n";
            return 1;
      }
      
//-----------------------------------------------------------------------------------------------------
      if(flag==0){
            cout << "START Grouping." << endl;
            Grouping(fin);
            cout << "END Grouping." << endl;
      }
      return 0;
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

void Grouping(ifstream &ifs)
{
       
//fileを行単位で読みこみ、chr・stでsortする
      ofstream fout("Group.gff");
      string lin;
      group gff;
      vector<string> I,gff_vec;
      vector<group> Box;
      vector<group> Box2;
      int st=0,ed=0;
      string chr,tool;
      string strand;
      int gro_num=1;

//一行目
      getline(ifs,lin);
      gff_vec.push_back(lin);
      I = split(lin,'\t');
      tool = I[1];st = stoi(I[3]);ed = stoi(I[4]);
      chr = I[0];strand = I[6];

      while(getline(ifs,lin)){
            I = split(lin,'\t');
            if(I[2] == "mRNA" || I[2] == "gene")
            {
                  gff = {st,ed,chr,tool,gff_vec};
                  if(strand=="+"){Box.push_back(gff);
                  }else{Box2.push_back(gff);}                  
                  tool = I[1];
                  st = stoi(I[3]);ed = stoi(I[4]);
                  chr = I[0];strand = I[6];
                  gff_vec.clear();gff_vec.push_back(lin);
            }    
            else if(I[2] == "CDS"){gff_vec.push_back(lin);}
      }
      gff = {st,ed,chr,tool,gff_vec};
      if(strand=="+"){Box.push_back(gff);
      }else{Box2.push_back(gff);}                  
      sort(Box.begin(),Box.end());
      sort(Box2.begin(),Box2.end());
      Output(fout,Box,gro_num);
      Output(fout,Box2,gro_num);
}

void Output(ofstream &fout,vector<group> &Box,int &gro_num)
{
//Groupingして、mapping_baseの数を確認する(0,1,2>=)

      int ele_num=0,map_num=0,tmpst=Box[0].st,tmped=Box[0].ed;
      string tmpchr=Box[0].chr;

      for(int i=0;i<Box.size();i++){
            if(tmpchr == Box[i].chr)
            {
                  if(Box[i].st <= tmped ){
                        if(tmped < Box[i].ed ){tmped = Box[i].ed;}           
                        ele_num ++;
                        if(Box[i].tool == "mappingbase"){map_num++;}  
                  }else{
                              cout <<gro_num<<"\t"<<ele_num<<"\t"<<map_num<<"\n";
                              
                              for(int j=i-ele_num;j<i;j++){
                                    for(int k=0;k<Box[j].gff.size();k++){
                                          fout <<gro_num<<"\t"<<ele_num<<"\t"<<map_num<<"\t"<<Box[j].gff[k]<<"\n";
                                    }
                              }
                        //初期化
                        gro_num++;ele_num=1;map_num=0;
                        tmpst = Box[i].st;
                        tmped = Box[i].ed;  
                        if(Box[i].tool == "mappingbase"){map_num ++;}
                  }
            }
            else{
                        cout <<gro_num<<"\t"<<ele_num<<"\t"<<map_num<<"\n";

                        for(int j=i-ele_num;j<i;j++){
                              for(int k=0;k<Box[j].gff.size();k++){
                                    fout <<gro_num<<"\t"<<ele_num<<"\t"<<map_num <<"\t"<<Box[j].gff[k]<<"\n";
                              }
                        }
                  //初期化
                  gro_num++;
                  ele_num=1;
                  map_num=0;
                  tmpchr = Box[i].chr;
                  tmpst = Box[i].st;
                  tmped = Box[i].ed;  
                  if(Box[i].tool == "mappingbase"){map_num ++;}
            }            
      }
            cout <<gro_num<<"\t"<<ele_num<<"\t"<<map_num<<"\n";
            for(int j=Box.size()-ele_num;j<Box.size();j++){
                  for(int k=0;k<Box[j].gff.size();k++){
                        fout <<gro_num<<"\t"<<ele_num<<"\t"<<map_num <<"\t"<<Box[j].gff[k]<<"\n";
                  }
            }
}

