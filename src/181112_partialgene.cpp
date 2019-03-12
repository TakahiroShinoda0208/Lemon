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
#include "/data/takahiro/work/Annotation_tool/tool_dev/function.hpp"
using namespace std;

//------------------------------ struct宣言
struct ice
{
      int ma,de,au,sn,ho;
      vector<string> anno;
};

//------------------------------- 宣言
void Output(unordered_map<string,ice>hash,ofstream &fout);
void Complete(ifstream &fin,ofstream &fout);
void Partial(ifstream &fin,ifstream &fin2,ofstream &fout,int cov);


//------------------------------------ main function
int main(int argc,char**argv)
{
      ifstream fin,fin2; //Group gf
      ofstream fout; //output file
      int flag=0,cov=80;
      
      for(int i=0;i<argc;i++){
            string ss=argv[i];
            if(ss=="-i"){fin.open(argv[i+1]);}
            if(ss=="-i2"){fin2.open(argv[i+1]);}
            if(ss=="-o"){fout.open(argv[i+1]);}
            if(ss=="-c"){flag=1;}
            if(ss=="-p"){flag=2;cov=atoi(argv[i+1]);}
      }

//-------------------------------------- 入力ファイルが存在しなかった時の出力

            if(argc<=2){
                  cout << endl;
                  cout <<"\t"<<"Annotation Completetion tool"<<"\n\n\n";
                  cout <<"version 1.0" <<"\n";
                  cout <<"updated 2018/11/09"<<"\t"<<"TAKAHIRO SHINODA"<<"\n\n\n";
                  cout <<"---How to use---"<<"\n";
                  cout <<"[Complete mode]\t"<< "./a.out -i Group.gff -o outputfile -c " <<"\n\n";
                  cout <<"[Partial mode]\t"<< "./a.out -i Group.gff -i2 Result.gff -o outputfile -p 80" <<"\n\n";
                  cout <<"----------------"<<"\n\n\n";
                  return 1;
            }
             if(flag==1){
                   Complete(fin,fout);
             }else if(flag==2){
                   Partial(fin,fin2,fout,cov);
             }
             
}
//------------------------------------------------ Complete
void Complete(ifstream &fin,ofstream &fout){

//declare variable
      vector<string> Vec,lemon;
      string lin,melon,name;
      unordered_map<string,ice>hash;
      ice box;      

//最初の行
      
      getline(fin,lin);
      Vec = split(lin,'\t');
      melon=Vec[3]+"\tcandidate\t"+Vec[5]+"\t"+Vec[6]+"\t"+Vec[7]+"\t0\t"+Vec[9]+"\t"+Vec[10]+"\t"+Vec[0];
      lemon.push_back(melon);
      

//ファイルの処理
      while(getline(fin,lin)){
            Vec = split(lin,'\t');
            if(Vec[5] == "mRNA"){
                  auto itr = hash.find(name);
                  if( itr == hash.end() ) {
                        box = {0,0,0,0,0,lemon};
                  }else{ //すでに存在していれば
                        box = hash[name];
                  }
                  if(Vec[4].find("mapping") !=  string::npos){box.ma++;
                  }else if(Vec[4].find("denovo") !=  string::npos){box.de++;
                  }else if(Vec[4].find("AUGUSTUS") !=  string::npos){box.au++;
                  }else if(Vec[4].find("SNAP") !=  string::npos){box.sn++;
                  }else if(Vec[4].find("homology") !=  string::npos){box.ho++;
                  }
                  hash[name]=box;
                  name="";lemon.clear();
                  melon=Vec[3]+"\tcandidate\t"+Vec[5]+"\t"+Vec[6]+"\t"+Vec[7]+"\t0\t"+Vec[9]+"\t"+Vec[10]+"\t"+Vec[0]; 
                  lemon.push_back(melon);
            }else if(Vec[5] == "CDS"){
                  name += Vec[3]+Vec[6]+Vec[7]+Vec[9];
                  melon=Vec[3]+"\tcandidate\t"+Vec[5]+"\t"+Vec[6]+"\t"+Vec[7]+"\t0\t"+Vec[9]+"\t"+Vec[10]+"\t"+Vec[0]; 
                  lemon.push_back(melon);
            }
      }
//最終行
      auto itr = hash.find(name);
      if( itr == hash.end() ) {
            box = {0,0,0,0,0,lemon};
      }else{ //すでに存在していれば
            box = hash[name];
      }
      if(Vec[4].find("mapping") !=  string::npos){box.ma++;
      }else if(Vec[4].find("denovo") !=  string::npos){box.de++;
      }else if(Vec[4].find("AUGUSTUS") !=  string::npos){box.au++;
      }else if(Vec[4].find("SNAP") !=  string::npos){box.sn++;
      }else if(Vec[4].find("homology") !=  string::npos){box.ho++;
      }
      hash[name]=box;
      Output(hash,fout);
}

//----------------------------------------------------------------------------------------------------------------------------------
void Output(unordered_map<string,ice>hash,ofstream &fout)
{
      vector<string> gff,gfftmp;
      for(auto itr = hash.begin(); itr != hash.end(); ++itr) {
            int count=0,num=0;
            vector<string>Vec;
            
            if(itr->second.ma >=1 ||itr->second.de >=1){count++;}
            if(itr->second.au >=1 ||itr->second.sn >=1){count++;}
            if(itr->second.ho >=1){count++;}
            
            if(count >=2){
                  num++;
                  string tmp;
                  //mRNAの出力
                  Vec = split(itr->second.anno[0],'\t');
                  tmp =  Vec[0]+"\t"+Vec[1]+"\t"+Vec[2]+"\t"+Vec[3]+"\t"+Vec[4]+"\t"+ItoS(count)+"\t"+Vec[6]+"\t"+Vec[7]+"\t"+"ID=comp_group_num_"+Vec[8]+"_"+ItoS(num)+".mrna1;Parent=group_num_"+Vec[8]+"_"+ItoS(num);
                  gff.push_back(tmp);

                  //CDSの出力
                  for(int i=1;i<itr->second.anno.size();i++){
                        Vec = split(itr->second.anno[i],'\t');
                        tmp = Vec[0]+"\t"+Vec[1]+"\t"+Vec[2]+"\t"+Vec[3]+"\t"+Vec[4]+"\t"+ItoS(count)+"\t"+Vec[6]+"\t"+Vec[7]+"\tID=comp_group_num_"+Vec[8]+"_"+ItoS(num)+".mrna1.cds"+ItoS(i)+";Parent=group_num_"+Vec[8]+"_"+ItoS(num)+".mrna1";
                        gff.push_back(tmp);
                  }
            }
      
            // if(count >=2){
            //       num++;
            //       string tmp;
            //       //gene・mRNAの出力
            //       Vec = split(itr->second.anno[0],'\t');
            //       fout << Vec[0]<<"\t"<<Vec[1]<<"\tgene\t"<<Vec[3]<<"\t"<<Vec[4]<<"\t"<<count<<"\t"<<Vec[6]<<"\t"<<Vec[7]<<"\t"<<"ID=comp_group_num_"<<Vec[8]<<"_"<<num<<"\n";
            //       fout << Vec[0]<<"\t"<<Vec[1]<<"\t"<<Vec[2]<<"\t"<<Vec[3]<<"\t"<<Vec[4]<<"\t"<<count<<"\t"<<Vec[6]<<"\t"<<Vec[7]<<"\t"<<"ID=comp_group_num_"<<Vec[8]<<"_"<<num<<".mrna1;Parent="<<"group_num_"<<Vec[8]<<"_"<<num<<"\n";

            //       //exonの出力
            //       for(int i=1;i<itr->second.anno.size();i++){
            //             Vec = split(itr->second.anno[i],'\t');
            //             fout << Vec[0]<<"\t"<<Vec[1]<<"\texon\t"<<Vec[3]<<"\t"<<Vec[4]<<"\t"<<count<<"\t"<<Vec[6]<<"\t"<<Vec[7]<<"\t"<<"ID=comp_group_num_"<<Vec[8]<<"_"<<num<<".mrna1.exon"<<i<<";Parent=group_num_"<<Vec[8]<<"_"<<num<<".mrna1\n";
            //       }
            //       //CDSの出力
            //       for(int i=1;i<itr->second.anno.size();i++){
            //             Vec = split(itr->second.anno[i],'\t');
            //             fout << Vec[0]<<"\t"<<Vec[1]<<"\t"<<Vec[2]<<"\t"<<Vec[3]<<"\t"<<Vec[4]<<"\t"<<count<<"\t"<<Vec[6]<<"\t"<<Vec[7]<<"\t"<<"ID=comp_group_num_"<<Vec[8]<<"_"<<num<<".mrna1.cds"<<i<<";Parent=group_num_"<<Vec[8]<<"_"<<num<<".mrna1\n";
            //       }
            // }
      }
      gfftmp=Grouping(gff);
      longest(gfftmp,fout);
}

//-----------------------------------------
void Partial(ifstream &fin,ifstream &fin2,ofstream &fout,int cov){
      
      unordered_map<string,int>hash2,hash3,lemonade;//hash2はmRNA確認用
      unordered_map<string,ice>hash;//hashはCDS格納用
      vector<string> Vec,emp,cds,gff,gfftmp;
      vector<pair<string,vector<string> > > lemon; //代表配列の格納
      string name,mname,line; //mname , name
      ice box;

//-------------------------------Result.gff
      cout <<"Start Input result.gff\n";
      
//----1行目
      getline(fin2,line);
//----本体      
      while(getline(fin2,line)){
            Vec = split(line,'\t');
            if(Vec[5-3] == "mRNA"){
                  hash3[mname]=1;
                  mname="";
            }else if(Vec[5-3] == "CDS"){
                  mname += Vec[3-3]+Vec[6-3]+Vec[7-3]+Vec[9-3]+Vec[10-3];
            }
      }
//----最終行
      hash3[mname]=1;
      mname="";
      cout <<"Finish Input result.gff\n";
//-------------------------------Group.gff
      cout <<"Start Input Group.gff\n";
//----1行目
      getline(fin,line);
      cds.push_back(line);

//CDSとintron情報を格納するphase
      while(getline(fin,line)){
            Vec = split(line,'\t');
            if(Vec[5] == "mRNA"){
                  auto itr = hash2.find(mname);
                  auto itr2 = hash3.find(mname);
                  if( itr == hash2.end() && itr2 == hash3.end() ) {
                        lemon.push_back(make_pair(mname,cds));
                  }
                  hash[name]=box;
                  hash2[mname]=1;

                  mname="";cds.clear();
                  cds.push_back(line);

            }else if(Vec[5] == "CDS"){

                  name = Vec[3]+Vec[6]+Vec[7]+Vec[9]+Vec[10]; //hash用
                  mname += Vec[3]+Vec[6]+Vec[7]+Vec[9]+Vec[10]; //hash2用
                  cds.push_back(line); //hash2のCDS用
                  
                  auto itr = hash.find(name);
                  if( itr == hash.end() ) {
                        box = {0,0,0,0,0,emp};
                  }else{ //すでに存在していれば
                        box = hash[name];
                  }
                  if(Vec[4].find("mapping") !=  string::npos){box.ma=1;
                  }else if(Vec[4].find("denovo") !=  string::npos){box.de=1;
                  }else if(Vec[4].find("AUGUSTUS") !=  string::npos){box.au=1;
                  }else if(Vec[4].find("SNAP") !=  string::npos){box.sn=1;
                  }else if(Vec[4].find("homology") !=  string::npos){box.ho=1;
                  }
                  hash[name]=box;
                  name="";
            }
      }
//----最終行
      auto itr = hash2.find(mname);
      if( itr == hash2.end() ) {lemon.push_back(make_pair(mname,cds));}
      hash[name]=box;
      mname="";cds.clear();
      cds.push_back(line);
      cout <<"Finish Input Group.gff\n";

//-------------------------------hashの中身整理
      cout <<"Start Calculate hash\n";
//hashの中身を計算する(2つ以上のevidenceによって支持されている//2tool以上に支持されている?)
      for(auto itr = hash.begin(); itr != hash.end(); ++itr) {
            int count=0,num=0;
            vector<string>Vec;
            count = (itr->second.ma + itr->second.de + itr->second.au + itr->second.sn + itr->second.ho);
            
            // if(itr->second.ma >=1 ||itr->second.de >=1){count++;}
            // if(itr->second.au >=1 ||itr->second.sn >=1){count++;}
            // if(itr->second.ho >=1){count++;}      
            lemonade[itr->first]=count;
      }
      cout <<"FInish Calculate hash\n";
      cout <<"Start Calculate lemon\n";
//lemonの代表配列を一つずつ検討していく
      for(int i=0;i<lemon.size();i++){ //mrnaの行以外を処理する
//CDSあ２つ以上のevidenceによって確認されているかを確認
            int lem1=0,lem2=0;
            for(int j=1;j<lemon[i].second.size();j++){ //mrnaの行以外を処理する
                  Vec = split(lemon[i].second[j],'\t');
                  name = Vec[3]+Vec[6]+Vec[7]+Vec[9]+Vec[10];
                  lem1+= stoi(Vec[7])-stoi(Vec[6])+1;
                  if(lemonade[name] >= 2){lem2+= stoi(Vec[7])-stoi(Vec[6])+1;}
            }
            //条件を満たすCDSlengthが全体の8割を超えた時出力
            if((double)lem2/(double)lem1*100 >= cov ){ //変数に指定できる
                  gff.insert(gff.end(), lemon[i].second.begin(), lemon[i].second.end());
            }
      }
      cout <<"Finish Calculate lemon\n";
//被っている配列からlongestを出力
      longest(gff,fout);
}
