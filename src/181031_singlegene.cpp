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

///宣言
vector<string> split(const string &str, char sep);
void Calculate(vector< pair<int , string > > &CDSrep,vector<string> &CDSlist,int gro_num,ofstream &fout,int cov); 
void Out(int gro_num, pair<int , string > &Single,int score,ofstream &fout);


//------------------------------------ main function
int main(int argc,char**argv)
{
      ifstream fin; //gff result
      ofstream fout; //output file
      int cov=70;
      
      for(int i=0;i<argc;i++){
            string ss=argv[i];
            if(ss=="-gff"){fin.open(argv[i+1]);}
            if(ss=="-o"){fout.open(argv[i+1]);}
            if(ss=="-cov"){cov=atoi(argv[i+1]);}
      }

//入力ファイルが存在しなかった時の出力
      if(argc<=2){            
            cout << endl;
            cout <<"\t"<<"Single exon predict tool"<<"\n\n\n";
            cout <<"version 1.0" <<"\n";
            cout <<"updated 2018/10/31"<<"\t"<<"TAKAHIRO SHINODA"<<"\n\n\n";
            cout <<"Single exon predict tool\n";
            cout <<"---How to use---"<<"\n";
            cout << "[Search all combination]"<<"\t"<< "./a.out -gff rna.gff -o outputfile -cov 70" <<"\n\n";
            cout <<"----------------"<<"\n\n\n";
            return 1;
      }
      
      //Groupごとに計算を行う。
      cout << "Start Single exon CDS prediction\n";

      //declare variable
      string lin;
      vector<string> Vec,CDSlist;
      vector< pair<int , string > > CDSrep;
      int gro_num=0,flag=100;

//ファイル一行目の処理
      getline(fin,lin);
      Vec = split(lin,'\t');
      gro_num=stoi(Vec[0]);

//ファイルの処理
      while(getline(fin,lin)){
            Vec = split(lin,'\t');
            //Group numが異なる場合、ST/CDS/EDのfilteringを実行
            if(gro_num!=stoi(Vec[0])){
                  //calculate
                  sort(CDSrep.begin(),CDSrep.end());
                  CDSrep.erase(unique(CDSrep.begin(), CDSrep.end()), CDSrep.end());
                  //cout << CDSrep.size() << "\t"<<CDSlist.size()<<"\n";
                  Calculate(CDSrep,CDSlist,gro_num,fout,cov);
                  //initialize
                  gro_num=stoi(Vec[0]);
                  CDSrep.clear();
                  CDSlist.clear();
            }else{
                  if(Vec[5] =="CDS"){
                        CDSlist.push_back(lin); //後ほど処理する
                        int st,ed=0;
                        if(stoi(Vec[6])>stoi(Vec[7])){st=stoi(Vec[7]);ed=stoi(Vec[6]);
                        }else{st=stoi(Vec[6]);ed=stoi(Vec[7]);}
                        CDSrep.push_back(make_pair(st+ed,Vec[3]+"\t"+Vec[6]+"\t"+Vec[7]+"\t"+Vec[9]));
                  }
            }
      }
//ファイル最終行の処理
      sort(CDSrep.begin(),CDSrep.end());
      CDSrep.erase(unique(CDSrep.begin(), CDSrep.end()), CDSrep.end());
      Calculate(CDSrep,CDSlist,gro_num,fout,cov);
      
      cout << "Finish Single exon CDS prediction\n";
      return 0;
}

//---------------------------------------------------Split
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


//---------------------------------------------------Calculate Single exon gene
void Calculate(vector< pair<int , string > > &CDSrep,vector<string> &CDSlist,int gro_num,ofstream &fout,int cov){

//delcare variable
      int dp[CDSrep.size()+1][7]; //len RNAseq Abinitio Homology
      vector<string> TMP;

//initialize table      
      for(int i=0;i<CDSrep.size();i++){
            TMP = split(CDSrep[i].second,'\t');
            int st,ed,len=0;
            if(stoi(TMP[1])>stoi(TMP[2])){st=stoi(TMP[2]);ed=stoi(TMP[1]);
            }else{st=stoi(TMP[1]);ed=stoi(TMP[2]);}
            dp[i][0] = ed-st+1;
            dp[i][1] = st;
            dp[i][2] = ed;
            dp[i][3] = 0;
            dp[i][4] = 0;
            dp[i][5] = 0;
            dp[i][6] = 0;
      }

//make table
      for(int i=0;i<CDSlist.size();i++){
            TMP = split(CDSlist[i],'\t');
            int st,ed,len=0;
            if(stoi(TMP[6])>stoi(TMP[7])){st=stoi(TMP[7]);ed=stoi(TMP[6]);
            }else{st=stoi(TMP[6]);ed=stoi(TMP[7]);}
            len=ed-st+1;
            
            for(int j=0;j<CDSrep.size();j++){
                  //filtering条件
                  if(stoi(TMP[6]) >= dp[j][1] && stoi(TMP[7]) <= dp[j][2]){
                        if((dp[j][1]-stoi(TMP[6]))%3 == 0 && (dp[j][2]-stoi(TMP[7]))%3 == 0){
                              if((double)len/(double)dp[j][0]*100 >= (double)cov){
                                    //RNAseq・Abinitio・Homologyに分類
                                    if(TMP[4].find("mapping") !=  string::npos ||TMP[4].find("denovo") !=  string::npos){
                                          dp[j][3] ++;
                                    }else if(TMP[4].find("AUGUSTUS") !=  string::npos || TMP[4].find("SNAP") !=  string::npos){
                                          dp[j][4] ++;
                                    }else if(TMP[4].find("homology") !=  string::npos){
                                          dp[j][5] ++;
                                    }else{
                                          
                                          cout << "認識できない名前がついています"<<TMP[4]<<"\n";

                                          return;
                                    }
                              }
                        }
                  }
            }
      }

//declare
      vector < pair < int , pair < int,string > > > Can;
      vector < pair < int , pair < int,string > > > Last;
      
//Judge whether output or not
      for(int i=0;i<CDSrep.size();i++){
            int sum=0;
            for(int j=3;j<=5;j++){
                  if(dp[i][j]>=1){sum++;}
            }

            //Scoring
            if(sum==3){
                  CDSrep[i].first=dp[i][0];
                  Can.push_back(make_pair(sum,CDSrep[i]));
                  //}else if(sum==2 && dp[i][0] >= 300){
            }else if(sum==2){
                  CDSrep[i].first=dp[i][0];
                  Can.push_back(make_pair(sum,CDSrep[i]));
                  //}else if(sum==1 && dp[i][0] >= 100){
            }else{}
            // if(sum==1){
            //       CDSrep[i].first=dp[i][0];
            //       Can.push_back(make_pair(sum,CDSrep[i]));
            //}
      }
      //ここ要確認
      sort(Can.begin(),Can.end());
      // for(int i=0;i<Can.size();i++){
      //       cout << Can[i].first <<"\t"<< Can[i].second.first<<"\t"<< Can[i].second.second<<"\n";
      // }
            
      if(Can.size()!=0){            
            Last.push_back(Can[Can.size()-1]);
            for(int j=Can.size()-2;j>=0;j--){
                  for(int k=0;k<Last.size();k++){
                        if(Last[k]==Can[j]){
                              vector<string>tmp1,tmp2;
                              tmp1 = split(Last[k].second.second,'\t');
                              tmp2 = split(Can[k].second.second,'\t');
                              if((stoi(tmp1[1])<=stoi(tmp2[1]) && stoi(tmp2[1])<=stoi(tmp1[2]))||(stoi(tmp1[1])<=stoi(tmp2[2]) && stoi(tmp2[2])<=stoi(tmp1[2]))){
                              }else{
                                    Last.push_back(Can[0]); 
                              }
                        }
                  }
            }
//Outで出力
            for(int i=0;i<Last.size();i++){
                  Out(gro_num,Last[i].second,Last[i].first,fout);
            }
      }
}
//----------------------------------------------------------------------------------------------------------------------------------

void Out(int gro_num, pair<int , string > &Single,int score,ofstream &fout)
{
      int num=1;
      vector<string> TMP;
      TMP = split(Single.second,'\t');
      fout <<TMP[0]<<"\tcandidate\tgene\t"<<TMP[1]<<"\t"<<TMP[2]<<"\t"<<score<<"\t"<<TMP[3]<<"\t."<<"\t"<<"ID=single_group_num_"<<gro_num<<"_"<<num<<"\n";
      fout <<TMP[0]<<"\tcandidate\tmRNA\t"<<TMP[1]<<"\t"<<TMP[2]<<"\t"<<score<<"\t"<<TMP[3]<<"\t."<<"\t"<<"ID=single_group_num_"<<gro_num<<"_"<<num<<".mrna1;Parent="<<"single_group_num_"<<gro_num<<"_"<<num<<"\n";
      fout <<TMP[0]<<"\tcandidate\texon\t"<<TMP[1]<<"\t"<<TMP[2]<<"\t"<<score<<"\t"<<TMP[3]<<"\t0"<<"\t"<<"ID=single_group_num_"<<gro_num<<"_"<<num<<".mrna1.exon"<<num<<";Parent=single_group_num_"<<gro_num<<"_"<<num<<".mrna1\n";
      fout <<TMP[0]<<"\tcandidate\tCDS\t"<<TMP[1]<<"\t"<<TMP[2]<<"\t"<<score<<"\t"<<TMP[3]<<"\t0"<<"\t"<<"ID=single_group_num_"<<gro_num<<"_"<<num<<".mrna1.cds"<<num<<";Parent=single_group_num_"<<gro_num<<"_"<<num<<".mrna1\n";

      // cout <<TMP[0]<<"\tcandidate\tgene\t"<<TMP[1]<<"\t"<<TMP[2]<<"\t"<<score<<"\t"<<TMP[3]<<"\t."<<"\t"<<"ID=group_num_"<<gro_num<<"_"<<num<<"\n";
      // cout <<TMP[0]<<"\tcandidate\tmRNA\t"<<TMP[1]<<"\t"<<TMP[2]<<"\t"<<score<<"\t"<<TMP[3]<<"\t."<<"\t"<<"ID=group_num_"<<gro_num<<"_"<<num<<".mrna1;Parent="<<"group_num_"<<gro_num<<"_"<<num<<"\n";
      // cout <<TMP[0]<<"\tcandidate\texon\t"<<TMP[1]<<"\t"<<TMP[2]<<"\t"<<score<<"\t"<<TMP[3]<<"\t0"<<"\t"<<"ID=group_num_"<<gro_num<<"_"<<num<<".mrna1.exon"<<num<<";Parent=group_num_"<<gro_num<<"_"<<num<<".mrna1\n";
      // cout <<TMP[0]<<"\tcandidate\tCDS\t"<<TMP[1]<<"\t"<<TMP[2]<<"\t"<<score<<"\t"<<TMP[3]<<"\t0"<<"\t"<<"ID=group_num_"<<gro_num<<"_"<<num<<".mrna1.cds"<<num<<";Parent=group_num_"<<gro_num<<"_"<<num<<".mrna1\n";
      
}
