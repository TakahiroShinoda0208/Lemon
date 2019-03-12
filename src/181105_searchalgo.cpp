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
#include "/data/takahiro/work/Annotation_tool/tool_dev/evmlike9.hpp" //得点
//#include "/data/tkazuki/annotation/annotationpipeline/score_function/evmlike10.hpp"
#include "/data/takahiro/work/Annotation_tool/tool_dev/function.hpp" //function関数
using namespace std;

//------------------------------------ declare function prototype
void CDS_selection(ifstream &ifs,ofstream &fout,unordered_map<string,string>hash,unordered_map<string,double> hashcds,unordered_map<string,double> hashintron,unordered_multimap<string,groupinter> hashinter,unordered_multimap<string,string> hashinter2,unordered_map<string, double> hashweight,int out);
void Combination(vector<pair<int , string > > CDS_list, int gro_num,ofstream &fout,unordered_map<string,string>hash,unordered_map<string,double> hashcds,unordered_map<string,double> hashintron,unordered_multimap<string,groupinter> hashinter,unordered_multimap<string,string> hashinter2,unordered_map<string, double> hashweight,int out);
void OUTPUT(int p,vector<vector<int> > inter_list,vector<vector<string> > correct_list,ofstream &fout,int gro_num,int &num);
void DP(vector<pair<int , string > > CDS_list, int gro_num,int &num,ofstream &fout,unordered_map<string,string>hash,unordered_map<string,double> hashcds,unordered_map<string,double> hashintron,unordered_multimap<string,groupinter> hashinter,unordered_multimap<string,string> hashinter2,unordered_map<string, double> hashweight ,string w,vector<long long int> &score_list,vector< vector<string> > &correct_list,vector< vector<int> > &inter_list);
void trace_back(vector<pair<int , string > > CDS_list, unordered_map <int,pair<long long int ,int> > maxdp,int t,vector <string> &correct,vector <int> &inter,int *intdp);


//------------------------------------ main function
int main(int argc,char**argv)
{
      ifstream fin; //gff result
      ifstream fin2; //genome fasta
      ifstream fin3; //weight file
      ofstream fout; //output file
      int out=0;

      for(int i=0;i<argc;i++){
            string ss=argv[i];
            if(ss=="-fa"){fin2.open(argv[i+1]);}
            if(ss=="-gff"){fin.open(argv[i+1]);}
            if(ss=="-w"){fin3.open(argv[i+1]);}
            if(ss=="-o"){fout.open(argv[i+1]);}
            if(ss=="-a"){out=1;}
            //if(ss=="-o"){fout.open(argv[i+1]);}
      }

//入力ファイルが存在しなかった時の出力
      if(argc<=2){            
            cout << endl;
            cout <<"\t"<<"Search algorithm"<<"\n\n\n";
            cout <<"version 1.5" <<"\n";
            cout <<"updated 2018/12/03"<<"\t"<<"TAKAHIRO SHINODA"<<"\n\n\n";
            cout <<"Scoringのバグを修正しました\n";
            cout <<"---How to use---"<<"\n";
            cout << "[Search all combination]"<<"\t"<< "./a.out -fa genome.fa -gff rna.gff -w weight.txt -o outputfile" <<"\n\n";
            cout << "[Option]"<<"\t"<< "-a :output both strand \t[default one strand]" <<"\n\n";
            cout <<"----------------"<<"\n\n\n";
            return 1;
      }
      
      //declare variable
      unordered_map<string,string>hash;     
      unordered_map<string, double> hashweight;
      unordered_map<string,double>hashcds,hashintron;
      unordered_multimap<string, groupinter> hashinter;
      unordered_multimap<string, string> hashinter2;
      

      cout << "Start making hash table\n";
      Genome(fin2,hash);
      cout << "End making hash table\n";
      cout << "Start CDS scoring\n";
      hashweight = Weight_hash(fin3);
      hashcds = Score_cds_hash(fin,hashweight);
      hashintron = Score_intron_hash(fin,hashweight);
      hashinter = Score_intergenic_hash(fin);
      hashinter2 = Score_intergenic_hash2(fin);
      cout << "End CDS scoring\n";
      cout << "Start CDS clustering\n";
      CDS_selection(fin,fout,hash,hashcds,hashintron,hashinter,hashinter2,hashweight,out);
      cout << "End CDS clustering\n";
      return 0;
}

//---------------------------------------------------filtering CDS which be used for search program
void CDS_selection(ifstream &ifs,ofstream &fout,unordered_map<string,string>hash,unordered_map<string,double> hashcds,unordered_map<string,double> hashintron,unordered_multimap<string,groupinter> hashinter,unordered_multimap<string,string> hashinter2,unordered_map<string, double> hashweight,int out)
{
//delcare variable
      string lin;
      vector<string> Vec;
      int gro_num=0,flag=100;
      pair<int , string> tmp;
      vector<pair<int , string > > CDS_list;
      vector< pair<int , string > > CDStmp;

//ファイル一行目の処理
      ifs.clear();
      ifs.seekg(0, ios_base::beg);
      getline(ifs,lin);
      Vec = split(lin,'\t');
      gro_num=stoi(Vec[0]);
      flag=3;//initialexonのflagに変更

//ファイルの処理
      while(getline(ifs,lin)){
            Vec = split(lin,'\t');
            //Group numが異なる場合、ST/CDS/EDのfilteringを実行
            if(gro_num!=stoi(Vec[0])){
                  if(Vec[5] =="mRNA"||Vec[5] =="gene"){
                        CDS_set(CDS_list,CDStmp,hash); //前行のCDSをCDS_list / ED_listに格納
                        CDS_filter(CDS_list); //CDS filtering
                        Combination(CDS_list,gro_num,fout,hash,hashcds,hashintron,hashinter,hashinter2,hashweight,out);
                        //削除
                        CDS_list.clear();
                        gro_num=stoi(Vec[0]);
                  }else{
                        cout << "想定外の文字が含まれています。\n";
                        return;
                  }
            }else{
                  if(Vec[5] =="CDS"){                  
                        pair<int,string> tmp;
                        
                        if(Vec[9] == "+"){tmp = make_pair(stoi(Vec[6])+stoi(Vec[7])+stoi(Vec[10]),Vec[3]+"\t"+Vec[6]+"\t"+Vec[7]+"\t"+Vec[10]+"\t"+Vec[9]+"\t"+"STbase"+"\t"+"EDbase"+"\t");
                        }else{tmp = make_pair(stoi(Vec[6])+stoi(Vec[7])+stoi(Vec[10])+1,Vec[3]+"\t"+Vec[6]+"\t"+Vec[7]+"\t"+Vec[10]+"\t"+Vec[9]+"\t"+"STbase"+"\t"+"EDbase"+"\t");}
                        CDStmp.push_back(tmp);
                  }else if(Vec[5] =="mRNA"||Vec[5] =="gene"){CDS_set(CDS_list,CDStmp,hash); //前行のCDSをCDS_list / ED_listに格納
                  }
            }
      }

//ファイル最終行の処理
      CDS_set(CDS_list,CDStmp,hash); //前行のCDSをCDS_list / ED_listに格納
      CDS_filter(CDS_list); //CDS filtering
      Combination(CDS_list,gro_num,fout,hash,hashcds,hashintron,hashinter,hashinter2,hashweight,out);
}

void Combination(vector<pair<int , string > > CDS_list,int gro_num,ofstream &fout,unordered_map<string,string>hash,unordered_map<string,double> hashcds,unordered_map<string,double> hashintron,unordered_multimap<string,groupinter> hashinter,unordered_multimap<string,string> hashinter2,unordered_map<string, double> hashweight,int out)
{
      //declare variable
      vector<string> TMP;
      int num=0;
      //decalre tmporary variable
      vector<pair<int , string > > CDS_list_p,CDS_list_m;
      //+ , - で max scoreを記録するほうを出力
      vector<long long int> score_list;
      vector<vector<int> > inter_list;
      vector<vector<string> > correct_list;

      //+ , - でCDSの情報を分け、それぞれで回す
      for(int i=0;i<CDS_list.size();i++){
            TMP =split(CDS_list[i].second,'\t');
            if(TMP[4] == "+"){CDS_list_p.push_back(CDS_list[i]);
            }else if(TMP[4] == "-"){CDS_list_m.push_back(CDS_list[i]);
            }else{cout <<"strand 情報が記載されていません。\n"; return;
            }
      }
      
      //+strand再帰関数を用いて全解探索
      if(CDS_list_p.size()!=0){DP(CDS_list_p,gro_num,num,fout,hash,hashcds,hashintron,hashinter,hashinter2,hashweight,"+", score_list,correct_list,inter_list);}
      
      //-strand再帰関数を用いて全解探索
      if(CDS_list_m.size()!=0){
            reverse(CDS_list_m.begin(),CDS_list_m.end());
            DP(CDS_list_m,gro_num,num,fout,hash,hashcds,hashintron,hashinter,hashinter2,hashweight,"-", score_list, correct_list,inter_list);
      }

      //出力
      if(correct_list.size()!=0){
            if(out ==0){
                  int q=-1;
                  if(correct_list.size()==1){q=0;}
                  else{if(score_list[0] > score_list[1]){q=0;}else{q=1;}}
                  OUTPUT(q,inter_list,correct_list,fout,gro_num,num);
            }else if(out ==1){
                  for(int p=0;p<correct_list.size();p++){
                        OUTPUT(p,inter_list,correct_list,fout,gro_num,num);      
                  }
            }
      }
      //if(num>=1){fout << "\n";}
}
void OUTPUT(int p,vector<vector<int> > inter_list,vector<vector<string> > correct_list,ofstream &fout,int gro_num,int &num)
{
      //+strandと-strandのmaxscoreを比較して、scoreの大きいほうを出力
      vector <string> tmp;
      string st,ed;
      vector<int> o;

      //内側intergenicの認識
      o.push_back(0);
      for(int j=0;j<inter_list[p].size();j++){if(inter_list[p][j] != -1){o.push_back(j);}}
      o.push_back(inter_list[p].size());
      
      for(int j=0;j<o.size()-1;j++){
            num++;
            //cout << o[j+1] << "\t"<<correct_list[p].size() << "\n";
            tmp = split(correct_list[p][o[j+1]-1],'\t');
            if(tmp[4]=="+"){ed = tmp[2];}else{ed = tmp[1];}
            tmp = split(correct_list[p][o[j]],'\t');
            if(tmp[4]=="+"){st = tmp[1];}else{st = tmp[2];}
            cout <<gro_num<<"\tdone\t"<< "\n";
            if(tmp[4]=="+"){
                  fout <<tmp[0]<<"\tcandidate\tgene\t"<<st<<"\t"<<ed<<"\t.\t"<<tmp[4]<<"\t"<<"."<<"\t"<<"ID=group_num_"<<gro_num<<"_"<<num<<"\n";
                  fout <<tmp[0]<<"\tcandidate\tmRNA\t"<<st<<"\t"<<ed<<"\t.\t"<<tmp[4]<<"\t"<<"."<<"\t"<<"ID=group_num_"<<gro_num<<"_"<<num<<".mrna1;Parent="<<"group_num_"<<gro_num<<"_"<<num<<"\n";
            }
            else{
                  fout <<tmp[0]<<"\tcandidate\tgene\t"<<ed<<"\t"<<st<<"\t.\t"<<tmp[4]<<"\t"<<"."<<"\t"<<"ID=group_num_"<<gro_num<<"_"<<num<<"\n";
                  fout <<tmp[0]<<"\tcandidate\tmRNA\t"<<ed<<"\t"<<st<<"\t.\t"<<tmp[4]<<"\t"<<"."<<"\t"<<"ID=group_num_"<<gro_num<<"_"<<num<<".mrna1;Parent="<<"group_num_"<<gro_num<<"_"<<num<<"\n";
            }
            int n=0;
            for(int k=o[j];k<o[j+1];k++){
                  n++;
                  tmp = split(correct_list[p][k],'\t');
                  fout <<tmp[0]<<"\tcandidate\texon\t"<<tmp[1]<<"\t"<<tmp[2]<<"\t.\t"<<tmp[4]<<"\t"<<tmp[3]<<"\t"<<"ID=group_num_"<<gro_num<<"_"<<num<<".mrna1.exon"<<n<<";Parent=group_num_"<<gro_num<<"_"<<num<<".mrna1\n";
            }
            n=0;
            for(int k=o[j];k<o[j+1];k++){
                  n++;
                  tmp = split(correct_list[p][k],'\t');
                  fout <<tmp[0]<<"\tcandidate\tCDS\t"<<tmp[1]<<"\t"<<tmp[2]<<"\t.\t"<<tmp[4]<<"\t"<<tmp[3]<<"\t"<<"ID=group_num_"<<gro_num<<"_"<<num<<".mrna1.cds"<<n<<";Parent=group_num_"<<gro_num<<"_"<<num<<".mrna1\n";
            }
      }
}

void DP(vector<pair<int , string > > CDS_list, int gro_num,int &num,ofstream &fout,unordered_map<string,string>hash,unordered_map<string,double> hashcds,unordered_map<string,double> hashintron,unordered_multimap<string,groupinter> hashinter,unordered_multimap<string,string> hashinter2,unordered_map<string, double> hashweight ,string w,vector<long long int> &score_list,vector< vector<string> > &correct_list,vector< vector<int> > &inter_list)
{
      //-----declare variable
      //最高scoreを記録するpathは1つに絞り込む scoreが一緒の場合、CDS長が長い方を採用、それも一緒の場合、start positionが早い方を出力
      unordered_map <int,pair<long long int ,int > > maxdp;
      int inidp[CDS_list.size() + 1],intdp[CDS_list.size() + 1],lendp[CDS_list.size() + 1]; //initialかどうかを判定するinidpとintergenicかどうかを判定するintdpを設定
      long long int dp[CDS_list.size() + 1][CDS_list.size() + 1];
      long long int ldp[CDS_list.size() + 1][CDS_list.size() + 1]; //CDSの長さを保存するためのtable
      string key;
      vector<string> tmp;
      int ter0=0,ter1=0,t=-1,path=-1; //tはtracebackするときに使用する最初のi
      long long int max=0,intergenic=0;

      //Groupの始まりと終わりのpositionを定義
      tmp=split(CDS_list[0].second,'\t');
      if(tmp[4]=="+"){ter0=stoi(tmp[1]);}else{ter1=stoi(tmp[2]);}
      tmp=split(CDS_list[CDS_list.size()-1].second,'\t');
      if(tmp[4]=="+"){ter1=stoi(tmp[2]);}else{ter0=stoi(tmp[1]);}

      //----------------------Initialize DP table
      for(int i=0;i<CDS_list.size();i++){
            //-----initialize table
            for(int j=0;j<i;j++){dp[j][i]=-1;ldp[j][i]=-1;}
            //-----initialize seq
            inidp[i]=-1;intdp[i]=-1;
            dp[CDS_list.size()][i]=-1;
            ldp[CDS_list.size()][i]=-1;
            
            //-----input exon score(ここに高橋くんのscoring関数をmergeする)
            tmp=split(CDS_list[i].second,'\t');
            key= ItoS(gro_num)+"_"+tmp[1]+"_"+tmp[2]+"_"+w+"_"+tmp[3];
            dp[i][i]=(long long int)hashcds[key];
            ldp[i][i]= (stoi(tmp[2])-stoi(tmp[1])+1);            
            lendp[i]=ldp[i][i];
            
            //-----initial exonの場合、Intergenicのscoreを保存しておく
            if(tmp[7]=="0" || tmp[7]=="3"){   //iniは1 or single
                  inidp[i]=1;
                  if(tmp[4]=="+"){
                        if(stoi(tmp[1])!=ter0){dp[CDS_list.size()][i]=(long long int)hashcds[key]+ (long long int)Score_intergenic2(gro_num, ter0, stoi(tmp[1])-1, hashweight, hashinter,hashinter2);}
                        else{dp[CDS_list.size()][i]=(long long int)hashcds[key];}
                  }else{
                        if(stoi(tmp[2])!=ter1){dp[CDS_list.size()][i]=(long long int)hashcds[key]+ (long long int)Score_intergenic2(gro_num, stoi(tmp[2])+1, ter1, hashweight, hashinter,hashinter2);}
                        else{dp[CDS_list.size()][i]=(long int)hashcds[key];}
                  }
            }
      }

      
      //Group内にSingle exon geneがただ１つのみの場合
      if(CDS_list.size()==1){t=0;max=dp[CDS_list.size()][0];}
      
      //----------------------Making DP table
      for(int i=0;i<CDS_list.size();i++){

            long long int score=0;
            int path=-1;
            //初期化
            if(dp[CDS_list.size()][i] != -1){
                  score=dp[CDS_list.size()][i];inidp[i] =1;intdp[i]=-1;
            }else{score=dp[i][i];}

            //探索
            for(int j=0;j<i;j++){
                  int st=0,ed=0;
                  int flag = filter_fuc(CDS_list[j].second,CDS_list[i].second,hash,st,ed); //計算するなら1,5,しないなら0
                                    
                  if(flag==1 ||flag==5){
                        key= ItoS(gro_num)+"_"+ ItoS(st+1)+"_"+ItoS(ed-1)+"_"+w;
                        //cout <<flag<<"\t"<<key<<"\t"<<hashintron[key] << "\n";
                        if(hashintron[key] != 0 ||flag==5){                              
                              //scoreの計算
                              if(flag ==1){
                                    dp[j][i] = maxdp[j].first + dp[i][i] +(long long int)hashintron[key];
                                    ldp[j][i] = lendp[j] + ldp[i][i];
                              }
                              else if(flag == 5){
                                    dp[j][i] = maxdp[j].first + dp[i][i] +(long long int)Score_intergenic2(gro_num, st, ed, hashweight, hashinter,hashinter2);
                                    ldp[j][i] = lendp[j] + ldp[i][i];
                              }else{cout <<"flagが指定された範囲を超えて使用されています。\n";}

                              //max scoreが更新される場合、store information about it
                              if(inidp[j]==1){
                                    if(score < dp[j][i] ||(score == dp[j][i] && lendp[i] < ldp[j][i])){
                                          path=j;score = dp[j][i];inidp[i] = inidp[j];lendp[i] =ldp[j][i];
                                          if(flag==5){intdp[i]=1;}else{intdp[i]=-1;}
                                    }}}}}
            maxdp[i]=make_pair(score,path);
            
            //terminal exonの場合、残りをIntergenicと仮定したときの値を保存しておく。
            tmp=split(CDS_list[i].second,'\t');
            if(tmp[7]=="3"||tmp[7]=="2"){
                  long long int predp = -1;
                  if(tmp[4]=="+"){
                        if(stoi(tmp[2]) != ter1){predp = maxdp[i].first +  (long long int)Score_intergenic2(gro_num,stoi(tmp[2])+1, ter1, hashweight, hashinter,hashinter2);
                        }else{predp = maxdp[i].first;}
                  }else{
                        if(stoi(tmp[1]) != ter0){predp = maxdp[i].first +  (long long int)Score_intergenic2(gro_num,ter0,stoi(tmp[1])-1, hashweight, hashinter,hashinter2);
                        }else{predp = maxdp[i].first;}
                  }
                  
                  if(inidp[i]==1){
                        if(t!=-1){
                              if(max < predp ||(max == predp && lendp[t] <  lendp[i])){
                                    t=i;max = predp;
                              }
                        }else{ //最初のtはここで無条件に格納する
                              t=i;max = predp;
                        }
                  }
            }
      }
      
//-----traceback
            //declare variable
            vector <string> correct;
            vector <int> inter;
            int st=0,ed=0,tt=0;            
            if(t != -1){                                    
                  trace_back(CDS_list,maxdp,t,correct,inter,intdp);
                  reverse(correct.begin(),correct.end());
                  reverse(inter.begin(),inter.end());
                  //格納
//-----max scoreと、Intergenicの場合のscore比較
                  vector<int>o;
                  o.push_back(0);
                  for(int j=0;j<inter.size();j++){
                        if(inter[j] != -1){o.push_back(j);}
                  }
                  o.push_back(inter.size());
                  //1. この時にIntergenicの値から内側Intergenic分引
                  for(int j=0;j<o.size()-1;j++){
                        tmp = split(correct[o[j+1]-1],'\t');
                        if(tmp[4]=="+"){ed = stoi(tmp[2]);}else{ed = stoi(tmp[1]);}
                        tmp = split(correct[o[j]],'\t');
                        if(tmp[4]=="+"){st = stoi(tmp[1]);}else{st = stoi(tmp[2]);}
                        if(st>ed){tt=st;st=ed;ed=tt;}                        
                        intergenic += (long long int)Score_intergenic(gro_num,st,ed, hashweight, hashinter);      
                  }
                  
            //max内のIntergenicの値を引
                  for(int j=0;j<o.size()-2;j++){
                        tmp = split(correct[o[j]],'\t');
                        if(tmp[4]=="+"){st = stoi(tmp[2]);}
                        else{ed = stoi(tmp[1]);}
                        tmp = split(correct[o[j+1]],'\t');
                        if(tmp[4]=="+"){ed = stoi(tmp[1]);}else{st = stoi(tmp[2]);}
                        if(st>ed){tt=st;st=ed;ed=tt;}
                        max -= (long long int)Score_intergenic2(gro_num,st,ed, hashweight, hashinter,hashinter2);     
                  }
                  
                  //2. 外側Intergenicの値をmaxから引
                  st=0;ed=0;tt=0;
                  tmp = split(correct[0],'\t');
                  if(tmp[4]=="+"){st = stoi(tmp[1]);}else{ed = stoi(tmp[2]);}
                  tmp = split(correct[correct.size()-1],'\t');
                  if(tmp[4]=="+"){ed = stoi(tmp[2]);}else{st = stoi(tmp[1]);}
                  if(st>ed){tt=st;st=ed;ed=tt;}                        
                  max -= (long long int)Score_intergenic2(gro_num,ter0,st, hashweight, hashinter,hashinter2);
                  max -= (long long int)Score_intergenic2(gro_num,ed,ter1, hashweight, hashinter,hashinter2);
                  if((double)max < (double)intergenic){t=-1;}
                  if(t!=-1){
                        score_list.push_back(max);
                        correct_list.push_back(correct);
                        inter_list.push_back(inter);
                  }
            }           
}

//trace back(maxdp[t].scoreから遡る)
void trace_back(vector<pair<int , string > > CDS_list, unordered_map <int,pair<long long int ,int> > maxdp,int t,vector <string> &correct,vector <int> &inter,int *intdp){
      vector<string> TMP;
      TMP=split(CDS_list[t].second,'\t');
      correct.push_back(CDS_list[t].second);
      inter.push_back(intdp[t]);      
      if(maxdp[t].second != -1){trace_back(CDS_list,maxdp,maxdp[t].second,correct,inter,intdp);}
};

