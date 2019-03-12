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
void Grouping(ifstream &ifs,ifstream &ifs2);
void Output(ofstream &fout,vector<group> &Box,int &gro_num,unordered_map<string,int> whash,int wtotal);
void Subgroup(int &gro_num,int ele_num,vector<group> &Box,int i,int gst,int ged,unordered_map<string,int> whash,int wtotal,unordered_map<string,int> whash2,ofstream &fout);


//main関数
int main(int argc,char**argv)
{
      ifstream fin,fin2; //fin1=gff result,fin2=weight.txt
      int num = 1;
      int flag=0;
      for(int i=0;i<argc;i++){
            string ss=argv[i];
            if(ss=="-f"){fin.open(argv[i+1]);}
            if(ss=="-w"){fin2.open(argv[i+1]);}
      }

//入力ファイルが存在しなかった時の出力
      if(argc<=2){            
            cout << endl;
            cout <<"\t"<<"annotation_grouping tool"<<"\n\n\n";
            cout <<"version 1.1" <<"\n";
            cout <<"updated 2018/11/22"<<"\t"<<"TAKAHIRO SHINODA"<<"\n\n\n";

            cout <<"---How to use---"<<"\n";
            cout << "[Group mode]"<<"\t\t\t"<< "0. Group -f rna.gff -w weight.txt" <<"\n";
            cout <<"------------"<<"\n\n\n";
            return 1;
      }
      
//-----------------------------------------------------------------------------------------------------
      if(flag==0){
            cout << "START Grouping." << endl;
            Grouping(fin,fin2);
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

void Grouping(ifstream &ifs,ifstream &ifs2)
{
//----------------------fileを行単位で読みこみ、chr・stでsortする
      ofstream fout("Group.gff");
      group gff;
      vector<string> I,gff_vec;
      vector<group> Box,Box2;
      int st=0,ed=0,gro_num=1,wtotal=0;
      string lin,chr,tool,strand;


      
//----------------------weight fileの読み込み
      unordered_map<string,int> whash;
      while(getline(ifs2,lin)){
            I = split(lin,'\t');
            if(I.size()!=2){
                  cout << "weight fileが認識できませんでした。\n";
                  return;
            }
            whash[I[0]]=stoi(I[1]);
            wtotal+=stoi(I[1]);
      } 
      cout <<wtotal<<"\n";
      

//----------------------gfffileの読み込み
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
      if(Box.size()!=0){Output(fout,Box,gro_num,whash,wtotal);}
      if(Box2.size()!=0){Output(fout,Box2,gro_num,whash,wtotal);}
      
}

void Output(ofstream &fout,vector<group> &Box,int &gro_num,unordered_map<string,int> whash,int wtotal)
{
//Groupingして、mapping_baseの数を確認する(0,1,2>=)

      int ele_num=0,map_num=0,tmpst=Box[0].st,tmped=Box[0].ed;
      string tmpchr=Box[0].chr;
      int ele_w=0;
      unordered_map<string,int> whash2;
      
      for(int i=0;i<Box.size();i++){
            if(tmpchr == Box[i].chr)
            {
                  if(Box[i].st <= tmped ){
                        if(tmped < Box[i].ed ){tmped = Box[i].ed;}           
                        ele_num ++;
                        whash2[Box[i].tool]=whash[Box[i].tool];
                        //if(Box[i].tool == "mappingbase"){map_num++;}  
                  }else{
                        for(auto itr = whash2.begin(); itr != whash2.end(); ++itr) {ele_w+=itr->second;}
                        if(((double)ele_w/(double)wtotal*(double)100) >=50){ //重みの合計が半分を超える場合のみ出力
                              cout <<gro_num<<"\t"<<ele_num<<"\t"<<0<<"\tmaster"<<"\n";
                              Subgroup(gro_num,ele_num,Box,i,tmpst,tmped,whash,wtotal,whash2,fout);
                        }
                        //初期化
                        ele_num=1;map_num=0;ele_w=0;whash2.clear();
                        tmpst = Box[i].st;tmped = Box[i].ed;  
                        whash2[Box[i].tool]=whash[Box[i].tool];
                  }
                  
            }else{
                  for(auto itr = whash2.begin(); itr != whash2.end(); ++itr) {ele_w+=itr->second;}
                  if(((double)ele_w/(double)wtotal*(double)100) >=50){ //重みの合計が半分を超える場合のみ出力
                        cout <<gro_num<<"\t"<<ele_num<<"\t"<<0<<"\tmaster"<<"\n";
                        Subgroup(gro_num,ele_num,Box,i,tmpst,tmped,whash,wtotal,whash2,fout);
                  }
                  //初期化
                  ele_num=1;map_num=0;ele_w=0;whash2.clear();
                  tmpst = Box[i].st;tmped = Box[i].ed;tmpchr = Box[i].chr;
                  whash2[Box[i].tool]=whash[Box[i].tool];
            }            
      }
      
      for(auto itr = whash2.begin(); itr != whash2.end(); ++itr) {ele_w+=itr->second;}
      if(((double)ele_w/(double)wtotal*(double)100) >=50){ //重みの合計が半分を超える場合のみ出力
            cout <<gro_num<<"\t"<<ele_num<<"\t"<<0<<"\tmaster"<<"\n";
            int i=Box.size();
            Subgroup(gro_num,ele_num,Box,i,tmpst,tmped,whash,wtotal,whash2,fout);
      }
}

void Subgroup(int &gro_num,int ele_num,vector<group> &Box,int i,int gst,int ged,unordered_map<string,int> whash,int wtotal,unordered_map<string,int> whash2,ofstream &fout)
{
      //-step0.---------------変数を宣言する・初期化
      int flag=0,tmpn=0,tbit=0;
      vector< pair<int,int> > correct;
      vector<string> tmp;      
      boost::dynamic_bitset<> pret(1),pre1(1),pre2(1),pre3(1),pre4(1),pre5(1);
      pre1.resize(ged-gst+3);pre2.resize(ged-gst+3);pre3.resize(ged-gst+3);pre4.resize(ged-gst+3);pre5.resize(ged-gst+3);pret.resize(ged-gst+3);
      whash2.clear();
      
      
      //-step1.---------------Subgroup内のmRNA領域を elementに応じて格納する。
      for(int j=i-ele_num;j<i;j++){
            pret.set();pret.set(0,0);pret.set(pret.size()-1,0);
            tmp=split(Box[j].gff[0],'\t');
            int mst=stoi(tmp[3]),med=stoi(tmp[4]);
            
            pret >>=(ged-gst)-(med-mst)+1;
            pret <<=(mst-gst)+1;
            //cout << pret[0] <<"\t"<<pret[1]<<"\t"<<pret[ged-gst+1]<<"\t"<<pret[ged-gst+1+1]<<"\n";
            if(Box[j].tool.find("mapping")  !=  string::npos){pre1=(pre1|pret);
            }else if(Box[j].tool.find("denovo")  !=  string::npos){pre2=(pre2|pret);
            }else if(Box[j].tool.find("AUGUSTUS")  !=  string::npos){pre3=(pre3|pret);
            }else if(Box[j].tool.find("SNAP")  !=  string::npos){pre4=(pre4|pret);
            }else if(Box[j].tool.find("homology")  !=  string::npos){pre5=(pre5|pret);
            }
      }

      //-step2.---------------masterをpositionごとに確認し、閾値の4分の1以下の値でgroupを切断(必ず採用されない値)
      for(int j=1;j<=ged-gst+1;j++){//groupのはじからはじまで
            tbit = whash["mappingbase"]*pre1[j]+ whash["denovobase"]*pre2[j]+ whash["AUGUSTUS"]*pre3[j]+ whash["SNAP"]*pre4[j]+ whash["homology"]*pre5[j];
            if(flag==0 && tbit>(wtotal/4)){
                  tmpn=j;flag=1;
            }else if(flag==1 && tbit<=(wtotal/4)){
                  correct.push_back(make_pair(tmpn,j-1));
                  flag=0;tmpn=0;
            }
      }
      
      if(tmpn!=0){correct.push_back(make_pair(tmpn,ged-gst+1));}

      //-step3.---------------correct中に格納されているペアに属するmRNAを判定して、それを1グループとして出力
      int pele_num =0,ele_w=0,prest=0,trust=0; //prestは一つ前の終了点を保存する
      vector<group> cBox;

      for(int j=0;j<correct.size();j++){
            for(int k=i-ele_num;k<i;k++){
                   //新しく定義した領域内に入るかどうか
                  if((Box[k].st < correct[j].second+gst-1 && correct[j].first+gst-1 <=Box[k].st) || (Box[k].ed <= correct[j].second+gst-1 && correct[j].first+gst-1 <Box[k].ed) || (Box[k].st <= correct[j].first+gst-1 && correct[j].second+gst-1 <=Box[k].ed)){
                        pele_num ++;
                        cBox.push_back(Box[k]);
                        whash2[Box[k].tool]=whash[Box[k].tool];
                  }else if(correct[j].second+gst-1 < Box[k].st){break;}
            }      
            for(auto itr = whash2.begin(); itr != whash2.end(); ++itr) {ele_w+=itr->second;}
            if(((double)ele_w/(double)wtotal*(double)100) >=50){ //重みの合計が半分を超える場合のみ出力
                  cout <<gro_num<<"\t"<<pele_num<<"\t0\n";
                  for(int l=0;l<cBox.size();l++){
                        //mRNAの行を出力
                        tmp = split(cBox[l].gff[0],'\t');
                        if(tmp[6]=="+"){                              
                              fout <<gro_num<<"\t"<<pele_num<<"\t0\t"<<cBox[l].gff[0]<<"\n";
                              for(int m=1;m<cBox[l].gff.size();m++){
                                    //breakpoint内に属する一つ外側のExon・Intronまで出力する
                                    tmp = split(cBox[l].gff[m],'\t');
                                    //+の時genomeの上流
                                    if((correct[j].first+gst-1)<stoi(tmp[4]) && flag==0){
                                          flag=1;
                                          if(m!=1){fout <<gro_num<<"\t"<<pele_num<<"\t0\t"<<cBox[l].gff[m-1]<<"\n";}
                                    }
                                    if(flag==1){fout <<gro_num<<"\t"<<pele_num<<"\t0\t"<<cBox[l].gff[m]<<"\n";}
                                    if((correct[j].second+gst-1)<stoi(tmp[3]) && flag==1 && m!=cBox[l].gff.size()-1){
                                          flag=0;
                                          break;
                                    }     
                              }
                        }
                        else if(tmp[6]=="-"){
                              fout <<gro_num<<"\t"<<pele_num<<"\t0\t"<<cBox[l].gff[0]<<"\n";
                              for(int m=1;m<cBox[l].gff.size();m++){
                                    //breakpoint内に属する一つ外側のExon・Intronまで出力する
                                    tmp = split(cBox[l].gff[m],'\t');
                                    //+の時genomeの上流
                                    if(stoi(tmp[3])<(correct[j].second+gst-1) && flag==0){
                                          flag=1;
                                          if(m!=1){fout <<gro_num<<"\t"<<pele_num<<"\t0\t"<<cBox[l].gff[m-1]<<"\n";}
                                    }
                                    if(flag==1){fout <<gro_num<<"\t"<<pele_num<<"\t0\t"<<cBox[l].gff[m]<<"\n";}
                                    if(stoi(tmp[4])<(correct[j].first+gst-1) && flag==1 && m!=cBox[l].gff.size()-1){
                                          flag=0;
                                          break;
                                    }     
                              }
                        } 
                  }
                  gro_num++;
                  //初期化
                  pele_num=0;ele_w=0;whash2.clear();cBox.clear();
            }
      }
}


