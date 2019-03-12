//遺伝子の間に重ならないように(遺伝子領域(intron含む))挟んでいく
//ただし、1gene 1transcript専用

#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <unordered_map>
#include "/data/tkazuki/onihitode/function.hpp"

using namespace std;

struct group
{
	int gspos;
	int gepos;
};

int main(int argc, char*argv[])
{
  //パラメータ数の取得
	if(argc != 4)
	{
		cout << "error:パラメーターの数が違います" << endl;
		cout << "./a.out [結果.gff] [間に挟みたい遺伝子.gff] [fout]" << endl;
		return 0;
	}
//ソートしたいファイルの読み込み
	ifstream ifs(argv[1]);
	if(ifs.fail())
	{
		cout << "error:gmapfilter file not open" << endl;
		return 0;
	}
//ソートしたいファイルの読み込み
	ifstream ifs2(argv[2]);
	if(ifs2.fail())
	{
		cout << "error:makergene(hawaii) file not open" << endl;
		return 0;
	}
//出力ファイル
	ofstream fout(argv[3]);
	if(fout.fail())
	{
		cout << "error:fout file not open" << endl;
		return 0;
	}

	string zero;

	string lin;

//makergene(hawaii).gffを読み込む
	vector<string> I;
	vector<group> G;
	group groups1;
	string scf;
	int spos, epos, spos2, epos2;
	int checkp = 0;
	unordered_map<string, vector<group> > hash;
	while(getline(ifs,lin) )
	{
		I = Split(lin,'\t');
		if(I.size() > 8)
		{
			if(I[2] == "CDS")
			{
				scf = I[0];
				spos = StoI(I[3]);
				epos = StoI(I[4]);
				groups1 = {spos, epos};
				auto itr = hash.find(scf);
				if( itr != hash.end() )
				{
					G = hash[scf];
					G.push_back(groups1);
					hash[scf] = G;
				}
				else
				{
					G.clear();
					G.push_back(groups1);
					hash[scf] = G;
				}
			}
		}
		fout << lin << endl;
	}

//fileを行単位で読む
	while(getline(ifs2,lin) )
	{
		I = Split(lin, '\t');
		if(I.size() > 8)
		{
			if(I[2] == "mRNA")
			{
				scf = I[0];
				spos = StoI(I[3]);
				epos = StoI(I[4]);
				checkp = 0;
				auto itr = hash.find(scf);
				if( itr != hash.end() )
				{
					G = hash[scf];
					for(int i = 0; i < G.size(); i++)
					{
						spos2 = G[i].gspos;
						epos2 = G[i].gepos;
						if(epos2 >= spos && epos >= spos2)
						{
							checkp = 1;
							break;
						}
					}
				}
				if(checkp == 0)
				{
                                      fout << lin << endl;
                                      cout << lin << endl;
				}
			}
			if(I[2] == "CDS" && checkp == 0)
			{
				fout << lin << endl;
				cout << lin << endl;
			}
		}
	}
	return 0;
}
