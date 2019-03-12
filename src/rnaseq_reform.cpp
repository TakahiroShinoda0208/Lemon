#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <algorithm>

using namespace std;

vector<string> Split(const string &s, char delim);
int StoI(string str);
string ItoS(int number);
string Replaced_V(string str, string replacestr, int replace);

int main(int argc, char*argv[])
{
  //パラメータ数の取得
	if(argc != 3)
	{
		cout << "error:パラメーターの数が違います" << endl;
		cout << "./a.out [transdecoder.gff3] [output.gff3]" << endl;
		return 0;
	}
//ファイルの読み込み
	ifstream ifs(argv[1]);
	if(ifs.fail())
	{
		cout << "error:fin file not open" << endl;
		return 0;
	}
//出力ファイル
	ofstream fout(argv[2]);
	if(fout.fail())
	{
		cout << "error:fout file not open" << endl;
		return 0;
	}

	string zero;

	string lin;

	vector<string> I, E, C;
	vector<int> P;
	string tmplin, tmplin2, genelin1, genelin2, mrnalin, sposs, eposs;
	string exonlin1, exonlin2, exonlin3, id0, id, id2, info, counts;
	int spos, epos, pos;
	int count = 0;
	int checkp = 0;
//fileを行単位で読む
	while(getline(ifs,lin) )
	{
		I = Split(lin, '\t');
		if(I.size() > 8)
		{
			if(I[2] == "mRNA")
			{
				if(checkp == 1)
				{
					spos = *min_element(P.begin(), P.end());
					sposs = ItoS(spos);
					epos = *max_element(P.begin(), P.end());
					eposs = ItoS(epos);

					tmplin = Replaced_V(genelin2, sposs, 4);
					tmplin2 = Replaced_V(tmplin, eposs, 5);
					fout << tmplin2 << endl;

					tmplin = Replaced_V(mrnalin, sposs, 4);
					tmplin2 = Replaced_V(tmplin, eposs, 5);
					fout << tmplin2 << endl;

					for(int i = 0; i < E.size(); i++)
					{
						fout << E[i] << endl;
					}
					for(int i = 0; i < C.size(); i++)
					{
						fout << C[i] << endl;
					}
				}
				count = 0;
				E.clear();
				C.clear();
				P.clear();

				string::size_type x1 = I[8].find("Parent=");
				string::size_type x2 = I[8].find(";", x1+1);
				id0 = I[8].substr(x1+7, x2-x1-7);
				info = "ID=" + id0 + ";Name=ORF";
				genelin1 = Replaced_V(lin, "gene", 3);
				genelin2 = Replaced_V(genelin1, info, 9);
				
				mrnalin = lin;
			}			
			else if(I[2] == "CDS")
			{
				count++;
				counts = ItoS(count);
				pos = StoI(I[3]);
				P.push_back(pos);
				pos = StoI(I[4]);
				P.push_back(pos);

				string::size_type a1 = I[8].find(";");
				id = I[8].substr(7, a1-7);
				string::size_type a2 = I[8].find(";Parent=");
				id2 = I[8].substr(a2);
				info = "ID=" + id + ".exon" + counts + id2;
				exonlin1 = Replaced_V(lin, "exon", 3);
				exonlin2 = Replaced_V(exonlin1, ".", 8);
				exonlin3 = Replaced_V(exonlin2, info, 9);
				E.push_back(exonlin3);
				C.push_back(lin);
				checkp = 1;
			}
		}
	}
	spos = *min_element(P.begin(), P.end());
	sposs = ItoS(spos);
	epos = *max_element(P.begin(), P.end());
	eposs = ItoS(epos);

	tmplin = Replaced_V(genelin2, sposs, 4);
	tmplin2 = Replaced_V(tmplin, eposs, 5);
	fout << tmplin2 << endl;

	tmplin = Replaced_V(mrnalin, sposs, 4);
	tmplin2 = Replaced_V(tmplin, eposs, 5);
	fout << tmplin2 << endl;

	for(int i = 0; i < E.size(); i++)
	{
		fout << E[i] << endl;
	}
	for(int i = 0; i < C.size(); i++)
	{
		fout << C[i] << endl;
	}
	return 0;
}


//split関数 使い方 Split(文字列,diliminator('\t'など)）
vector<string> Split(const string &s, char delim)
{
	vector<string> elems;
	stringstream ss(s);
	string item;
	while (getline(ss, item, delim))
	{
		if (!item.empty())
		{
			elems.push_back(item);
		}
	}
	return elems;
}

//stringをintに変える関数
int StoI(string str)
{
	int number;
	istringstream iss(str);
	iss >> number;
	return number;
}

//intをstringに変える関数
string ItoS(int number)
{
	stringstream ss;
	ss << number;
	return ss.str();
}

//あるtabのやつを入れ替える
string Replaced_V(string str, string replacestr, int replace)
{
	string output;
	vector<string> I;
	I = Split(str, '\t');
	for(int i = 0; i < replace - 1; i++)
	{
		output += I[i];
		output += '\t';
	}
	output += replacestr;
	for(int i = replace; i < I.size(); i++)
	{
		output+= '\t';
		output+= I[i];
	}
	return(output);
}
