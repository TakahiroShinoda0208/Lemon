//ver.2 intergenicのスコアの考え方を変更
#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <unordered_map>

using namespace std;

//関数
vector<string> Split(const string &s, char delim);
int StoI(string str);
string ItoS(int number);
double StoD(string str);

//構造体
struct tgroup
{
	double weight;
	int spos;
	int epos;
};

//weight fileの読み込みhash
unordered_map<string, double> Weight_hash(ifstream &file1)
{
	file1.clear();
	file1.seekg(0, ios_base::beg);
	string lin;
	vector<string> I;
	unordered_map<string, double> hash;
	while(getline(file1, lin) )
	{
		I = Split(lin, '\t');
		hash[I[0]] = StoD(I[1]);
	};
	return hash;
}

//CDS用のscore hash
//file1:grouping file
unordered_map<string, double> Score_cds_hash(istream &file1, unordered_map<string, double> hashweight)
{
	file1.clear();
	file1.seekg(0, ios_base::beg);
	string lin;
	vector<string> I;
	double weight;
	string key, strand, tool, key2;
	int cdsspos, cdsepos, groupnum;
	tgroup cdsinfo;
	unordered_multimap<int, tgroup> cdshash;
	unordered_map<string, int> cdshash2;
	vector<string> K;
	//grouping fileのinput
	while(getline(file1, lin) )
	{
		I = Split(lin, '\t');
		if(I[5] == "mRNA")
		{
			groupnum = StoI(I[0]);
			weight = hashweight[I[4]];
			strand = I[9];
			tool = I[4];
		}
		else if(I[5] == "CDS")
		{
			cdsspos = StoI(I[6]);
			cdsepos = StoI(I[7]);

			//同じものはtoolごとに1回しか出力しない　key2 groupnum_tool_spos_epos_strand_frame
			key2 = I[0] + "_" + tool + "_" + I[6] + "_" + I[7] + "_" + strand + "_" + I[10];
			auto itr = cdshash2.find(key2);
			if( itr == cdshash2.end() )
			{
				//groupの情報の格納 key groupnum_cdsspos_cdsepos_strand_frame
				key = I[0] + "_" + I[6] + "_" + I[7] + "_" + strand + "_" + I[10];
				K.push_back(key);
				cdsinfo = {weight, cdsspos, cdsepos};
				cdshash.insert(pair<int, tgroup>(groupnum, cdsinfo));
				cdshash2[key2] = 0;
			}
		}
	}

	//score格納用のhash作成
	double tmpscore, score;
	unordered_map<string, double> hash;
	for(int i = 0; i < K.size(); i++)
	{
		I = Split(K[i], '_');
		groupnum = StoI(I[0]);
		cdsspos = StoI(I[1]);
		cdsepos = StoI(I[2]);
		score = 0;

		auto range = cdshash.equal_range(groupnum);
		for (auto iterator = range.first; iterator != range.second; iterator++)
		{
			auto target = *iterator;
			//case1	- < - > -
			if(target.second.spos <= cdsspos && cdsepos <= target.second.epos)
			{
				tmpscore = (double)(cdsepos - cdsspos + 1) * target.second.weight;
				score += tmpscore;
			}
			//case2	 < - >
			else if(cdsspos <= target.second.spos && target.second.epos <= cdsepos)
			{
				tmpscore = (double)(target.second.epos - target.second.spos + 1) * target.second.weight;
				score += tmpscore;
			}
			//case3 - < - >
			else if(target.second.spos <= cdsspos && cdsspos <= target.second.epos)
			{
				tmpscore = (double)(target.second.epos - cdsspos + 1) * target.second.weight;
				score += tmpscore;
			}
			//case4 < - > -
			else if(target.second.spos <= cdsepos && cdsepos <= target.second.epos)
			{
				tmpscore = (double)(cdsepos - target.second.spos + 1) * target.second.weight;
				score += tmpscore;
			}
		}

		hash[K[i]] = score;
	}
	return hash;
}

//intron用のscore hash
//file1:weight file, file2:grouping file
unordered_map<string, double> Score_intron_hash(istream &file1, unordered_map<string, double> hashweight)
{
	file1.clear();
	file1.seekg(0, ios_base::beg);
	string lin, tool;
	vector<string> I;
	bool check = false;
	double weight;
	string strand, frame, key, key2, intronsposs, introneposs;
	int cdsspos, cdsepos, intronspos, intronepos, groupnum;
	tgroup introninfo;
	unordered_multimap<int, tgroup> intronhash;
	unordered_map<string, int> intronhash2;
	vector<string> K;
	//grouping fileのinput
	while(getline(file1, lin) )
	{
		I = Split(lin, '\t');
		if(I[5] == "mRNA")
		{
			groupnum = StoI(I[0]);
			weight = hashweight[I[4]];
			strand = I[9];
			tool = I[4];
			check = false;
		}
		else if(I[5] == "CDS")
		{
			//first exon
			if(check == false)
			{
				cdsspos = StoI(I[6]);
				cdsepos = StoI(I[7]);
				if(strand == "+")
				{
					intronspos = cdsepos + 1;
				}
				else
				{
					intronepos = cdsspos - 1;
				}
				check = true;
			}
			//それ以外
			else
			{
				cdsspos = StoI(I[6]);
				cdsepos = StoI(I[7]);
				if(strand == "+")
				{
					intronepos = cdsspos - 1;
				}
				else
				{
					intronspos = cdsepos + 1;
				}
				//同じtoolから1回しか出力しない
				intronsposs = ItoS(intronspos);
				introneposs = ItoS(intronepos);
				key2 = I[0] + "_" + tool + "_" + intronsposs + "_" + introneposs + "_" + strand;

				auto itr = intronhash2.find(key2);
				if( itr == intronhash2.end() )
				{
					//key groupnum_intronspos_intronepos_strand
					key = I[0] + "_" + intronsposs + "_" + introneposs + "_" + strand;
					K.push_back(key);

					introninfo = {weight, intronspos, intronepos};
					intronhash.insert(pair<int, tgroup>(groupnum, introninfo));
					
					intronhash2[key2] = 0;
				}

				if(strand == "+")
				{
					intronspos = cdsepos + 1;
				}
				else
				{
					intronepos = cdsspos - 1;
				}
			}
		}
	}

	//score格納用のhash作成
	double tmpscore, score;
	unordered_map<string, double> hash2;
	for(int i = 0; i < K.size(); i++)
	{
		I = Split(K[i], '_');
		groupnum = StoI(I[0]);
		intronspos = StoI(I[1]);
		intronepos = StoI(I[2]);
		score = 0;

		auto range = intronhash.equal_range(groupnum);
		for (auto iterator = range.first; iterator != range.second; iterator++)
		{
			auto target = *iterator;
			if(intronspos == target.second.spos && intronepos == target.second.epos)
			{
				tmpscore = (double)(intronepos - intronspos + 1) * target.second.weight;
				score += tmpscore;
			}
		}
		hash2[K[i]] = score;
	}
	return hash2;
}


//intergenic用の構造体
struct groupinter
{
	int gspos;
	int gepos;
};

//intergenic領域計算用のhash	file1:group
//key:groupnum_stat	value:計算用データの格納(spos, epos)
unordered_multimap<string, groupinter> Score_intergenic_hash(istream &file1)
{
	file1.clear();
	file1.seekg(0, ios_base::beg);

	string lin, key;
	int spos, epos;
	vector<string> I;
	groupinter groups2;
	unordered_multimap<string, groupinter> hash;

	//group file
	while(getline(file1, lin) )
	{
		I = Split(lin, '\t');	
		key = I[0] + "_" + I[4];
		if(I[5] == "mRNA")
		{
			spos = StoI(I[6]);
			epos = StoI(I[7]);
			groups2 = {spos, epos};
			hash.insert(pair<string, groupinter>(key, groups2));
		}
	}
	return hash;
}

//intergenic領域計算用のhash2	file1:group
//key:groupnum_stat	value:group内に存在するtool名
unordered_multimap<string, string> Score_intergenic_hash2(istream &file1)
{
	file1.clear();
	file1.seekg(0, ios_base::beg);

	string lin;
	unordered_multimap<string, string> hash;
	vector<string> I;
	int checkp = 0;
	//group file
	while(getline(file1, lin) )
	{
		I = Split(lin, '\t');
		size_t count = hash.count(I[0]);
		//keyが存在する場合
		if(count > 0)
		{
			checkp = 0;
			auto range = hash.equal_range(I[0]);
			for (auto iterator = range.first; iterator != range.second; iterator++)
			{
				auto target = *iterator;
				if(target.second == I[4])
				{
					checkp = 1;
				}
			}
			if(checkp == 0)
			{
				hash.insert(pair<string, string>(I[0], I[4]));
			}
		}
		//keyが存在しない場合
		else
		{
			hash.insert(pair<string, string>(I[0], I[4]));
		}
	}
	return hash;
}

//intergenic領域の計算(全toolでの計算)
double Score_intergenic(int groupnum, int spos, int epos, 
		unordered_map<string, double> hashweight, 
		unordered_multimap<string, groupinter> hash)
{
	string groupnums, key;
	int spos2, epos2, vsize, count;
	double weight, tmpscore;
	double score = 0;
	vector<bool> V;
	groupnums = ItoS(groupnum);
	for(auto itr = hashweight.begin(); itr != hashweight.end(); ++itr)
	{
		key = groupnums + "_" + itr->first;
		//空vector作成
		V.clear();
		for(int i = spos; i <= epos; i++)
		{
			V.push_back(false);
		}
		//cout << "check" << '\t'  << V.size() << endl;
		vsize = V.size();
		//あったらfalse -> trueに
		auto range = hash.equal_range(key);
		for (auto iterator = range.first; iterator != range.second; iterator++)
		{
			auto target = *iterator;
			//cout << target.second.gspos << '\t' << target.second.gepos << '\t' << spos << '\t' << epos << endl;
			spos2 = target.second.gspos - spos;
			epos2 = target.second.gepos - spos;
			//cout << spos2 << '\t' << epos2 << endl;
			if(spos2 >= 0)
			{
				for(int i = spos2; i <= epos2; i++)
				{
					if(i < vsize)
					{
						V[i] = true;
					}
				}
			}
			else if(epos2 >= 0)
			{
				for(int i = 0; i <= epos2; i++)
				{
					if(i < vsize)
					{
						V[i] = true;
					}
				}
			}
		}

		count = 0;
		//計算
		for(int i = 0; i < V.size(); i++)
		{
			if(V[i] == false)
			{
				count++;
			}
		}
		weight = hashweight[itr->first];
		tmpscore = (double)count * weight;
		//cout << "tmpscore" <<  '\t' << tmpscore << endl;
		score += tmpscore;
	}
	return score;
}


//intergenic領域の計算(group内に存在するtoolでの計算)
double Score_intergenic2(int groupnum, int spos, int epos, 
		unordered_map<string, double> hashweight, 
		unordered_multimap<string, groupinter> hash,
		unordered_multimap<string, string> hash2)
{
	string groupnums, key;
	int spos2, epos2, vsize, count;
	double weight, tmpscore;
	double score = 0;
	vector<bool> V;
	groupnums = ItoS(groupnum);
	for(auto itr = hashweight.begin(); itr != hashweight.end(); ++itr)
	{
		auto range = hash2.equal_range(groupnums);
		for (auto iterator = range.first; iterator != range.second; iterator++)
		{
			auto target = *iterator;
			if(target.second == itr->first)
			{
				key = groupnums + "_" + itr->first;
				//空vector作成
				V.clear();
				for(int i = spos; i <= epos; i++)
				{
					V.push_back(false);
				}
				//cout << "check" << '\t'  << V.size() << endl;
				vsize = V.size();
				//あったらfalse -> trueに
				auto range = hash.equal_range(key);
				for (auto iterator = range.first; iterator != range.second; iterator++)
				{
					auto target = *iterator;
					//cout << target.second.gspos << '\t' << target.second.gepos << '\t' << spos << '\t' << epos << endl;
					spos2 = target.second.gspos - spos;
					epos2 = target.second.gepos - spos;
					//cout << spos2 << '\t' << epos2 << endl;
					if(spos2 >= 0)
					{
						for(int i = spos2; i <= epos2; i++)
						{
							if(i < vsize)
							{
								V[i] = true;
							}
						}
					}
					else if(epos2 >= 0)
					{
						for(int i = 0; i <= epos2; i++)
						{
							if(i < vsize)
							{
								V[i] = true;
							}
						}
					}
				}

				count = 0;
				//計算
				for(int i = 0; i < V.size(); i++)
				{
					if(V[i] == false)
					{
						count++;
					}
				}
				weight = hashweight[itr->first];
				tmpscore = (double)count * weight;
				//cout << "tmpscore" <<  '\t' << tmpscore << endl;
				score += tmpscore;
			}
		}
	}
	return score;
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

//stringをdoubleに変える関数
double StoD(string str)
{
	double number;
	istringstream iss(str);
	iss >> number;
	return number;
}

// //intをstringに変える関数
// string ItoS(int number)
// {
// 	stringstream ss;
// 	ss << number;
// 	return ss.str();
// }
