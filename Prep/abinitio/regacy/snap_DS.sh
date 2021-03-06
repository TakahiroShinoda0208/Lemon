export PATH="/data/yuasa/tools/snap/:$PATH" #SNAPが使えるようにパスを通す。

PWD=`pwd`

####実行前の注意####
#Augustusのsecond.gb.train.trainのファイルを使うので、そのファイルが生成されてから実行してください。
####事前入力####
genome=  #名前をgffに合わせ、mtDNA配列を除いてリピートをマスクしたゲノム配列(もしも名前がgffと一致していない時、「Scaffold名+空白+冗長な情報」のような形を取っていた場合は「python /data/yuasa/script/tag_name_simple.py genome.fasta > out.fasta」を実行して見てください。)
gff_p= #学習セットのgffファイル
species= #最終的に出力するファイル名
work= #作成したい作業ディレクトリ名
size= #学習セットのサイズ(Augustusで使ったもの)
################

second_dir=`cat PASS_information_for_snap_t${size}.txt` #Augustus_t1000.shのdictで指定したディレクトリのパス

mkdir -p SNAP
cd SNAP

python /data/yuasa/script/snap/down_size.py ${second_dir}/second.gb.train.train $gff_p
num_locus=`grep "LOCUS" ${second_dir}/second.gb.train.train |wc -l`
gff=${PWD}/down_size${num_locus}.gff

#出力先のディレクトリ作成
mkdir $work
cd $work

#入力データ作成のためのディレクトリ作成
mkdir -p make_data_set
cd make_data_set

#ゲノムに含まれる配列名でリスト作成
python /data/yuasa/script/snap/make_fasta_list.py $genome > genome.txt

echo "##FASTA" > tmp.txt

#配列名ごとにMaker独自のgffフォーマットを作成
#gffを配列ごとに作成、sortして、末尾に配列をくっ付ける。
while read line
do
	#gffから対象の配列だけ抜き出す。
	python /data/yuasa/script/snap/first_gff_picker.py $gff $line > ${line}.st.tmp.gff
	#mRNAの位置でソートするために、同じ遺伝子データを一行に加工
	python /data/yuasa/script/snap/second_gff_liner.py ${line}.st.tmp.gff > ${line}.one_line.tmp.gff
	#ソート
	sort -nk 4 ${line}.one_line.tmp.gff > ${line}.sort.tmp.gff
	#一行にしていたデータを元に戻す
	python /data/yuasa/script/snap/third_gff_reset.py ${line}.sort.tmp.gff > ${line}.en.tmp.gff
	#genome配列から対象の配列のみを抽出
	python /data/yuasa/script/snap/fasta_pickup.py $genome $line > ${line}.fasta
	#Maker独自のgffフォーマットを作り出す
	cat ${line}.en.tmp.gff tmp.txt ${line}.fasta > ${line}.gff
	#中間ファイルの掃除
	rm *.tmp.gff
	rm *.fasta
	#Maker形式のgffからzff形式に変換（中村さんのディレクトリのMakerの中のスクリプト使用）
	/data/nyuta/tools/maker/bin/maker2zff ${line}.gff
	#中間gffの削除
	rm ${line}.gff
	#名前が混同しないようにファイル名genomeを配列名に変更
	mv genome.ann ${line}.ann
	mv genome.dna ${line}.dna
	#rm genome.seqs2keep
done < ./genome.txt

#ファイル加工段階で用いた一時的なファイルの消去
rm tmp.txt
rm genome.txt

cd ../

#SNAPの遺伝子予測の結果（HMMとgffなど）を出力するディレクトリの作成
mkdir snap_prediction

#配列ごとにばらけていたファイルを統合し、最終ファイルを作り出す
cat make_data_set/*.ann > snap_prediction/genome.ann
cat make_data_set/*.dna > snap_prediction/genome.dna

#データセット作成用のディレクトリを削除
rm -r make_data_set

cd snap_prediction/

#学習セットの遺伝子の特徴を記述[スキップ可]
fathom genome.ann genome.dna -gene-stats > gene-stats.log 2>&1
#学習セットとして不適当な可能性のある遺伝子のエラーログをつける。[スキップ可]
fathom genome.ann genome.dna -validate > validate.log 2>&1
#一領域一遺伝子となるようにカットする。
fathom genome.ann genome.dna -categorize 1000  > categorize.log 2>&1
#方向をplusに統一する。
fathom uni.ann uni.dna -export 1000 -plus> uni-plus.log 2>&1
#遺伝子予測用のファイル構築
forge export.ann export.dna > forge.log 2>&1
#学習
hmm-assembler.pl ${species} . > ${species}.hmm
#遺伝子予測(gff作成;出力ファイルはCDS )
snap ${species}.hmm $genome -gff > ${species}.gff

#SNAPのgffをEVMのgffに変換（高橋くん作成）
g++ /data/tkazuki/tool/EVidenceModeler-1.1.1/EvmUtils/misc/snap_evm.cpp -std=c++0x -O3
./a.out ${species}.gff ${species}.evmtmp.gff
perl /data/tkazuki/tool/EVidenceModeler-1.1.1/EvmUtils/misc/SNAP_to_GFF3.pl ${species}.evmtmp.gff > ${species}.evm.gff
rm ${species}.evmtmp.gff

#評価gffファイル生成
awk '$3=="gene" || $3=="CDS"' ${species}.evm.gff > ${species}.evaluation.gff
sleep 10s

cd ../

#最終結果をまとめる
mkdir -p final_result
cd final_result
ln -s ../snap_prediction/${species}.hmm ./
ln -s ../snap_prediction/${species}.evm.gff ./
ln -s ../snap_prediction/${species}.evaluation.gff

cd ../../


#主要出力ファイル(final_result以下)
#	${species}.hmm：パラメーターファイル、Makerのインプットとして使用
#	${species}.evm.gff：SNAPの遺伝子予測の結果、EVMのインプットファイルとして使用
#	${species}.evaluation.gff：高橋くん作成評価ツール用ファイル

#修正2018/7/19(湯淺)：SNAPに入力するインプットファイル作成のための中間ディレクトリmake_data_setの削除(draftゲノムでscaffoldの数が多くなると中間ディレクトリのサイズが大きくなってしまうため。)
#修正日：2018/6/29(湯淺)
#作成日：2018/5/28
#2018/5/29：snapのgffからEVMのgffに変換する高橋くんのC++スクリプトを追加
#pythonはver.2.7で作成
