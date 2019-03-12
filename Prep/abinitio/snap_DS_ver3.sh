export PATH="/data/yuasa/tools/snap/:$PATH" #SNAPが使えるようにパスを通す。

PWD=`pwd`
####実行前の注意####
#Augustusのsecond.gb.train.trainのファイルを使うので、そのファイルが生成されてから実行してください。
####事前環境設定####
script=/data/yuasa/script/abinitio/script #scriptディレクトリーのパス(/data/yuasa/script/abinitio/script以下にあり)
maker=/data/yuasa/tools/maker/bin #maker内のbinディレクトリー(内部のmaker2zffを使用している。makerをインストールしていない場合は/data/yuasa/script/abinitio/scriptでも可能。)
evm=/data/yuasa/tools/EVidenceModeler-1.1.1/EvmUtils/misc #EVM内部のEvmUtils/misc以下にあるSNAP_to_GFF3.plを使用している。EVMをインストールする必要あり。
####事前入力情報####
genome= #名前をgffに合わせ、mtDNA配列を除いてリピートをマスクしたゲノム配列(もしも名前がgffと一致していない場合は「python /data/yuasa/script/tag_name_simple.py genome.fasta > out.fasta」を実行して見てください。)
gff_p= #学習セットのgffファイル
species= #最終的に出力するファイル名
work= #SNAP以下に作成したい作業ディレクトリ名
size= #Augustusのサンプル数と同じ
####################
second_dir=`cat PASS_information_for_snap_t${size}.txt` #Augustus_t1000.shのdictで指定したディレクトリのパス

mkdir -p SNAP
cd SNAP

python ${script}/down_size.py ${second_dir}/second.gb.train.train $gff_p
num_locus=`grep "LOCUS" ${second_dir}/second.gb.train.train |wc -l`
gff=${PWD}/down_size${num_locus}.gff

#出力先のディレクトリ作成
mkdir $work
cd $work

#入力データ作成のためのディレクトリ作成
mkdir -p make_data_set
cd make_data_set

#ゲノムに含まれる配列名でリスト作成
awk '{print $1}' $gff > genome.pre.txt
python ${script}/same_line_remover.py genome.pre.txt > genome.txt
rm genome.pre.txt

#改行を取り除いたfastaファイルの作成
python ${script}/rmk_fasta_maker_2018NOV11.py $genome

echo "##FASTA" > tmp.txt

#配列名ごとにMaker独自のgffフォーマットを作成
#gffを配列ごとに作成、sortして、末尾に配列をくっ付ける。
while read line
do
	#gffから対象の配列だけ抜き出す。
	python ${script}/first_gff_picker.py $gff $line > ${line}.st.tmp.gff
	#mRNAの位置でソートするために、同じ遺伝子データを一行に加工
	python ${script}/second_gff_liner.py ${line}.st.tmp.gff > ${line}.one_line.tmp.gff
	#ソート
	sort -nk 4 ${line}.one_line.tmp.gff > ${line}.sort.tmp.gff
	#一行にしていたデータを元に戻す
	python ${script}/third_gff_reset.py ${line}.sort.tmp.gff > ${line}.en.tmp.gff
	#genome配列から対象の配列のみを抽出
	python ${script}/fasta_pickup_2018NOV11.py $genome $line > ${line}.fasta
	#Maker独自のgffフォーマットを作り出す
	cat ${line}.en.tmp.gff tmp.txt ${line}.fasta > ${line}.gff
	#中間ファイルの掃除
	rm *.tmp.gff &
	#Maker形式のgffからzff形式に変換
	${maker}/maker2zff ${line}.gff
	#/data/nyuta/tools/maker/bin/maker2zff ${line}.gff
	#中間gffの削除
	rm ${line}.gff &
	#名前が混同しないようにファイル名genomeを配列名に変更
	mv genome.ann ${line}.ann &
	mv genome.dna ${line}.dna &
	#rm genome.seqs2keep
	wait
done < ./genome.txt

#ファイル加工段階で用いた一時的なファイルの消去
rm tmp.txt &
rm genome.txt &
rm *.fasta &

cd ../

#SNAPの遺伝子予測の結果（HMMとgffなど）を出力するディレクトリの作成
mkdir snap_prediction

#配列ごとにばらけていたファイルを統合し、最終ファイルを作り出す
cat make_data_set/*.ann > snap_prediction/genome.ann &
cat make_data_set/*.dna > snap_prediction/genome.dna &
wait

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
g++ ${script}/snap_evm.cpp -std=c++0x -O3
./a.out ${species}.gff ${species}.evmtmp.gff
perl ${evm}/SNAP_to_GFF3.pl ${species}.evmtmp.gff > ${species}.evm.gff
rm ${species}.evmtmp.gff

#評価gffファイル生成
awk '$3=="gene" || $3=="CDS"' ${species}.evm.gff > ${species}.evaluation.gff
sleep 10s

cd ../

#最終結果をまとめる
mkdir -p final_result
cd final_result
ln -s ../snap_prediction/${species}.gff ./
ln -s ../snap_prediction/${species}.hmm ./
ln -s ../snap_prediction/${species}.evm.gff ./
ln -s ../snap_prediction/${species}.evaluation.gff

cd ../../


#主要出力ファイル(final_result以下)
#.gff：SNAPの生結果
#.hmm：パラメーターファイル、Makerのインプットとして使用
#.evm.gff：SNAPの遺伝子予測の結果、EVMのインプットファイルとして使用
#.evaluation.gff：高橋くん作成評価ツール用ファイル

#修正2018/7/19(湯淺)：SNAPに入力するインプットファイル作成のための中間ディレクトリmake_data_setの削除(draftゲノムでscaffoldの数が多くなると中間ディレクトリのサイズが大きくなってしまうため。)
#修正日：2018/6/29(湯淺)
#作成日：2018/5/28
#2018/5/29：snapのgffからEVMのgffに変換する高橋くんのC++スクリプトを追加

#内部で動いているプログラム
#python ver.2.7
#snap(version 2006-07-28)
#Maker-2.31.10
#EVM-1.1.1