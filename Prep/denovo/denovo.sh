###事前にインストールが必要なツール
#Transdecoder
#Trinity
#Oases
#GMAP

###必要なユーティリティ
ROOT=$(dirname `which $0`)
ORF_finder=${ROOT}/scripts/ORF_finder
Gmap_native=${ROOT}/scripts/gmap_native_to_format_converter.pl
Cufflinks_gff3=${ROOT}/scripts/cufflinks_gtf_to_alignment_gff3.pl
Cufflinks_fasta=${ROOT}/scripts/cufflinks_gtf_genome_to_cdna_fasta.pl
FASPRO=${ROOT}/scripts/FASPRO
Cdna_align=${ROOT}/scripts/cdna_alignment_orf_to_genome_orf.pl

###--USAGE--------------------------------------------------------------
#実行コマンドでどのフェーズまで実行するのか指定する
if [ $# -lt 3 ]; then
    echo "$1 ...genome.fasta $2 ...denovo.fasta $3 ...denovo2.fasta $4 ...thread "
    exit 1
fi
genome=$1 #reference配列情報
denovo1=$2 #trinityで構築した配列
denovo2=$3 #oasesで構築した配列
t=$4 #thread

ln -s $genome ./genome.fa
cp $denovo1 ./denovo1.fasta
cp $denovo2 ./denovo2.fasta

#fileが存在するかどうか
[ -f $genome ] || abort "$genome file not exist."
[ -f $denovo1 ] || abort "$denovo1 file not exist."
[ -f $denovo2 ] || abort "$denovo2 file not exist."

### -------------------- ホモロジー検索でhitするprotein情報をpredict時に使用

mkdir denovo1
mkdir denovo2

cd denovo1
#denovo1
$ORF_finder ../denovo1.fasta denovo1 120 false

cd ../denovo2
$ORF_finder ../denovo2.fasta denovo2 120 false

### -------------------- cdhit
cd ../
cat denovo1/denovo1.cds denovo2/denovo2.cds > denovo.cds
cd-hit-est -i ./denovo.cds -c 1.0 -o denovo.cd.fasta -T 10 -M 20000 &

### -------------------- GMAPでmapping

gmap_build -D . -d genome_gmap_DB genome.fa >gmap_build.lob 2>&1
gmap -S -t $t -n 1 -D . -d ./genome_gmap_DB denovo.cd.fasta > denovo.cd.fasta.gmap  2>gmap.stderr


#ファイル変換
$Gmap_native denovo.cd.fasta.gmap GTF >denovo.cd.fasta.gtf
$Cufflinks_gff3 denovo.cd.fasta.gtf > denovo.cd.fasta.gff3
$Cufflinks_fasta denovo.cd.fasta.gtf ./genome.fa >denovo_genome.fa

#Nを含む遺伝子配列を除去する
$FASPRO -f denovo_genome.fa -rmN > denovo_genome.rmN.fa
$ORF_finder denovo_genome.rmN.fa prefinal_min300 300 false
$Cdna_align prefinal_min300.gff3 denovo.cd.fasta.gff3 denovo_genome.rmN.fa >final_min300.gff3
