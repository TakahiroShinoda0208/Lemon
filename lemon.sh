###動かす前に指定
T=5
genome="/data/takahiro/work/Lemon/test/testset/genome.fa"
weight="/data/takahiro/work/Lemon/weights.txt"

###各種Inputファイル (絶対pathで指定してください)
mapping="/data/takahiro/work/Lemon/test/testset/mapping.gff"
denovo="/data/takahiro/work/Lemon/test/testset/denovo.gff"
augustus="/data/takahiro/work/Lemon/test/testset/augustus.gff"
snap="/data/takahiro/work/Lemon/test/testset/snap.gff"
homology="/data/takahiro/work/Lemon/test/testset/homology.gff"
mapp_high="/data/takahiro/work/Lemon/test/testset/mapping.hippv.gff"
augu_high="/data/takahiro/work/Lemon/test/testset/augustus.hippv.gff"
homo_high="/data/takahiro/work/Lemon/test/testset/homology.hippv.gff"
###---------------------------------------------------------------------

###内部で動かす各フェーズのツール
#ROOT=$(dirname $( readlink -f $0 ))
ROOT=$(dirname `which $0`)
phase0=${ROOT}/scripts/phase0.sh
phase1=${ROOT}/scripts/phase1.sh
phase2=${ROOT}/scripts/phase2.sh
phase3=${ROOT}/scripts/phase3.sh
phase4=${ROOT}/scripts/phase4.sh
phase5=${ROOT}/scripts/phase5.sh

###--USAGE--------------------------------------------------------------
#実行コマンドでどのフェーズまで実行するのか指定する
if [ $# -lt 1 ]; then
    echo "ファイルのパスを指定してください"
    echo "Annotationをどこまで実行するか指定してください(5が全階層実行です。)"
    echo "[Mode Option]"
    echo "0 ...phase0 1 ...phase1 2 ...phase2 3 ...phase3 4 ...phase4 5 ...Polish"
    echo ""
    echo "Annotationを各フェーズだけ流したい時、どのフェーズを実行するか指定してください"
    echo "尚、実行するにはそのフェーズより前の処理が全て終了している必要があります。"
    echo "[Mode Option]"
    echo "6 ...phase0 7 ...phase1 8 ...phase2 9 ...phase3 10 ...phase4 11 ...Polish"
    echo ""
    echo ""
    exit 1
fi

#コマンドが失敗した時にエラー出力して強制終了する関数
function abort
{
        echo "$@" 1>&2
        exit 1
}
#fileが存在するかどうか
[ -f $mapping ] || abort "$mapping file not exist."
[ -f $denovo ] || abort "$denovo file not exist."
[ -f $augustus ] || abort "$augustus file not exist."
[ -f $snap ] || abort "$snap file not exist."
[ -f $homology ] || abort "$homology file not exist."
[ -f $mapp_high ] || abort "$mapp_high file not exist."
[ -f $augu_high ] || abort "$augu_high file not exist."
[ -f $homo_high ] || abort "$homo_high file not exist."

###---------------------------------------------------------------------

### -----phase0-----
if ( test $1 -ge 0 && test $1 -lt 6 ) || test $1 -eq 6; then
$phase0 $mapping 2 mappingbase
$phase0 $denovo 2 denovobase
$phase0 $augustus 3 AUGUSTUS
$phase0 $snap 4 SNAP
$phase0 $homology 1 homology
cat *.fix5.gff > all.gff
cat *.fix5.gff.rmsin > all.gff.rmsin
cat *.fix5.gff.sin > all.gff.sin

$phase0 $mapp_high 5 mappingbasehippv
$phase0 $augu_high 5 AUGUSTUShippv
$phase0 $homo_high 6 homologyhippv

mkdir data
mv *.gff* data/
fi

### -----phase1-----
if ( test $1 -ge 1 && test $1 -lt 6 ) || test $1 -eq 7; then
echo "start phase1-----------------------------"
mkdir lemon
cd lemon
$phase1 ../data/all.gff $T $genome $weight
echo "finish phase1-----------------------------"
fi

### -----phase2-----
if ( test $1 -ge 2 && test $1 -lt 6 ) || test $1 -eq 8; then
echo "start phase2-----------------------------"
$phase2 ../data/all.gff.sin phase1.gff
echo "finish phase2-----------------------------"
fi

### -----phase3-----
if ( test $1 -ge 3 && test $1 -lt 6 ) || test $1 -eq 9; then
echo "start phase3-----------------------------"
$phase3 phase2.gff ../data/mappingbasehippv.fix5.gff.rmsin ../data/AUGUSTUShippv.fix5.gff.rmsin ../data/homologyhippv.fix5.gff.rmsin
echo "finish phase3-----------------------------"
fi

### -----phase4-----
if ( test $1 -ge 4 && test $1 -lt 6 ) || test $1 -eq 10; then
echo "start phase4-----------------------------"
$phase4 ../data/all.gff.rmsin phase3.gff
echo "finish phase4-----------------------------"
fi

### -----phase5-----
if ( test $1 -ge 5 && test $1 -lt 6 ) || test $1 -eq 11; then
echo "start phase5-----------------------------"
$phase5 phase4.gff
echo "finish phase5-----------------------------"
fi

### ----delete不必要なデータ----
rm Group.gff
rm phase4.80.candidate
rm single.50.gff
