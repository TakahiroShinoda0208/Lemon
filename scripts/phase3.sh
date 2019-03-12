ROOT=$(dirname $( readlink -f $0 ))
Geneadd=${ROOT}/Geneadd_v2 

if [ $# -lt 4 ]; then
    echo "1 ...phase2.gff 2~4 ... add.gff.rmsin"
    echo "入れたい順番にファイルを入れてください"
    echo "4つ未満の場合、空ファイルを入れてください。"
    exit 1
fi

$Geneadd $1 $2 tmp
$Geneadd tmp $3 tmp2
$Geneadd tmp2 $4 tmp3
mv tmp3 phase3.gff
rm tmp tmp2