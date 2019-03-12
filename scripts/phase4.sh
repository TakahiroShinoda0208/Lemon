ROOT=$(dirname $( readlink -f $0 ))
Geneadd=${ROOT}/Geneadd_v2 
Gtool=${ROOT}/Group_181015
Com=${ROOT}/Partialgene_181112

if [ $# -lt 2 ]; then
    echo "1 ...all.gff(全result.rmsinの結合ファイル) 2 ...phase3.gff"
    echo "入れたい順番にファイルを入れてください"
    exit 1
fi

###Grouping
$Gtool -f $1
###Predict Single CDS gene
$Com -i Group.gff -i2 $2 -o phase4.80.candidate -p 80

##Add gene
$Geneadd $2 phase4.80.candidate phase4.gff

