ROOT=$(dirname $( readlink -f $0 ))
Geneadd=${ROOT}/Geneadd_v2
Gtool=${ROOT}/Group_181015
Sin=${ROOT}/Singlegene_181031

if [ $# -lt 2 ]; then
    echo "1 ...all.single.gff(全result.sinの結合ファイル) 2 ...phase1.gff"
    exit 1
fi

###Grouping
$Gtool -f $1
###Predict Single CDS gene
$Sin -gff Group.gff -o single.50.gff -cov 50

##Add gene
$Geneadd $2 single.50.gff phase2.gff

