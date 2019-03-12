ROOT=$(dirname $( readlink -f $0 ))
FIX=${ROOT}/ORFframefix
Pol=${ROOT}/Polish_181207

if [ $# -lt 1 ]; then
    echo "1 ...phase4.gff"
    exit 1
fi

$FIX -i phase4.gff -o tmp.gff
$Pol -f tmp.gff
rm tmp.gff
