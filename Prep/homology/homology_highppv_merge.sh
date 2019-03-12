#homology.shを動かすのに必要なscriptの位置
SCRIPT=/data/tkazuki/annotation/homologypipeline/homologyscript/

#コマンドが失敗した時にエラー出力して強制終了する関数
function abort
{
	echo "$@" 1>&2
	exit 1
}

#正しい引数かどうか
if [ $# -ne 2 ]; then
	echo "error::実行するには3個の引数が必要です。" 1>&2
	echo "./homology_highppv_merge.sh [catしたspalnresult_over90.gff3] [output.gff3]" 1>&2
	exit 1
fi

cp $SCRIPT/function.hpp .

#fileが存在するかどうか
[ -f $1 ] || abort "$1 file not exist."

echo "start"
echo "$0 $1 $2"

g++ $SCRIPT/spaln_reform.cpp -std=c++0x -O3
./a.out $1 homologytmp.gff3

awk '$3=="mRNA" || $3=="CDS"' homologytmp.gff3 > homologytmp2.gff3

$SCRIPT/181109_gff_editor -f homologytmp2.gff3 -l
mv longest.gff $2

rm a.out homologytmp.gff3 homologytmp2.gff3 function.hpp
echo "finish"
