ROOT=$(dirname $( readlink -f $0 ))
Value=$2 #何分割したいか
Genome=$3
Weight=$4
Gtool=${ROOT}/Group_181130
Tool=${ROOT}/Searchalgo_181105
Edit=${ROOT}/181109_gff_editor

#------------Grouping
$Gtool -f $1 -w $Weight

#------------Scoring
#巨大クラスター特定個別対応
mkdir Part
A=`expr $Value / 4`
B=`expr $Value - $A`
cut -f1,2 Group.gff |sort -nrk2|uniq |head -n$A|cut -f1 > Part/part.txt
cd Part
cat part.txt|while read line
do
awk -v val=$line '{if($1 != val){print $0}}' ../Group.gff > ../tmp
awk -v val=$line '{if($1 == val){print $0}}' ../Group.gff > ./${line}.pre
nohup $Tool -fa $Genome -gff ${line}.pre -w $Weight -o ${line}.gff.result &
mv -f ../tmp ../Group.gff
done    
cd ../

#Group数読み込み
Gronum=`tail -n 1 Group.gff |cut -f1`
IN=`expr $Gronum / $B`

#分割パート
echo "start separate"

for f in `seq 0 $IN $Gronum`
do
    awk -v val=$f -v val2=`expr $f + $IN` '{if($1 > val && $1 <= val2){print $0}} ' Group.gff > ${f}.pre
done
echo "finish separate"

#統合ツールを回すパート
for f in `seq 0 $IN $Gronum`
do
    nohup $Tool -fa $Genome -gff ${f}.pre -w $Weight -o ${f}.gff.result &
done

#jobが完了したか確認
while [ `ps aux|grep "takahiro"|grep "Searchalgo"|wc -l` -gt 1 ]
do
    echo `ps aux|grep "takahiro"|grep "Searchalgo"|wc -l`
    sleep 10s
done
    	    
#------------Combining
#統合パート
echo "start Combining"    
cat *.gff.result Part/*.gff.result > Result.gff
#不要ファイル削除
rm -rf *.pre *.gff.result Part
awk '{if($3=="mRNA"||$3=="CDS"){print $0}}' Result.gff > tmp
$Edit -f tmp -num 1 -v
mv filtered.gff phase1.gff
rm tmp Result.gff
