 Denovo-based gene prediction protocol -- Itoh-lab annotation pipeline ver1.0
 
 Version: 2.0.0

 Author: Takahiro Shinoda

 ---- 概要 ----
     denovo-based遺伝子予測プロトコルは以下の2プログラムで構成されています。
         
    1. trinity.sh
       Trinityで、gene predictionを行う。fastaが出力される。
    
    2. oases.sh
       Oasesで、gene predictionを行う。fastaが出力される。
    
    3. denovo.sh
       1.2.で予測したdenovo結果をgenomeにmappingして、gffファイルを出力する。    

      ※BUSCOの結果が低い場合、denovo.shのOR_finderの代わりにTransdecoderを使う
      ※Transdecoderは次のように動かしてください
	TransDecoder.LongOrfs -m 40 -t ./CDS.fasta
     	TransDecoder.Predict --cpu $t -t ./CDS.fasta

 ---- 使用ツール ----
     各プログラム内部で以下のツールが動いています。
    
    1. trinity.sh
        - trinity
	※trinityはversionごとに出力数が大きく異なることが分かっています。
	※推薦は2.8.4 (2019/03/05現時点)
	https://www.stock-app.jp/teams/itoh-lab/search/145801/edit?q=trinity&type=all

    2. oases.sh
        - oases
	※defaultで回しています。

    3. denovo.sh
        - Transdecoder
	- Gmap
	- CD-hit-est
	- ORF_finder(別)
     また、自作スクリプト、seqkit、Gmap、Cufflinks TransDecoder付属ユーティリティが使われていますが、これらはscriptディレクトリ内に格納されています。

 ---- 使用方法 ----
    
    0. 全般
        - 自作ツールのソースコードは/src/にあります。

    1. trinity.sh
        - 以下の引数が必要です。
            $1: output file name directory <STRING>
            $2: RNA_seq_read1.fastq <STRING>
            $3: RNA_seq_read2.fastq <STRING>
            $4: the number of threads <INT>
            $5: the ammount of memory <INT> (ex. 100G)

    2. oases.sh
        - 以下の引数が必要です。
            $1: output file name directory <STRING>
            $2: RNA_seq_read1.fastq <STRING>
            $3: RNA_seq_read2.fastq <STRING>

    3. denovo.sh
        - 以下の引数が必要です。
            $1: genome.fasta <STRING> (mapping.shの出力結果)
            $2: denovo1.fasta <STRING>
            $3: denovo2.fasta <STRING>
            $4: the number of threads <INT>
          
 ---- Version更新履歴 ----
