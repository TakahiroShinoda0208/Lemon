 Mapping-based gene prediction protocol -- Itoh-lab annotation pipeline ver1.0
 
 Version: 2.2.2

 Author: Fumiya Kobayashi

 - 概要 -
     mapping-based遺伝子予測プロトコルは以下の4プログラムで構成されています。
    
    1. mapping.sh
        mapping -> gene predictionを行う。.gtfが出力される。
    
    2. to_merge.sh
        EVM, MAKERに与える統合用遺伝子データセットを作成する。
    
    3. to_learn_1st.sh
        Augustus, SNAPに与える学習モデル用遺伝子データセットを作成する（第1段階）。
    
    4. to_learn_2nd.sh
        Augustus, SNAPに与える学習モデル用遺伝子データセットを作成する（第2段階）。
        Repeat Maskerの結果を必要とする。

 - 使用ツール -
     各プログラム内部で以下のツールが動いています。
    
    1. mapping.sh
        - hisat2
        - hisat2-build
        - samtools
        - stringtie

    3. to_learn_1st.sh
        - gffread
        - cd-hit
    
     また、自作スクリプト、seqkit、TransDecoder付属ユーティリティが使われていますが、これらはutilディレクトリ内に格納されています。

 - 使用方法 -
    
    0. 全般
        - utilディレクトリは絶対パス参照のため基本的に変更不要ですが、別環境で実行する場合はscript冒頭のutil=以下を書き換えて下さい。
        - utilディレクトリをコピーする場合は、同一階層にPerlLibディレクトリを置いて下さい
        - 自作ツールのソースコードはutil/src/にあります。
        - samtoolsのversionが古いと正常に動作しない可能性があります。

    1. mapping.sh
        - 以下の引数が必要です。
            $1: genome.fasta <STRING>
            $2: RNA_seq_read1.fastq <STRING>
            $3: RNA_seq_read2.fastq <STRING>
            $4: prefix in output file name <STRING>
            $5: the number of threads <INT>

        - 適宜、使用ツールのパスを書き換えて下さい。
            HISAT2=${hisat2}
            HISAT2BUILD=${hisat2-build}
            SAMTOOLS=${samtools}
            STRINGTIE=${stringtie}

        - 以下のファイルが作成されます。
            ${4}.gtf
    
    2. to_merge.sh
        - 以下の引数が必要です。
            $1: stringtie.gtf <STRING> (mapping.shの出力結果)
            $2: genome.fasta <STRING>>
            $3: prefix in output file name <STRING>

        - 変数$MINを書き換えることで最小ORF長を変更できます（defaultでは90base）
        
        - 以下のファイルが作成されます。
            ${3}_merge.gff3（EVM, MAKERに渡すファイル）

    3. to_learn_1st.sh
        - version 1.x.0ではto_merge.shの結果を必要としていましたが、本versionでは独立に動作するように変更しました。
        
        - 以下の引数が必要です。
            $1: stringtie.gtf <STRING> (mapping.shの出力結果)
            $2: genome.fasta <STRING>>
            $3: prefix in output file name <STRING>
        
        - 適宜、使用ツールのパスを書き換えて下さい。
            GFFREAD=${gffread}
            CDHIT=${cd-hit}
         
        - 変数$MINを書き換えることで最小ORF長を変更できます（defaultでは300base）

        - 以下のファイルが作成されます。
            ${5}_learn_1st.gff3

    4. to_learn_2nd.sh
        - 以下の引数が必要です。
            $1: learn_1st.gff3 <STRING>
            $2: genome.fasta <STRING>
            $3: genome.fa.out <STRING>
                genome.fa.outはRepeatMaskerの出力結果に含まれるファイルです。
            $4: prefix in output file name <STRING>

        - 以下のファイルが作成されます。
            ${4}_learn_2nd.gff3（Augustus, SNAPの学習モデル用遺伝子セット）
 
 - Version更新履歴 -
    v2.2.2 (2019/02/20)
        - hisat2, samtoolsを/scratch以下で実行するように変更
    v2.2.1
        - genome配列にlower caseが使われていた場合にORFが正しく検出できていなかった問題を修正

    v2.2.0
        - 学習用データセット作成の際に、mRNA配列全体とrepeat領域が少しでも重なる配列を除外する -> exon領域とrepeat領域が少しでも重なる配列を除外するに変更

    v2.1.0
        - 統合用データセットに不完全なORF配列が出力される問題を修正

    v2.0.0
        - ORF最小長を調整
        - TransDecoderとORFfinderで比較評価
        - バグ修正
