 ab initio protocol -- Itoh-lab annotation pipeline ver1.0
 
 Version: 1.0.0

 Author: Hideaki Yuasa (湯淺 英知)

 - 概要 -
     ab initioによる遺伝子予測は以下の2プログラムで構成されています。
    
    1. Augustus_DS_ver2.sh
        Augustusによるab initio
	最終出力結果はAugustus_abinitio.gtf
    
    2. snap_DS_ver3.sh
        SNAPによるab initio
	最終出力ファイルはfinal_resultに以下4つできる
	・.gff：SNAPの生結果
	・.hmm：パラメータファイル→Makerのインプットファイルとして使用可能
	・.evm.gff：SANPの結果をEVM用の入力ファイルに変換したもの
	・.evaluation.gff：高橋くん(2019年3月卒業)作成の評価ツール用に変換したもの

 - 使用ツール -
     各プログラム内部で以下のツールが動いています。
    
    1. Augustus_DS_ver2.sh
        - python2系
        - Augustus

    3. snap_DS_ver3.sh
        - python2系 (3系では動きません)
        - SNAP
	- MAKER (インストールしなくても使用可能です。このファイルの62行目を確認してください。)
	- EVM (EVM用の出力結果が必要なければインストールする必要はありません。このファイルの63行目を確認してください。)
    
     自作スクリプトと必要なMAKERのプログラムはscriptディレクトリの下にあります。

 - 使用方法 -
    
    1. Augustus_DS_ver2.sh
        - スクリプト内の変数に以下の情報を書き込む必要があります。
            script: scriptディレクトリのパス
            spec: Augustus内部に作成されるパラメーターディレクトリの名前($AUGUSTUS_CONFIG_PATH/speciesの下に同じ名前ができないように注意)
            genome: シンプルリピートと複雑性の低いリピート以外のリピートをマスクしたゲノム配列(fasta形式)
            gff: gff形式の学習セット(配列の名前はgenomeのfastaと対応している必要がある)
            size: 学習セットサイズ(処理の過程で評価用に100遺伝子抜くのとAugusutus内部のフィルターをかけるので、遺伝子数が減ることに注意)
　　　　　　optimize_num: thred数(並列化を行うにはperlのモジュールであるParallel::ForkManagerをインストールする必要があり)
	    nice_value: Nice値の変更
		
        - 必要に応じて以下のログファイルを確認してください。
	    train.err: (エラーは多少出てきますが、)かなり多くのエラーが見られた場合は--stopCodonExcludedFromCDS=のオプションを追加してみてください。
            test_result1.txt & test_result2.txt: 1よりも2の方が値が向上していれば最適化が上手く行われています。
        
	- 並列化するとI/Oを食うのでScrach上で行うことを推奨します。
    
    2. snap_DS_ver3.sh
        - Augustus_DS_ver2.shと同じディレクトリで実行してください。

	- Augustus_DS_ver2.shを先に実行してsecond.gb.train.trainのファイルが生成されたのを確認してから実行してください。(学習セットを共有しているためです。)

        - スクリプト内の変数に以下の情報を書き込む必要があります。
            script: scriptディレクトリのパス
            maker: makerのbinディレクトリーのパス(必要なプログラムはscript以下にもあるので、$scriptと同じディレクトリパスでも可)
            evm: EVM内部のEvmUtils/miscのパス(.evm.gffの生成に必要だが、必要なければこの行と109行目から113行目と126行目をコメントアウト)
	    genome: シンプルリピートと複雑性の低いリピート以外のリピートをマスクしたゲノム配列(fasta形式)
            gff_p: gff形式の学習セット(配列の名前はgenomeのfastaと対応している必要がある)
	    species: 最終的に出力するファイル名
	    work: 作業ディレクトリー名
            size: Augustus_DS_ver2.shと同じにする
 
 - Version更新履歴 -
    v1.0.0
        - 2019/MAR/07作成
