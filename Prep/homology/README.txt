
--homology.sh, homology_highppv_merge.sh の使い方--


-注意事項

	homology.sh, homology_highppv_merge.shは少なくとも自分のディレクトリにコピーしてください。

	できれば、僕が勝手にフォルダの構造を変えて、このツールを動かすために必要なスクリプトのリンクが崩れる可能性もありますので、
	このフォルダ内にあるhomologyscriptフォルダも自分のディレクトリ内にコピーして、
	homology.sh, homology_highppv_merge.sh 内の
		SCRIPT=/data/tkazuki/annotation/homologypipeline/homologyscript/
	を自分のコピーしたやつにきりかえることをおすすめします。

	できれば、入力に用いるfastaのscaffold名、gene名などの名前は単純にしてください。バグの原因になる可能性があります。	
	ある程度修正プログラムは入れていますが, "|"とかが入るとバグが起きる例も確認されてます。
	例)
	>lcl|gene00001 [chromosome=I]
	→ >gene00001
	
	読み書きが激しいので、/scratch上で動かすのが良いかと思います

-動かし方

0. spaln, mafftのインストール
	以下のようにpathを通してください
	例)
		export PATH=$PATH:/data/tkazuki/tool/spaln/spaln-master/bin/
		export PATH=$PATH:/data/tkazuki/tool/mafft-7.305-with-extensions/bin/

1. 近縁種生物のタンパク質配列のダウンロード
	ensembleなどからタンパク質配列をダウンロード
	ミトコンドリア配列などは各自で抜いてください

2. homology.shを動かす
	./homology.sh [genome.fa] [related species protein.faa]
			[spalndatabase] [thread] [超近縁 -> 1,　近縁 -> 0] [prefix]

		[genome.fa]
			遺伝子予測するgenomeのfastaのpathを指定してください

		[related species protein.faa]
			近縁生物のタンパク質配列のpathを指定してください

		[spalndatabase]
			spaln toolのspaln-master/table/のなかから
			自分の予測したい生物種のパラメータを取得してきてください
			詳細はspaln-master/table/の中のgnm2tabに書いてあり、
			自分の生物種もしくは最近縁の生物種のパラメータを選ぶと良いと思います。
				例)NematodC, InsectDm
		
		[thread]
			spalnを動かす際に用いるthread数です。
			serverの空き具合等をみて設定してください。

		[超近縁 -> 1,　近縁 -> 0]
			超近縁な生物種(オニヒトデでいうアカオニヒトデ的な生物)はalignmentの50%のfilterをかけてます
			近縁生物種はalignmentの30%のfilterをかけてます

		[prefix]
			出力とかに使います。自分のファイルが上書きされないように設定してください

3.出力file
	[prefix]_spalnresult.gff
		spalnの出力結果にfilterを掛ける前の結果です。
		spalnの生の結果に終止コドンを後ろに加えるという加工をした結果になってます。
	[prefix]_spalnresult_filter.gff
		spalnの結果にidentityなどのfilterをかけた結果です。
		EVMにはこちらの方を統合することをおすすめします。
		篠田さんのannotation pipelineのphase1に用いるのもこちらのデータです。
			(篠田さんの方の場合はcatで各生物種のgffを繋げる形になるかと思います)
	[prefix]_spalnresult_over90.gff
		篠田さんのannotation pipelineに用いる高精度遺伝子の結果です。
		identity90%のfilterがかかってます。
	alignmentresult.txt
		中間fileです。その遺伝子と近縁タンパク質がどれくらいのidentityかが書いてあります。
		消していただいても構いません。

例)
	nohup ./homology.sh c_elegance.fa c_briggsae.faa NematodC 8 1 cbriggsae &
	nohup ./homology.sh drosophila.fa c_elegance.faa InsectDm 8 0 drosophila &


4.(篠田さん作成annotation pipelineを用いる場合)
	homologyの高精度遺伝子セットが必要になります
	4.1
		各生物種ごとにでた[prefix]_spalnresult_over90.gffをcatでつなげる
	4.2
		./homology_highppv_merge.sh [catした[prefix]_spalnresult_over90.gff] [出力file名]

例)
	cat acp_spalnresult_over90.gff spu_spalnresult_over90.gff scu_spalnresult_over90.gff > merge.gff3
	./homology_highppv_merge.sh merge.gff3 longorf.gff3

変更
-2018 7/17　spalnのバグのような出力結果をフィルターするように修正しました
-2018 12/4  篠田さんannotation pipelineに使えるように修正しました
