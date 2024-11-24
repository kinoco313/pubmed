from Bio import Entrez

Entrez.email = "n23m620c@mail.cc.niigata-u.ac.jp"

def main():
    keyword = '((CVD) AND (egfr) AND ((slope) OR (change) OR (decline) OR (increase)) ) NOT (surrogate)'
    record = get_article_pmids(keyword,max_results=10000)

def get_article_pmids(keyword: str, max_results: int=10) -> list[str]:
    """キーワードに関連する論文のPMIDを返す関数

    Args:
        keyword (str): 検索したいキーワード
        max_results (int, optional): 取得したい件数の最大値 デフォルトは10件。最大で10000件まで可能

    Returns:
        list[str]: 長さが最大でmax_resules件のリスト
    """
    stream = Entrez.esearch(db="pubmed",
                            term=keyword,
                            retmax=str(max_results))
    record = Entrez.read(stream)
    return record['IdList']



if __name__ == "__main__":
    main()
    
    