from Bio import Entrez
import xml.etree.ElementTree as ET

Entrez.email = "example@gmail.com"

def main():
    keyword = '((CVD) AND (egfr) AND ((slope) OR (change) OR (decline) OR (increase)) ) NOT (surrogate)'
    record = get_article_pmids(keyword,max_results=10000)
    result = get_absts(record[0])
    print(result)
    
    # xmlの構造を見てみたい
    # with open('test_result.xml', 'wb') as f:
    #     f.write(result)

def get_article_pmids(keyword: str, max_results: int=10) -> list[str]:
    """キーワードに関連する論文のPMIDを取得する関数

    Args:
        keyword (str): 検索したいキーワード
        max_results (int, optional): 取得したい件数の最大値 デフォルトは10件。最大で10000件まで可能

    Returns:
        list[str]: 長さが最大でmax_resulesのリスト
    """
    stream = Entrez.esearch(db="pubmed", term=keyword, retmax=str(max_results))
    record = Entrez.read(stream)
    return record['IdList']


def get_absts(pmid: str) -> str:
    """論文のアブストを取得する関数

    Args:
        pmid (str): アブストを取得したい論文のPMID

    Returns:
        str: 論文のアブスト
    """
    handle = Entrez.efetch(db="pubmed", id=pmid, retmode="xml")
    xml = handle.read()
    PubmedArticleSet = ET.fromstring(xml)
    AbstractText = (
        PubmedArticleSet
        .find("PubmedArticle")
        .find("MedlineCitation")
        .find("Article")
        .find("Abstract")
        .find("AbstractText")
        .text
    )
        
    return AbstractText


if __name__ == "__main__":
    main()
    
    