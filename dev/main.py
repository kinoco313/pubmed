from Bio import Entrez
import xml.etree.ElementTree as ET
import deepl
from datetime import date


Entrez.email = "example@gmail.com"

def main():
    keyword = '((New England Journal of Medicine[Journal] OR BMJ[Journal] OR The Lancet[Journal] OR JAMA[Journal] OR Annals of Internal Medicine[Journal] OR Kidney International[Journal] OR Journal of the American Society of Nephrology[Journal] OR American Journal of Kidney Diseases[Journal] OR Clinical Journal of the American Society of Nephrology[Journal] OR Nephrology Dialysis Transplantation[Journal])) AND (((CVD) AND (egfr) AND ((slope) OR (change) OR (decline) OR (increase)) ) NOT (surrogate))'
    pmids = search_pmids(keyword,max_results=10000)
    root = fetch_article(pmid=pmids[0])
    articule_title = extract_article_title(root)
    
    print(articule_title)
    

def search_pmids(keyword: str, max_results: int=10) -> list[str]:
    """キーワードに関連する論文のPMIDを検索する関数

    Args:
        keyword (str): 検索したいキーワード
        max_results (int, optional): 取得したい件数の最大値 デフォルトは10件。最大で10000件まで可能

    Returns:
        list[str]: 長さが最大でmax_resultsのリスト
    """
    stream = Entrez.esearch(db="pubmed", term=keyword, retmax=str(max_results))
    record = Entrez.read(stream)
    return record['IdList']


def fetch_article(pmid: str) -> ET.Element:
    """指定されたPMIDからXML形式のデータを取得する関数

    Args:
        pmid (str): アブストを取得したい論文のPMID

    Returns:
        ET.Element: XML形式データの根要素
    """
    handle = Entrez.efetch(db="pubmed", id=pmid, retmode="xml")
    xml = handle.read()
    # xmlの構造を見てみたい
    # with open('test_result.xml', 'wb') as f:
    #     f.write(xml)
    
    root = ET.fromstring(xml)
    
    return root


def extract_article_title(root: ET.Element) -> str:
    """論文のタイトルタグのテキストを抽出する関数

    Args:
        root (ET.Element): XML形式データの根要素

    Returns:
        str: 論文のタイトル
    """
    
    ArticleTitle = (
        root
        .find("PubmedArticle")
        .find("MedlineCitation")
        .find("Article")
        .find("ArticleTitle")
        .text
    )
    
    return ArticleTitle
    

def extract_abst(root: ET.Element) -> str:
    """XMLデータからアブストタグのテキストを抽出する関数

    Args:
        root (ET.Element): XML形式データの根要素

    Returns:
        str: 論文の英文アブスト
    """
    
    AbstractText = (
        root
        .find("PubmedArticle")
        .find("MedlineCitation")
        .find("Article")
        .find("Abstract")
        .find("AbstractText")
        .text
    )
        
    return AbstractText

def translate_to_ja(abst: str) -> str:
    """アブストを日本語に要約する関数

    Args:
        abst (str): アブスト（英文）

    Returns:
        str: アブスト（日本語訳）
    """
    auth_key = "your-api-key"
    translator = deepl.Translator(auth_key)
    result = translator.translate_text(abst, target_lang="JA")
    return result.text




if __name__ == "__main__":
    main()
    
    