from Bio import Entrez
import xml.etree.ElementTree as ET
import deepl
from openpyxl import Workbook
from openpyxl.styles import Alignment

Entrez.email = "example@gmail.com"

def main():
    # キーワードからPMIDを取得
    keyword = '((New England Journal of Medicine[Journal] OR BMJ[Journal] OR The Lancet[Journal] OR JAMA[Journal] OR Annals of Internal Medicine[Journal] OR Kidney International[Journal] OR Journal of the American Society of Nephrology[Journal] OR American Journal of Kidney Diseases[Journal] OR Clinical Journal of the American Society of Nephrology[Journal] OR Nephrology Dialysis Transplantation[Journal])) AND (((CVD) AND (egfr) AND ((slope) OR (change) OR (decline) OR (increase)) ) NOT (surrogate))'
    max_results: int = 10 # ワークシートの行数に使いたいので変数にしておく
    pmids = search_pmids(keyword,max_results)
    
    # ワークブックを開き、ワークシートをアクティベートする
    wb = Workbook()
    ws = wb.active
    
    # 1行目に項目名を記述する
    ws.append(["PMID", "Title", "Jounal", "PubYear", "Abst_en"]) # 後でAbst_jaも追加予定
    
    # forループで各情報を追加していく
    for pmid in pmids:
        
        # 情報源(XML)
        root = fetch_article(pmid)
        
        # 各情報
        articule_title = extract_article_title(root)
        jounal_title = extract_journal_title(root)
        pub_year = extract_pubdate(root)
        abst_en = extract_abst(root)
        
        # evid_tableに追加
        ws.append([pmid, articule_title, jounal_title, pub_year, abst_en])
        
    # 上下中央揃え、左揃え、折り返して全体を表示
    alignment = Alignment(horizontal="left", vertical="center", wrap_text=True)
    for row in ws.iter_rows(min_row=1, min_col=1, max_row=max_results+1, max_col=6):
        for cell in row:
            cell.alignment = alignment
    
    # 列幅を指定
    column_width = {"A":10, "B":50, "C":50, "D":10, "E":200}
    for column, width in column_width.items():
            ws.column_dimensions[column].width = width
    
    # ワークブックを保存
    wb.save("test.xlsx")
    
    

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
    with open('test_result.xml', 'wb') as f:
        f.write(xml)
    
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

def extract_journal_title(root: ET.Element) -> str:
    """雑誌のタイトルを抽出する関数

    Args:
        root (ET.Element): XML形式データの根要素

    Returns:
        str: 雑誌のタイトル
    """
    
    Title = (
        root
        .find("PubmedArticle")
        .find("MedlineCitation")
        .find("Article")
        .find("Journal")
        .find("Title")
        .text
    )
    
    return Title


def extract_pubdate(root: ET.Element) -> str:
    """論文の出版された年を抽出する関数

    Args:
        root (ET.Element): XML形式データの根要素

    Returns:
        str: 論文が出版された年
    """
    
    PubDate = (
        root
        .find("PubmedArticle")
        .find("MedlineCitation")
        .find("Article")
        .find("Journal")
        .find("JournalIssue")
        .find("PubDate")
    )
    
    Year = PubDate.find("Year").text

    return Year
    

def extract_abst(root: ET.Element) -> str:
    """XMLデータからアブストタグのテキストを抽出する関数

    Args:
        root (ET.Element): XML形式データの根要素

    Returns:
        str: 論文の英文アブスト
    """
    
    AbstractText_list = (
        root
        .find("PubmedArticle")
        .find("MedlineCitation")
        .find("Article")
        .find("Abstract")
        .findall("AbstractText")
    )
    
    AbstractText_str = [abst.text for abst in AbstractText_list]
    abst_text_joined = "".join(AbstractText_str)
        
    return abst_text_joined

def translate_into_ja(abst: str) -> str:
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
    
    