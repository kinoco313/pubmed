from logging import getLogger, StreamHandler, INFO
import xml.etree.ElementTree as ET

from Bio import Entrez
import deepl
from openpyxl import Workbook
from openpyxl.styles import Alignment

# 検索用の設定
Entrez.email = "example@gmail.com"
keyword = "(New England Journal of Medicine[Journal] OR BMJ[Journal] OR The Lancet[Journal] OR JAMA[Journal] OR Annals of Internal Medicine[Journal] OR Kidney International[Journal] OR Journal of the American Society of Nephrology[Journal] OR American Journal of Kidney Diseases[Journal] OR Clinical Journal of the American Society of Nephrology[Journal] OR Nephrology Dialysis Transplantation[Journal]) AND cvd AND egfr AND (slope OR change)"

# ログ出力の設定
logger = getLogger(__name__)
logger.setLevel(INFO)
s_handler = StreamHandler()
s_handler.setLevel(INFO)
logger.addHandler(s_handler)

def main():
    # キーワードからPMIDを取得
    pmids = search_pmids(keyword, max_results=10)

    # 論文情報をまとめて取得
    articles = fetch_articles(pmids)

    # # ワークブックを作成
    # wb = Workbook()
    # ws = wb.active
    # header = ["PMID", "Title", "Journal", "PubYear", "Abst_en"]
    # ws.append(header)

    # # データをExcelに書き込む
    # for article in articles:
    #     ws.append([article["pmid"], article["title"], article["journal"], article["pub_year"], article["abstract"]])

    # # Excelの書式設定
    # alignment = Alignment(horizontal="left", vertical="center", wrap_text=True)
    # for row in ws.iter_rows(min_row=1, max_row=100, max_col=len(header)):
    #     for cell in row:
    #         cell.alignment = alignment

    # column_width = {"A": 10, "B": 50, "C": 50, "D": 10, "E": 100, "F": 100}
    # for col, width in column_width.items():
    #     ws.column_dimensions[col].width = width

    # file_path = "evidence_table.xlsx"
    # wb.save(file_path)
    

def search_pmids(keyword: str, max_results: int = 10) -> list[str]:
    """キーワードに合致する論文のPMIDをリストとして返す関数

    Args:
        keyword (str): 検索ワード できれば検索式が望ましい
        max_results (int, optional): 探してくる論文の最大件数 デフォルトは10 10~30ぐらいがちょうどいい

    Returns:
        list[str]: PMIDのリスト。
    """
    stream = Entrez.esearch(db="pubmed", term=keyword, retmax=str(max_results))
    record = Entrez.read(stream)
    # 取得できているか確認
    logger.info(f"{len(record["IdList"])}件のPMIDを取得")
    return record['IdList']

def fetch_articles(pmids: list[str]) -> list[dict]:
    """pmidリストから論文情報を取得し、辞書のリストとして返す関数

    Args:
        pmids (list[str]): PMIDのリスト

    Returns:
        list[dict]: 論文情報を辞書として格納したリスト
    """
    handle = Entrez.efetch(db="pubmed", id=",".join(pmids), retmode="xml")
    xml = handle.read()
    root = ET.fromstring(xml)
    logger.info("XMLデータの取得に成功")

    articles = []
    for article in root.findall(".//PubmedArticle"):
        title = article.findtext(".//ArticleTitle")
        journal = article.findtext(".//Journal/Title")
        pub_year = article.findtext(".//PubDate/Year")
        abstract = "".join([abst.text for abst in article.findall(".//AbstractText")])
        pmid = article.findtext(".//PMID", default="")
        articles.append({
            "pmid": pmid,
            "title": title,
            "journal": journal,
            "pub_year": pub_year,
            "abstract": abstract
        })
    logger.info("論文情報の取得を完了")
        
    return articles



if __name__ == "__main__":
    main()
