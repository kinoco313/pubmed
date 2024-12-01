import json
from logging import getLogger, StreamHandler, INFO
import xml.etree.ElementTree as ET

from Bio import Entrez
from openai import OpenAI
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
    pmids = search_pmids(keyword, max_results=20)

    # 論文情報をまとめて取得
    articles = fetch_articles(pmids)
    
    # エクセルファイルとして保存
    save_evidence_table(articles)


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
    logger.info(f"{len(record["IdList"])}件のPMIDを取得しました!")
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
    logger.info("XMLデータの取得に成功しました!")

    articles = []
    for article in root.findall(".//PubmedArticle"):
        pmid = article.findtext(".//PMID")
        title = article.findtext(".//ArticleTitle")
        journal = article.findtext(".//Journal/Title")
        pub_year = article.findtext(".//PubDate/Year")
        summary_dic = {"PMID": pmid, "Title": title, "Journal": journal, "PubYear": pub_year}
        abstract = "".join([abst.text for abst in article.findall(".//AbstractText")])
        
        # 辞書形式でさらに要約されたアブストを結合
        logger.info(f"{title[:10]}...を要約しています・・・")    
        summary_dic |= summarise_abst(abstract)
        articles.append(summary_dic)
            
    return articles

def summarise_abst(abst: str) -> dict[str, str]:
    """アブストから目的や結果などを要約する関数

    Args:
        abst (str): 英文アブスト

    Returns:
        dict[str, str]: キーに項目名、バリューに内容
    """
    
    client = OpenAI()
    prompt = """
    あなたは{# 役割}です。{# 形式}にしたがって{# 入力文}を要約してください。

    # 役割
    英語と日本語が得意な臨床研究の専門家

    # 形式
    ・JSON形式
    ・キーは絶対に「"Objective", "Participants", "Exposure", "Comparison", "Outcome", "Results"」
    ・バリューは各キーに該当する内容を{# 入力文}から抽出し、日本語で100字以内に要約する
    
    # 入力文
    """
    prompt += abst
    

    completion = client.chat.completions.create(
        model="gpt-4o-mini",
        messages=[
            {"role": "system", "content": "You are a helpful assistant."},
            {
                "role": "user",
                "content": prompt
            }
        ]
    )
    # 先頭の```jsonと最後の```を取り除く
    res = completion.choices[0].message.content
    json_str = res[7:-3]
    return json.loads(json_str)

def save_evidence_table(articles: list[dict]) -> None:
    wb = Workbook()
    ws = wb.active
    header = ["PMID", "Title", "Journal", "PubYear", "Objective", "Participants", "Exposure", "Comparison", "Outcome", "Results"]
    ws.append(header)

    # データをExcelに書き込む
    for article in articles:
        tmp = []
        for _, value in article.items():
            tmp.append(value)
        ws.append(tmp)

    # Excelの書式設定
    alignment = Alignment(horizontal="left", vertical="center", wrap_text=True)
    for row in ws.iter_rows(min_row=1, max_row=100, max_col=len(header)):
        for cell in row:
            cell.alignment = alignment

    column_width = {"A": 10, "B": 50, "C": 50, "D": 10, "E": 100, "F": 100, "G":100, "H":100, "I":100, 'J':100}
    for col, width in column_width.items():
        ws.column_dimensions[col].width = width

    file_path = "evidence_table.xlsx"
    try:
        wb.save(file_path)
        logger.info("ファイルの保存が完了しました!")
    except PermissionError:
        logger.exception("ファイルを閉じてください")


if __name__ == "__main__":
    main()
