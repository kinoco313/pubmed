from Bio import Entrez
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor, as_completed
import deepl
from openpyxl import Workbook
from openpyxl.styles import Alignment

Entrez.email = "example@gmail.com"

def main():
    # キーワードからPMIDを取得
    keyword = '((New England Journal of Medicine[Journal] OR BMJ[Journal] OR The Lancet[Journal] OR JAMA[Journal] OR Annals of Internal Medicine[Journal] OR Kidney International[Journal] OR Journal of the American Society of Nephrology[Journal] OR American Journal of Kidney Diseases[Journal] OR Clinical Journal of the American Society of Nephrology[Journal] OR Nephrology Dialysis Transplantation[Journal])) AND (((CVD) AND (egfr) AND ((slope) OR (change) OR (decline) OR (increase)) ) NOT (surrogate))'
    max_results = 30
    pmids = search_pmids(keyword, max_results)

    # 論文情報をまとめて取得
    articles = fetch_articles(pmids)

    # 翻訳を並列処理で実施
    abstracts_ja = translate_abstracts_parallel([article["abstract"] for article in articles])

    # ワークブックを作成
    wb = Workbook()
    ws = wb.active
    ws.append(["PMID", "Title", "Journal", "PubYear", "Abst_en", "Abst_ja"])

    # データをExcelに書き込む
    for article, abst_ja in zip(articles, abstracts_ja):
        ws.append([article["pmid"], article["title"], article["journal"], article["pub_year"], article["abstract"], abst_ja])

    # Excelの書式設定
    alignment = Alignment(horizontal="left", vertical="center", wrap_text=True)
    for row in ws.iter_rows(min_row=1, max_row=max_results+1, max_col=6):
        for cell in row:
            cell.alignment = alignment

    column_width = {"A": 10, "B": 50, "C": 50, "D": 10, "E": 100, "F": 100}
    for col, width in column_width.items():
        ws.column_dimensions[col].width = width

    wb.save("test.xlsx")

def search_pmids(keyword: str, max_results: int = 10) -> list[str]:
    stream = Entrez.esearch(db="pubmed", term=keyword, retmax=str(max_results))
    record = Entrez.read(stream)
    return record['IdList']

def fetch_articles(pmids: list[str]) -> list[dict]:
    """PMIDリストからまとめて論文情報を取得"""
    handle = Entrez.efetch(db="pubmed", id=",".join(pmids), retmode="xml")
    xml = handle.read()
    root = ET.fromstring(xml)

    articles = []
    for article in root.findall(".//PubmedArticle"):
        try:
            title = article.findtext(".//ArticleTitle", default="")
            journal = article.findtext(".//Journal/Title", default="")
            pub_year = article.findtext(".//PubDate/Year", default="")
            abstract = "".join([abst.text for abst in article.findall(".//AbstractText") if abst.text])
            pmid = article.findtext(".//PMID", default="")
            articles.append({
                "pmid": pmid,
                "title": title,
                "journal": journal,
                "pub_year": pub_year,
                "abstract": abstract
            })
        except Exception as e:
            print(f"Error processing article: {e}")
    return articles

def translate_abstracts_parallel(abstracts: list[str], max_workers: int = 5) -> list[str]:
    """アブストラクトを並列で翻訳"""
    with open("apikey.txt") as f:
        auth_key = f.readline().strip()
    translator = deepl.Translator(auth_key)

    def translate(abst: str) -> str:
        return translator.translate_text(abst, target_lang="JA").text if abst else ""

    results = []
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_abst = {executor.submit(translate, abst): abst for abst in abstracts}
        for future in as_completed(future_to_abst):
            try:
                results.append(future.result())
            except Exception as e:
                print(f"Translation error: {e}")
                results.append("")

    return results

if __name__ == "__main__":
    main()
