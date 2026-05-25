# activity.py
import csv
import logging
import json
from src.core import common
from src.core.config import Config
from src.core.exceptions import APIError, FileProcessingError

logger = logging.getLogger("dbaasp_pipeline")


def get_activities(peptide_json):
    a = peptide_json.get("targetActivities")
    if a is None:
        a = peptide_json.get("activityAgainstTargetSpecies")
    return a if isinstance(a, list) else []


def get_unusual_amino_acids_map(peptide_json):
    unusual = peptide_json.get("unusualAminoAcids", [])
    if not unusual:
        return "{}"
    mapping = {}
    for aa in unusual:
        pos = aa.get("position")
        name = aa.get("modificationType", {}).get("name")
        if pos is not None and name:
            mapping[str(pos)] = name
    return json.dumps(mapping)


def get_unusual_amino_acids(peptide_json):
    unusual = peptide_json.get("unusualAminoAcids", [])
    if not unusual:
        return ""
    names = [aa.get("modificationType", {}).get("name", "") for aa in unusual if aa.get("modificationType")]
    return ", ".join(filter(None, names))


def format_article_citation(article: dict) -> str:
    """Format a single article dict into a readable citation string.

    Format: "Author1, Author2, ... | Title | Journal Year;Vol:Pages | PMID:xxx"
    """
    authors = article.get("authors", [])
    author_names = ", ".join(a.get("name", "") for a in authors if a.get("name"))

    title = (article.get("title") or "").strip()

    journal = (article.get("journal") or {}).get("name", "")
    year = str(article.get("year", ""))
    volume = str(article.get("volume", ""))
    pages = str(article.get("pages", ""))

    # Build journal reference part: "J Pept Sci 2017;23:45-55"
    journal_ref = journal
    if year:
        journal_ref += f" {year}"
    if volume:
        journal_ref += f";{volume}"
    if pages:
        journal_ref += f":{pages}"

    pubmed_id = (article.get("pubmed") or {}).get("pubmedId", "")
    pmid_part = f"PMID:{pubmed_id}" if pubmed_id else ""

    parts = [p for p in [author_names, title, journal_ref, pmid_part] if p]
    return " | ".join(parts)


def get_reference_citation(peptide_json: dict, ref_id: str) -> str:
    """Resolve a reference ID (1-based string index) to a formatted citation.

    The DBAASP API stores a short integer string (e.g. "1") in each activity
    row's `reference` field. This is a 1-based index into the peptide's
    top-level `articles` list. Falls back to all articles (joined by " ;; ")
    if the index is missing or out of range.
    """
    articles = peptide_json.get("articles", [])
    if not articles:
        return ""

    # Try to resolve by 1-based index
    if ref_id:
        try:
            idx = int(ref_id) - 1
            if 0 <= idx < len(articles):
                return format_article_citation(articles[idx])
        except (ValueError, TypeError):
            pass

    # Fallback: format all articles
    return " ;; ".join(format_article_citation(a) for a in articles)


def collect_activity_keys(all_peptides):
    exclude = {"id", "activity", "activityMeasureValue"}
    order, seen = [], set()
    for d in all_peptides:
        for a in get_activities(d):
            for k in a.keys():
                if not k or k.lower() in exclude:
                    continue
                if k not in seen:
                    seen.add(k)
                    order.append(k)
    return order


def flatten_value(v):
    if v is None:
        return ""
    if isinstance(v, dict):
        return v.get("name", str(v))
    return v


def run():
    logger.info("Starting activity data collection")

    try:
        ids = common.load_ids()
        if not ids:
            logger.warning(f"No IDs found in {Config.INPUT_PEPTIDES_CSV}")
            return

        logger.info(f"Processing {len(ids)} peptides for activity data")
        peptides = []
        failed_peptides = []

        for i, pid in enumerate(ids, start=1):
            print(f"[{i}/{len(ids)}] Fetching activity for peptide {pid} ...", end="", flush=True)
            try:
                peptides.append(common.fetch(pid))
                print(" ok")
            except APIError as e:
                logger.error(f"API error for peptide {pid}: {e}")
                failed_peptides.append(pid)
                print(f" fail (API error)")
            except Exception as e:
                logger.error(f"Unexpected error for peptide {pid}: {e}")
                failed_peptides.append(pid)
                print(f" fail (unexpected error)")

        if failed_peptides:
            logger.warning(f"Failed to fetch {len(failed_peptides)} peptides: {failed_peptides}")

        if not peptides:
            logger.error("No data fetched; nothing to write.")
            return

        logger.info(f"Successfully fetched data for {len(peptides)} peptides")
        activity_cols = collect_activity_keys(peptides)
        header = (
            ["Peptide ID", "N TERMINUS", "SEQUENCE", "C TERMINUS",
             "Unusual Amino Acids", "Unusual Amino Acids Map"]
            + activity_cols
            + ["reference_citation"]
        )

        try:
            with open(Config.OUTPUT_ACTIVITY_CSV, "w", newline="", encoding=Config.CSV_ENCODING) as f:
                w = csv.DictWriter(f, fieldnames=header)
                w.writeheader()

                for d in peptides:
                    base = {
                        "Peptide ID": str(d.get("id", "")),
                        "N TERMINUS": (d.get("nTerminus") or {}).get("name", ""),
                        "SEQUENCE": d.get("sequence", ""),
                        "C TERMINUS": (d.get("cTerminus") or {}).get("name", ""),
                        "Unusual Amino Acids": get_unusual_amino_acids(d),
                        "Unusual Amino Acids Map": get_unusual_amino_acids_map(d),
                    }
                    acts = get_activities(d)
                    if not acts:
                        w.writerow(base)
                        continue
                    for a in acts:
                        row = dict(base)
                        for k in activity_cols:
                            row[k] = flatten_value(a.get(k))
                        # Resolve full citation from articles list using 1-based reference index
                        row["reference_citation"] = get_reference_citation(d, str(a.get("reference", "")))
                        w.writerow(row)

            logger.info(f"Successfully wrote activity data to {Config.OUTPUT_ACTIVITY_CSV}")

        except IOError as e:
            raise FileProcessingError(f"Error writing activity file: {e}", filename=Config.OUTPUT_ACTIVITY_CSV)

    except (FileProcessingError, APIError) as e:
        logger.error(f"Activity collection failed: {e}")
        raise
    except Exception as e:
        logger.error(f"Unexpected error in activity collection: {e}")
        raise
