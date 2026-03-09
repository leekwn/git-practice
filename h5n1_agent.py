#!/usr/bin/env python3
"""
H5N1 Research Agent — fetch latest avian influenza papers and produce structured summaries

Usage:
  pip install openai httpx pydantic python-dateutil pdfminer.six tenacity
  export OPENAI_API_KEY=sk-...
  python h5n1_agent.py --query "(H5N1 OR \"avian influenza\")" --days 45 --limit 10 --out report.md

Notes:
- Sources: Europe PMC (journals + preprints incl. bioRxiv/medRxiv) and Crossref.
- Tries to capture abstracts and open-access full text/PDF when available.
- Summaries are structured (key findings, methods, limitations, impact) via OpenAI Responses API.
"""
from __future__ import annotations
import argparse
import base64
import dataclasses
import datetime as dt
import io
import json
import math
import os
import re
import sys
from dataclasses import dataclass
from typing import List, Optional, Dict, Any, Tuple

import httpx
from dateutil import parser as dateparser
from tenacity import retry, stop_after_attempt, wait_exponential
from pydantic import BaseModel, Field

try:
    from openai import OpenAI
except Exception:
    OpenAI = None

# ----------------------------- Data models -----------------------------

class Paper(BaseModel):
    id: str
    title: str
    doi: Optional[str] = None
    source: str
    url: str
    pdf_url: Optional[str] = None
    abstract: Optional[str] = None
    authors: List[str] = []
    journal: Optional[str] = None
    pub_date: Optional[str] = None  # ISO string
    oa: bool = False

class Summary(BaseModel):
    title: str
    doi: Optional[str]
    source: str
    url: str
    pub_date: Optional[str]
    study_type: str
    population_or_system: str
    geography: str
    methods: str
    sample_size: Optional[str]
    key_results: List[str]
    biosafety_or_public_health_impact: List[str]
    implications_for_policy_or_practice: List[str]
    limitations: List[str]
    evidence_level: str  # e.g., preprint/observational/RCT/systematic review
    confidence: str  # model’s self-rated confidence + rationale

# ----------------------------- Utilities -----------------------------

EP_API = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
CR_API = "https://api.crossref.org/works"
HEADERS = {"User-Agent": "H5N1-Research-Agent/1.0 (mailto:you@example.com)"}

@retry(stop=stop_after_attempt(3), wait=wait_exponential(multiplier=1, min=1, max=8))
async def fetch_json(client: httpx.AsyncClient, url: str, params: Dict[str, Any]) -> Dict[str, Any]:
    r = await client.get(url, params=params, headers=HEADERS, timeout=30)
    r.raise_for_status()
    return r.json()

@retry(stop=stop_after_attempt(3), wait=wait_exponential(multiplier=1, min=1, max=8))
async def fetch_bytes(client: httpx.AsyncClient, url: str) -> bytes:
    r = await client.get(url, headers=HEADERS, timeout=60)
    r.raise_for_status()
    return r.content

async def search_europe_pmc(query: str, since_date: str, limit: int = 20) -> List[Paper]:
    """Search Europe PMC for latest items. since_date in YYYY-MM-DD."""
    q = f"({query}) AND FIRST_PDATE:[{since_date} TO 3000-01-01]"
    params = {
        "query": q,
        "format": "json",
        "pageSize": str(min(100, limit)),
        "sort": "DATE_D",
    }
    async with httpx.AsyncClient() as client:
        data = await fetch_json(client, EP_API, params)
    hits = data.get("resultList", {}).get("result", [])
    papers: List[Paper] = []
    for it in hits[:limit]:
        doi = it.get("doi")
        url = it.get("fullTextUrlList", {}).get("fullTextUrl", [])
        pdf_url = None
        landing = it.get("pubUrl") or it.get("source")
        for u in url:
            if u.get("documentStyle") == "pdf":
                pdf_url = u.get("url")
                break
        papers.append(
            Paper(
                id=f"EPMC:{it.get('id')}",
                title=it.get("title", "").strip(),
                doi=doi,
                source="Europe PMC",
                url=it.get("pubUrl") or it.get("pmcid") and f"https://europepmc.org/abstract/MED/{it.get('pmid')}" or it.get("doi") and f"https://doi.org/{doi}" or landing or "",
                pdf_url=pdf_url,
                abstract=it.get("abstractText"),
                authors=[a.get("fullName") for a in it.get("authorList", {}).get("author", []) if a.get("fullName")],
                journal=it.get("journalTitle"),
                pub_date=it.get("firstPublicationDate"),
                oa=bool(pdf_url),
            )
        )
    return papers

async def search_crossref(query: str, since_date: str, limit: int = 20) -> List[Paper]:
    params = {
        "query": query,
        "filter": f"from-pub-date:{since_date}",
        "sort": "issued",
        "order": "desc",
        "rows": str(limit),
    }
    async with httpx.AsyncClient() as client:
        data = await fetch_json(client, CR_API, params)
    items = data.get("message", {}).get("items", [])
    papers: List[Paper] = []
    for it in items:
        doi = it.get("DOI")
        link = None
        for l in it.get("link", []) or []:
            if l.get("content-type") == "application/pdf":
                link = l.get("URL")
                break
        date_parts = it.get("issued", {}).get("'date-parts", it.get("issued", {}).get("date-parts")) or [[None]]
        ymd = date_parts[0]
        pub_date = None
        if ymd and any(ymd):
            parts = [str(p) for p in ymd if p is not None]
            pub_date = "-".join(parts)
        papers.append(
            Paper(
                id=f"CR:{doi}",
                title=(it.get("title") or [""])[0].strip(),
                doi=doi,
                source="Crossref",
                url=f"https://doi.org/{doi}" if doi else (it.get("URL") or ""),
                pdf_url=link,
                abstract=(it.get("abstract") or "").strip("<>"),
                authors=[" ".join([a.get("given", ""), a.get("family", "")]).strip() for a in it.get("author", [])],
                journal=(it.get("container-title") or [None])[0],
                pub_date=pub_date,
                oa=bool(link),
            )
        )
    return papers

# ----------------------------- PDF text extraction -----------------------------

def extract_text_from_pdf(pdf_bytes: bytes) -> str:
    """Best-effort text extraction using pdfminer.six."""
    from pdfminer.high_level import extract_text
    with io.BytesIO(pdf_bytes) as f:
        try:
            return extract_text(f)
        except Exception:
            return ""

# ----------------------------- LLM summarization -----------------------------

SUMMARY_SCHEMA: Dict[str, Any] = {
    "name": "paper_summary",
    "schema": {
        "type": "object",
        "properties": {
            "study_type": {"type": "string"},
            "population_or_system": {"type": "string"},
            "geography": {"type": "string"},
            "methods": {"type": "string"},
            "sample_size": {"type": "string"},
            "key_results": {"type": "array", "items": {"type": "string"}},
            "biosafety_or_public_health_impact": {"type": "array", "items": {"type": "string"}},
            "implications_for_policy_or_practice": {"type": "array", "items": {"type": "string"}},
            "limitations": {"type": "array", "items": {"type": "string"}},
            "evidence_level": {"type": "string"},
            "confidence": {"type": "string"},
        },
        "required": [
            "study_type","methods","key_results","limitations","evidence_level","biosafety_or_public_health_impact","implications_for_policy_or_practice"
        ],
        "additionalProperties": False,
    },
    "strict": True,
}

SYSTEM_PROMPT = (
    "You are an expert epidemiology analyst. Summarize each paper in crisp, factual bullet points. "
    "Cite only what is in the provided text. If claims are uncertain, say so. Output must follow the provided JSON schema."
)


def _truncate(text: str, max_chars: int = 30000) -> str:
    if not text:
        return ""
    if len(text) <= max_chars:
        return text
    return text[:max_chars]


def summarize_with_openai(client: OpenAI, paper: Paper, text: str) -> Summary:
    content = [
        {
            "role": "user",
            "content": [
                {"type": "input_text", "text": f"TITLE: {paper.title}\nDOI: {paper.doi}\nJOURNAL: {paper.journal}\nDATE: {paper.pub_date}\nURL: {paper.url}\nABSTRACT/FULLTEXT:\n{_truncate(text)}"}
            ],
        }
    ]

    resp = client.responses.create(
        model=os.getenv("OPENAI_SUMMARY_MODEL", "gpt-4o-mini"),
        input=content,
        instructions=SYSTEM_PROMPT,
        response_format={"type": "json_schema", "json_schema": SUMMARY_SCHEMA},
        temperature=0.2,
    )
    # Prefer parsed JSON if available; otherwise parse text
    try:
        payload = json.loads(resp.output_text)
    except Exception:
        # Fallback: attempt to find JSON block
        m = re.search(r"\{[\s\S]*\}$", resp.output_text)
        payload = json.loads(m.group(0)) if m else {}

    summary = Summary(
        title=paper.title,
        doi=paper.doi,
        source=paper.source,
        url=paper.url,
        pub_date=paper.pub_date,
        **payload,
    )
    return summary

# ----------------------------- Orchestration -----------------------------

async def harvest(query: str, days: int, limit: int) -> List[Paper]:
    since = (dt.date.today() - dt.timedelta(days=days)).isoformat()
    epm, cr = await search_europe_pmc(query, since, limit), await search_crossref(query, since, limit)
    # Deduplicate by DOI or title
    seen: Dict[str, Paper] = {}
    for p in epm + cr:
        key = (p.doi or p.title).lower()
        if key not in seen:
            seen[key] = p
        else:
            # prefer the one that has a PDF_url
            if not seen[key].pdf_url and p.pdf_url:
                seen[key] = p
    return list(seen.values())[:limit]

async def enrich_with_fulltext(papers: List[Paper]) -> List[Tuple[Paper, str]]:
    out: List[Tuple[Paper, str]] = []
    async with httpx.AsyncClient() as client:
        for p in papers:
            text = p.abstract or ""
            if p.pdf_url:
                try:
                    pdf = await fetch_bytes(client, p.pdf_url)
                    t = extract_text_from_pdf(pdf)
                    if len(t) > len(text):
                        text = t
                except Exception:
                    pass
            out.append((p, text))
    return out

async def run_agent(query: str, days: int, limit: int, out_path: str) -> str:
    papers = await harvest(query, days, limit)
    enriched = await enrich_with_fulltext(papers)

    if OpenAI is None:
        raise RuntimeError("openai Python package not installed. pip install openai")
    client = OpenAI()

    summaries: List[Summary] = []
    for p, text in enriched:
        try:
            summaries.append(summarize_with_openai(client, p, text or p.abstract or ""))
        except Exception as e:
            # fallback minimal summary
            summaries.append(
                Summary(
                    title=p.title,
                    doi=p.doi,
                    source=p.source,
                    url=p.url,
                    pub_date=p.pub_date,
                    study_type="unknown",
                    population_or_system="unknown",
                    geography="unknown",
                    methods="unknown",
                    sample_size=None,
                    key_results=["See abstract/full text at source."],
                    biosafety_or_public_health_impact=[],
                    implications_for_policy_or_practice=[],
                    limitations=[f"Summarization failed: {e}"],
                    evidence_level="unknown",
                    confidence="low: fallback summary",
                )
            )

    # Compose a markdown report
    lines = []
    today = dt.datetime.now().strftime("%Y-%m-%d")
    lines.append(f"# Avian Influenza (H5N1) — Latest Research Digest (as of {today})\n")
    lines.append(f"Query: `{query}` | Window: last {days} days | Items: {len(summaries)}\n")

    for i, s in enumerate(summaries, 1):
        lines.append(f"## {i}. {s.title}")
        meta = []
        if s.doi: meta.append(f"DOI: [{s.doi}](https://doi.org/{s.doi})")
        if s.pub_date: meta.append(f"Date: {s.pub_date}")
        if s.url: meta.append(f"Link: {s.url}")
        if meta:
            lines.append("**" + " | ".join(meta) + "**\n")
        lines.append(f"- **Study type:** {s.study_type}")
        if s.population_or_system:
            lines.append(f"- **Population/System:** {s.population_or_system}")
        if s.geography:
            lines.append(f"- **Geography:** {s.geography}")
        if s.sample_size:
            lines.append(f"- **Sample size:** {s.sample_size}")
        lines.append(f"- **Methods:** {s.methods}")
        lines.append("- **Key results:**\n" + "\n".join([f"  - {k}" for k in s.key_results]))
        if s.biosafety_or_public_health_impact:
            lines.append("- **Public health impact:**\n" + "\n".join([f"  - {k}" for k in s.biosafety_or_public_health_impact]))
        if s.implications_for_policy_or_practice:
            lines.append("- **Implications:**\n" + "\n".join([f"  - {k}" for k in s.implications_for_policy_or_practice]))
        if s.limitations:
            lines.append("- **Limitations:**\n" + "\n".join([f"  - {k}" for k in s.limitations]))
        lines.append(f"- **Evidence level:** {s.evidence_level}")
        lines.append(f"- **Confidence:** {s.confidence}\n")

    report = "\n".join(lines)
    with open(out_path, "w", encoding="utf-8") as f:
        f.write(report)
    return out_path

# ----------------------------- CLI -----------------------------

def main():
    ap = argparse.ArgumentParser(description="H5N1 Research Agent")
    ap.add_argument("--query", default="(H5N1 OR \"avian influenza\")")
    ap.add_argument("--days", type=int, default=45)
    ap.add_argument("--limit", type=int, default=8)
    ap.add_argument("--out", default="h5n1_digest.md")
    args = ap.parse_args()

    import asyncio
    out_path = asyncio.run(run_agent(args.query, args.days, args.limit, args.out))
    print(f"Saved report -> {out_path}")

if __name__ == "__main__":
    main()
