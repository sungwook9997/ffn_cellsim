"""External-literature adapter for Aleph's public RAG--TAG--CAG harness."""

from .ingest import IngestResult, ingest_corpus, verify_idempotent_ingest

__all__ = ["IngestResult", "ingest_corpus", "verify_idempotent_ingest"]
