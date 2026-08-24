#!/usr/bin/env python3
"""LLM backends for the TAG synthesis and answer-generation stages.

The Anthropic SDK path supports an explicit cached prefix. The ``claude -p``
fallback cannot express ``cache_control``; it emits a warning and runs as a
tool-free completion from an empty temporary directory so the cache bypass is
visible and the project tree cannot leak into the prompt.
"""
from __future__ import annotations

import json
import logging
import os
import shutil
import subprocess
import tempfile
import time
from dataclasses import dataclass
from typing import Any

DEFAULT_MODEL = "claude-opus-5"
DEFAULT_CACHE_TTL = "1h"
SUPPORTED_CACHE_TTLS = frozenset({"5m", "1h"})

LOGGER = logging.getLogger(__name__)


@dataclass(frozen=True)
class Completion:
    """One backend completion plus provider-reported accounting."""

    text: str
    backend: str
    model: str
    elapsed_s: float
    usage: dict[str, Any]


def _usage_dict(usage: Any) -> dict[str, Any]:
    """Convert an Anthropic SDK usage object to plain serialisable data."""
    if usage is None:
        return {}
    if hasattr(usage, "model_dump"):
        return usage.model_dump()
    if isinstance(usage, dict):
        return dict(usage)
    return {
        name: getattr(usage, name)
        for name in (
            "input_tokens",
            "cache_creation_input_tokens",
            "cache_read_input_tokens",
            "output_tokens",
        )
        if hasattr(usage, name)
    }


def _anthropic_completion(
    prompt: str,
    *,
    system: str,
    model: str,
    max_tokens: int,
    cache_prefix: str,
    cache_key: str,
    cache_ttl: str,
    api_key: str,
) -> Completion:
    """Call the Messages API, optionally marking the static user prefix."""
    import anthropic

    content: str | list[dict[str, Any]]
    if cache_prefix:
        if cache_ttl not in SUPPORTED_CACHE_TTLS:
            allowed = ", ".join(sorted(SUPPORTED_CACHE_TTLS))
            raise ValueError(f"cache_ttl must be one of {allowed}; got {cache_ttl!r}")
        content = [
            {
                "type": "text",
                "text": cache_prefix,
                "cache_control": {"type": "ephemeral", "ttl": cache_ttl},
            },
            {"type": "text", "text": prompt},
        ]
    else:
        content = prompt

    kwargs: dict[str, Any] = {
        "model": model,
        "max_tokens": max_tokens,
        "system": system or None,
        "messages": [{"role": "user", "content": content}],
    }
    # Opus 5 rejects non-default temperature values. Disable thinking instead:
    # these bounded SQL/grounded-answer completions were designed around visible
    # output budgets, not an adaptive hidden-thinking budget.
    if model == DEFAULT_MODEL:
        kwargs["thinking"] = {"type": "disabled"}
        kwargs["output_config"] = {"effort": "low"}

    client = anthropic.Anthropic(api_key=api_key)
    started = time.perf_counter()
    msg = client.messages.create(**kwargs)
    elapsed = time.perf_counter() - started
    usage = _usage_dict(msg.usage)
    LOGGER.info(
        "prompt-cache backend=anthropic model=%s key=%s ttl=%s "
        "input=%s create=%s read=%s output=%s elapsed_s=%.3f",
        model,
        cache_key[:12] if cache_key else "-",
        cache_ttl if cache_prefix else "-",
        usage.get("input_tokens", 0),
        usage.get("cache_creation_input_tokens", 0),
        usage.get("cache_read_input_tokens", 0),
        usage.get("output_tokens", 0),
        elapsed,
    )
    return Completion(
        text="".join(block.text for block in msg.content if block.type == "text"),
        backend="anthropic",
        model=model,
        elapsed_s=elapsed,
        usage=usage,
    )


def _claude_cli_completion(
    prompt: str,
    *,
    system: str,
    model: str,
    cache_prefix: str,
    cache_key: str,
) -> Completion:
    """Run the uncached Claude Code fallback as an isolated completion."""
    executable = shutil.which("claude")
    if executable is None:
        raise RuntimeError(
            "ANTHROPIC_API_KEY is unset and the `claude` fallback is not on PATH"
        )

    if cache_prefix:
        LOGGER.warning(
            "prompt-cache BYPASS backend=claude-cli key=%s: `claude -p` "
            "cannot send cache_control; resending %d static characters. "
            "Set ANTHROPIC_API_KEY to enable the SDK cache.",
            cache_key[:12] if cache_key else "-",
            len(system) + len(cache_prefix),
        )
    full_prompt = f"{cache_prefix}\n\n{prompt}" if cache_prefix else prompt
    minimal_system = system or "Follow the user's instructions exactly."
    cmd = [
        executable,
        "-p",
        "--tools",
        "",
        "--system-prompt",
        minimal_system,
        "--model",
        model,
        "--no-session-persistence",
        "--output-format",
        "json",
    ]
    with tempfile.TemporaryDirectory(prefix="tag-kb-claude-") as scratch:
        started = time.perf_counter()
        result = subprocess.run(
            cmd,
            input=full_prompt,
            text=True,
            capture_output=True,
            timeout=180,
            cwd=scratch,
        )
        elapsed = time.perf_counter() - started
    if result.returncode != 0:
        raise RuntimeError(f"claude -p failed: {result.stderr[:400]}")
    try:
        payload = json.loads(result.stdout)
    except json.JSONDecodeError as exc:
        raise RuntimeError(
            f"claude -p returned invalid JSON: {result.stdout[:400]}"
        ) from exc
    usage = payload.get("usage") or {}
    LOGGER.info(
        "completion backend=claude-cli model=%s cache=bypass input=%s "
        "create=%s read=%s output=%s elapsed_s=%.3f",
        model,
        usage.get("input_tokens", 0),
        usage.get("cache_creation_input_tokens", 0),
        usage.get("cache_read_input_tokens", 0),
        usage.get("output_tokens", 0),
        elapsed,
    )
    return Completion(
        text=str(payload.get("result", "")).strip(),
        backend="claude-cli",
        model=model,
        elapsed_s=elapsed,
        usage=usage,
    )


def complete(
    prompt: str,
    *,
    system: str = "",
    model: str | None = None,
    max_tokens: int = 1500,
    cache_prefix: str = "",
    cache_key: str = "",
    cache_ttl: str = DEFAULT_CACHE_TTL,
) -> Completion:
    """Complete one prompt through the SDK or explicit uncached fallback.

    Args:
        prompt: Per-request, non-cacheable prompt suffix.
        system: System instruction. It becomes part of the SDK cache prefix when
            ``cache_prefix`` is present.
        model: Anthropic model id. Defaults to Claude Opus 5.
        max_tokens: Maximum response tokens.
        cache_prefix: Stable content block to end with a cache breakpoint.
        cache_key: Diagnostic key identifying the prefix, normally a schema hash.
        cache_ttl: Anthropic ephemeral cache TTL, ``5m`` or ``1h``.

    Returns:
        Completion text and provider accounting.
    """
    selected_model = model or DEFAULT_MODEL
    api_key = os.environ.get("ANTHROPIC_API_KEY", "").strip()
    if api_key:
        return _anthropic_completion(
            prompt,
            system=system,
            model=selected_model,
            max_tokens=max_tokens,
            cache_prefix=cache_prefix,
            cache_key=cache_key,
            cache_ttl=cache_ttl,
            api_key=api_key,
        )
    return _claude_cli_completion(
        prompt,
        system=system,
        model=selected_model,
        cache_prefix=cache_prefix,
        cache_key=cache_key,
    )


def llm(
    prompt: str,
    system: str = "",
    model: str | None = None,
    max_tokens: int = 1500,
    temperature: float | None = None,
    *,
    cache_prefix: str = "",
    cache_key: str = "",
    cache_ttl: str = DEFAULT_CACHE_TTL,
) -> str:
    """Compatibility wrapper returning completion text only.

    ``temperature`` remains in the public signature for existing callers but is
    intentionally ignored: Opus 5 rejects non-default sampling temperatures.
    """
    if temperature not in (None, 1.0):
        LOGGER.debug(
            "ignoring temperature=%s because current Anthropic models reject "
            "non-default sampling temperatures",
            temperature,
        )
    return complete(
        prompt,
        system=system,
        model=model,
        max_tokens=max_tokens,
        cache_prefix=cache_prefix,
        cache_key=cache_key,
        cache_ttl=cache_ttl,
    ).text
