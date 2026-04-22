#!/usr/bin/env python3
"""
Unified LLM API Configuration for variant classification.

Provides a consistent interface for LLM API keys across all modules:
- literature_classifier.py (OpenAI GPT-4o)
- rag_llm_evidence_adjuster.py (DeepSeek)

Usage:
    from llm_config import get_llm_config

    # Get API key and provider
    config = get_llm_config()  # uses environment variables
    config = get_llm_config("openai")  # specifically OpenAI
    config = get_llm_config("deepseek")  # specifically DeepSeek

    if config["available"]:
        api_key = config["api_key"]
        provider = config["provider"]  # "openai" or "deepseek"
"""

import os
from typing import Optional
from dataclasses import dataclass
from enum import Enum


class LLMProvider(Enum):
    """Supported LLM providers."""
    OPENAI = "openai"
    DEEPSEEK = "deepseek"


@dataclass
class LLMConfig:
    """LLM configuration."""
    available: bool
    api_key: Optional[str] = None
    provider: Optional[str] = None
    model: Optional[str] = None


def get_openai_api_key() -> Optional[str]:
    """Get OpenAI API key from environment or parameter."""
    return os.environ.get("OPENAI_API_KEY") or os.environ.get("LLM_API_KEY")


def get_deepseek_api_key() -> Optional[str]:
    """Get DeepSeek API key from environment or parameter."""
    return os.environ.get("DEEPSEEK_API_KEY") or os.environ.get("LLM_API_KEY")


def get_llm_config(provider: Optional[str] = None) -> LLMConfig:
    """
    Get unified LLM configuration.

    Priority:
    1. If provider specified, use that provider's API key
    2. If OPENAI_API_KEY set, use OpenAI
    3. If DEEPSEEK_API_KEY set, use DeepSeek
    4. If LLM_API_KEY set, use OpenAI as default
    5. Otherwise, not available

    Args:
        provider: Optional specific provider ("openai" or "deepseek")

    Returns:
        LLMConfig with availability, api_key, provider, and model
    """
    if provider == "deepseek":
        api_key = get_deepseek_api_key()
        if api_key:
            return LLMConfig(
                available=True,
                api_key=api_key,
                provider="deepseek",
                model="deepseek-chat",
            )
        return LLMConfig(available=False)

    if provider == "openai" or provider is None:
        api_key = get_openai_api_key()
        if api_key:
            return LLMConfig(
                available=True,
                api_key=api_key,
                provider="openai",
                model="gpt-4o",
            )

    # Fallback to DeepSeek if OpenAI not available
    api_key = get_deepseek_api_key()
    if api_key:
        return LLMConfig(
            available=True,
            api_key=api_key,
            provider="deepseek",
            model="deepseek-chat",
        )

    return LLMConfig(available=False)
