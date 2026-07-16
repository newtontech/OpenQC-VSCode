"""
OpenQC AI Module

AI-powered features for quantum chemistry input file optimization,
generation, explanation, and debugging.
"""

from .parser import ResponseParser
from .prompts import PromptTemplates

__all__ = ["AIClient", "AIRequest", "AIResponse", "PromptTemplates", "ResponseParser"]


def __getattr__(name):
    """Load client types lazily so ``python -m openqc.ai.client`` stays warning-free."""
    if name in {"AIClient", "AIRequest", "AIResponse"}:
        from . import client

        return getattr(client, name)
    raise AttributeError(name)
