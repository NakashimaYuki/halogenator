# -*- coding: ascii -*-
"""Plugin framework for halogenator extensibility."""

from .base import PluginBase, PluginType
from .manager import PluginManager, plugin_manager
from .discovery import discover_plugins, load_plugin_from_file
from .registry import plugin_registry

__all__ = [
    'PluginBase',
    'PluginType', 
    'PluginManager',
    'plugin_manager',
    'discover_plugins',
    'load_plugin_from_file',
    'plugin_registry'
]