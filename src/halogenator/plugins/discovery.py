# -*- coding: ascii -*-
"""Plugin discovery system for automatic plugin loading."""

import os
import sys
import importlib
import importlib.util
import pkgutil
import logging
from pathlib import Path
from typing import List, Dict, Any, Optional, Type

from .base import PluginBase
from .registry import plugin_registry

LOG = logging.getLogger(__name__)


def discover_plugins(plugin_dirs: List[str] = None, 
                    package_prefixes: List[str] = None) -> Dict[str, Any]:
    """
    Discover and load plugins from various sources.
    
    Args:
        plugin_dirs: List of directories to search for plugin files
        package_prefixes: List of package prefixes to search for plugins
        
    Returns:
        Dictionary with discovery results and statistics
    """
    results = {
        'discovered': 0,
        'loaded': 0,
        'failed': 0,
        'errors': [],
        'plugins': []
    }
    
    # Default plugin directories
    if plugin_dirs is None:
        plugin_dirs = []
        
        # Add default plugin directories
        current_dir = Path(__file__).parent.parent
        default_dirs = [
            current_dir / "plugins" / "builtin",
            current_dir / "plugins" / "external",
            Path.cwd() / "plugins",
            Path.home() / ".halogenator" / "plugins"
        ]
        
        for plugin_dir in default_dirs:
            if plugin_dir.exists():
                plugin_dirs.append(str(plugin_dir))
    
    # Default package prefixes
    if package_prefixes is None:
        package_prefixes = [
            'halogenator_plugin_',
            'halogenator.plugins.',
            'halogenator_',
        ]
    
    # Discover from directories
    for plugin_dir in plugin_dirs:
        dir_results = discover_plugins_from_directory(plugin_dir)
        results['discovered'] += dir_results['discovered']
        results['loaded'] += dir_results['loaded']
        results['failed'] += dir_results['failed']
        results['errors'].extend(dir_results['errors'])
        results['plugins'].extend(dir_results['plugins'])
    
    # Discover from packages
    for prefix in package_prefixes:
        pkg_results = discover_plugins_from_packages(prefix)
        results['discovered'] += pkg_results['discovered']
        results['loaded'] += pkg_results['loaded']
        results['failed'] += pkg_results['failed']
        results['errors'].extend(pkg_results['errors'])
        results['plugins'].extend(pkg_results['plugins'])
    
    LOG.info(f"Plugin discovery complete: {results['loaded']} loaded, "
               f"{results['failed']} failed from {results['discovered']} discovered")
    
    return results


def discover_plugins_from_directory(plugin_dir: str) -> Dict[str, Any]:
    """
    Discover plugins from a directory.
    
    Args:
        plugin_dir: Directory to search for plugin files
        
    Returns:
        Dictionary with discovery results
    """
    results = {
        'discovered': 0,
        'loaded': 0,
        'failed': 0,
        'errors': [],
        'plugins': []
    }
    
    plugin_path = Path(plugin_dir)
    if not plugin_path.exists():
        LOG.debug(f"Plugin directory does not exist: {plugin_dir}")
        return results
    
    LOG.info(f"Discovering plugins in directory: {plugin_dir}")
    
    # Find all Python files
    for py_file in plugin_path.glob("**/*.py"):
        if py_file.name.startswith("__"):
            continue  # Skip __init__.py and __pycache__ etc.
        
        results['discovered'] += 1
        
        try:
            plugin_info = load_plugin_from_file(str(py_file))
            if plugin_info:
                results['loaded'] += 1
                results['plugins'].append(plugin_info)
            else:
                results['failed'] += 1
        except Exception as e:
            results['failed'] += 1
            error_msg = f"Failed to load plugin from {py_file}: {e}"
            results['errors'].append(error_msg)
            LOG.error(error_msg)
    
    return results


def discover_plugins_from_packages(package_prefix: str) -> Dict[str, Any]:
    """
    Discover plugins from installed packages.
    
    Args:
        package_prefix: Package prefix to search for
        
    Returns:
        Dictionary with discovery results
    """
    results = {
        'discovered': 0,
        'loaded': 0,
        'failed': 0,
        'errors': [],
        'plugins': []
    }
    
    LOG.info(f"Discovering plugins with package prefix: {package_prefix}")
    
    # Search for packages with the given prefix
    for finder, name, ispkg in pkgutil.iter_modules():
        if name.startswith(package_prefix):
            results['discovered'] += 1
            
            try:
                module = importlib.import_module(name)
                plugin_info = extract_plugins_from_module(module, name)
                if plugin_info:
                    results['loaded'] += 1
                    results['plugins'].extend(plugin_info)
                else:
                    results['failed'] += 1
            except Exception as e:
                results['failed'] += 1
                error_msg = f"Failed to load plugin package {name}: {e}"
                results['errors'].append(error_msg)
                LOG.error(error_msg)
    
    return results


def load_plugin_from_file(file_path: str) -> Optional[Dict[str, Any]]:
    """
    Load a plugin from a Python file.
    
    Args:
        file_path: Path to Python file containing plugin
        
    Returns:
        Plugin information dictionary or None if no plugins found
    """
    file_path = Path(file_path)
    module_name = f"plugin_{file_path.stem}_{id(file_path)}"
    
    try:
        # Load module from file
        spec = importlib.util.spec_from_file_location(module_name, file_path)
        if spec is None or spec.loader is None:
            LOG.warning(f"Could not create module spec for {file_path}")
            return None
        
        module = importlib.util.module_from_spec(spec)
        sys.modules[module_name] = module
        spec.loader.exec_module(module)
        
        # Extract plugins from module
        plugins = extract_plugins_from_module(module, str(file_path))
        
        if plugins:
            LOG.info(f"Loaded {len(plugins)} plugin(s) from {file_path}")
            return {
                'source': str(file_path),
                'type': 'file',
                'plugins': plugins
            }
        else:
            LOG.debug(f"No plugins found in {file_path}")
            return None
            
    except Exception as e:
        LOG.error(f"Failed to load plugin from {file_path}: {e}")
        raise


def extract_plugins_from_module(module, source_name: str) -> List[Dict[str, Any]]:
    """
    Extract plugin classes from a loaded module.
    
    Args:
        module: Loaded Python module
        source_name: Name/path of the source
        
    Returns:
        List of plugin information dictionaries
    """
    plugins = []
    
    # Look for classes that inherit from PluginBase
    for attr_name in dir(module):
        attr = getattr(module, attr_name)
        
        if (isinstance(attr, type) and 
            issubclass(attr, PluginBase) and 
            attr is not PluginBase):
            
            try:
                # Register the plugin class
                plugin_registry.register_class(attr)
                
                plugins.append({
                    'name': getattr(attr, '_plugin_name', attr.__name__),
                    'class': attr,
                    'source': source_name
                })
                
                LOG.debug(f"Found plugin class: {attr.__name__} in {source_name}")
                
            except Exception as e:
                LOG.error(f"Failed to register plugin class {attr.__name__}: {e}")
    
    # Look for plugin instances (for pre-instantiated plugins)
    for attr_name in dir(module):
        attr = getattr(module, attr_name)
        
        if isinstance(attr, PluginBase):
            try:
                # Register the plugin instance
                plugin_registry.register_instance(attr)
                
                plugins.append({
                    'name': attr.name,
                    'instance': attr,
                    'source': source_name
                })
                
                LOG.debug(f"Found plugin instance: {attr.name} in {source_name}")
                
            except Exception as e:
                LOG.error(f"Failed to register plugin instance {attr.name}: {e}")
    
    return plugins


def create_plugin_directory(plugin_dir: str):
    """
    Create a plugin directory with example files.
    
    Args:
        plugin_dir: Directory to create
    """
    plugin_path = Path(plugin_dir)
    plugin_path.mkdir(parents=True, exist_ok=True)
    
    # Create __init__.py
    init_file = plugin_path / "__init__.py"
    if not init_file.exists():
        init_file.write_text("""# -*- coding: ascii -*-
\"\"\"Custom halogenator plugins.\"\"\"
""")
    
    # Create example plugin
    example_file = plugin_path / "example_filter.py"
    if not example_file.exists():
        example_content = '''# -*- coding: ascii -*-
"""Example filter plugin for halogenator."""

from halogenator.plugins.base import FilterPlugin, plugin, PluginType


@plugin(PluginType.FILTER, "example_filter", "1.0.0")
class ExampleFilterPlugin(FilterPlugin):
    """Example plugin that filters products by molecular weight."""
    
    @property
    def description(self) -> str:
        return "Filter products by molecular weight threshold"
    
    def can_handle(self, context):
        return True  # Can handle any context
    
    def filter_products(self, products, context):
        # Example: filter by molecular weight if available
        threshold = self.config.get('max_mw', 500)
        
        filtered = []
        for product in products:
            mw = product.get('molecular_weight', 0)
            if mw <= threshold:
                filtered.append(product)
        
        return filtered
'''
        example_file.write_text(example_content)
    
    LOG.info(f"Created plugin directory: {plugin_dir}")


def reload_plugins():
    """Reload all plugins by clearing registry and rediscovering."""
    LOG.info("Reloading all plugins...")
    plugin_registry.clear()
    return discover_plugins()