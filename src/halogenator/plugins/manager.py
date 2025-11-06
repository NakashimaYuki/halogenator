# -*- coding: ascii -*-
"""Plugin manager for coordinating plugin lifecycle and execution."""

import logging
from typing import Dict, List, Any, Optional, Type, Callable, Tuple
from pathlib import Path

from .base import PluginBase, PluginType
from .registry import plugin_registry
from .discovery import discover_plugins

LOG = logging.getLogger(__name__)


class PluginManager:
    """
    Central manager for plugin system.
    
    Coordinates plugin discovery, loading, configuration, and execution.
    """
    
    def __init__(self, auto_discover: bool = True):
        """
        Initialize plugin manager.
        
        Args:
            auto_discover: Whether to automatically discover plugins on initialization
        """
        self._initialized = False
        self._config = {}
        
        if auto_discover:
            self.initialize()
    
    def initialize(self, plugin_dirs: List[str] = None, 
                  config: Dict[str, Any] = None):
        """
        Initialize the plugin system.
        
        Args:
            plugin_dirs: Additional plugin directories to search
            config: Plugin configuration dictionary
        """
        if self._initialized:
            LOG.warning("Plugin manager already initialized")
            return
        
        LOG.info("Initializing plugin manager...")
        
        # Store configuration
        self._config = config or {}
        
        # Discover and load plugins
        discovery_results = discover_plugins(plugin_dirs)
        
        # Configure plugins from config
        if config:
            self._apply_configuration(config)
        
        self._initialized = True
        
        LOG.info(f"Plugin manager initialized. "
                   f"Loaded {discovery_results['loaded']} plugins, "
                   f"failed {discovery_results['failed']}")
    
    def _apply_configuration(self, config: Dict[str, Any]):
        """Apply configuration to plugins."""
        plugins_config = config.get('plugins', {})
        
        for plugin_name, plugin_config in plugins_config.items():
            if isinstance(plugin_config, dict):
                # Configure plugin
                if 'enabled' in plugin_config:
                    if plugin_config['enabled']:
                        plugin_registry.enable_plugin(plugin_name)
                    else:
                        plugin_registry.disable_plugin(plugin_name)
                
                if 'config' in plugin_config:
                    plugin_registry.configure_plugin(plugin_name, plugin_config['config'])
    
    def get_plugins(self, plugin_type: PluginType, 
                   enabled_only: bool = True) -> List[PluginBase]:
        """
        Get plugins of a specific type.
        
        Args:
            plugin_type: Type of plugins to retrieve
            enabled_only: Whether to return only enabled plugins
            
        Returns:
            List of plugin instances
        """
        if enabled_only:
            return plugin_registry.get_enabled_plugins_by_type(plugin_type)
        else:
            return plugin_registry.get_plugins_by_type(plugin_type)
    
    def execute_filters(self, products: List[Dict[str, Any]], 
                       context: Dict[str, Any] = None) -> List[Dict[str, Any]]:
        """
        Execute filter plugins on products.
        
        Args:
            products: List of product dictionaries
            context: Context information
            
        Returns:
            Filtered products
        """
        context = context or {}
        filters = self.get_plugins(PluginType.FILTER)
        
        filtered_products = products
        
        for filter_plugin in filters:
            if filter_plugin.can_handle(context):
                try:
                    filtered_products = filter_plugin.filter_products(filtered_products, context)
                    LOG.debug(f"Filter '{filter_plugin.name}' processed {len(products)} -> {len(filtered_products)} products")
                except Exception as e:
                    LOG.error(f"Filter plugin '{filter_plugin.name}' failed: {e}")
        
        return filtered_products
    
    def execute_validators(self, data: Any, 
                          context: Dict[str, Any] = None) -> List[str]:
        """
        Execute validator plugins on data.
        
        Args:
            data: Data to validate
            context: Context information
            
        Returns:
            List of validation error messages
        """
        context = context or {}
        validators = self.get_plugins(PluginType.VALIDATOR)
        
        all_errors = []
        
        for validator in validators:
            if validator.can_handle(context):
                try:
                    errors = validator.validate(data, context)
                    all_errors.extend(errors)
                    LOG.debug(f"Validator '{validator.name}' found {len(errors)} errors")
                except Exception as e:
                    LOG.error(f"Validator plugin '{validator.name}' failed: {e}")
                    all_errors.append(f"Validator '{validator.name}' failed: {e}")
        
        return all_errors
    
    def execute_analyzers(self, data: Any, 
                         context: Dict[str, Any] = None) -> Dict[str, Any]:
        """
        Execute analyzer plugins on data.
        
        Args:
            data: Data to analyze
            context: Context information
            
        Returns:
            Combined analysis results
        """
        context = context or {}
        analyzers = self.get_plugins(PluginType.ANALYZER)
        
        combined_results = {}
        
        for analyzer in analyzers:
            if analyzer.can_handle(context):
                try:
                    results = analyzer.analyze(data, context)
                    combined_results[analyzer.name] = results
                    LOG.debug(f"Analyzer '{analyzer.name}' completed analysis")
                except Exception as e:
                    LOG.error(f"Analyzer plugin '{analyzer.name}' failed: {e}")
                    combined_results[analyzer.name] = {'error': str(e)}
        
        return combined_results
    
    def execute_constraints(self, mol: Any, history: List[Dict[str, Any]], 
                           context: Dict[str, Any] = None) -> Tuple[bool, Dict[str, Any]]:
        """
        Execute constraint plugins.
        
        Args:
            mol: Molecule object
            history: Substitution history
            context: Context information
            
        Returns:
            Tuple of (all_constraints_satisfied, violation_details)
        """
        context = context or {}
        constraints = self.get_plugins(PluginType.CONSTRAINT)
        
        all_satisfied = True
        all_violations = {}
        
        for constraint in constraints:
            if constraint.can_handle(context):
                try:
                    satisfied, violations = constraint.check_constraint(mol, history, context)
                    if not satisfied:
                        all_satisfied = False
                        all_violations[constraint.name] = violations
                    LOG.debug(f"Constraint '{constraint.name}' satisfied: {satisfied}")
                except Exception as e:
                    LOG.error(f"Constraint plugin '{constraint.name}' failed: {e}")
                    all_satisfied = False
                    all_violations[constraint.name] = {'error': str(e)}
        
        return all_satisfied, all_violations
    
    def export_data(self, data: Any, output_path: str, 
                   context: Dict[str, Any] = None) -> bool:
        """
        Export data using exporter plugins.
        
        Args:
            data: Data to export
            output_path: Output file path
            context: Context information
            
        Returns:
            True if export was successful
        """
        context = context or {}
        exporters = self.get_plugins(PluginType.EXPORTER)
        
        # Find exporter that supports the file format
        file_ext = Path(output_path).suffix.lower()
        
        for exporter in exporters:
            if (exporter.can_handle(context) and 
                any(fmt.lower() == file_ext for fmt in exporter.supported_formats)):
                try:
                    success = exporter.export(data, output_path, context)
                    if success:
                        LOG.info(f"Data exported using '{exporter.name}' to {output_path}")
                        return True
                    else:
                        LOG.warning(f"Exporter '{exporter.name}' failed to export data")
                except Exception as e:
                    LOG.error(f"Exporter plugin '{exporter.name}' failed: {e}")
        
        LOG.warning(f"No suitable exporter found for format: {file_ext}")
        return False
    
    def preprocess_data(self, data: Any, 
                       context: Dict[str, Any] = None) -> Any:
        """
        Preprocess data using preprocessor plugins.
        
        Args:
            data: Input data
            context: Context information
            
        Returns:
            Preprocessed data
        """
        context = context or {}
        preprocessors = self.get_plugins(PluginType.PREPROCESSOR)
        
        processed_data = data
        
        for preprocessor in preprocessors:
            if preprocessor.can_handle(context):
                try:
                    processed_data = preprocessor.preprocess(processed_data, context)
                    LOG.debug(f"Preprocessor '{preprocessor.name}' completed")
                except Exception as e:
                    LOG.error(f"Preprocessor plugin '{preprocessor.name}' failed: {e}")
        
        return processed_data
    
    def postprocess_data(self, data: Any, 
                        context: Dict[str, Any] = None) -> Any:
        """
        Postprocess data using postprocessor plugins.
        
        Args:
            data: Output data
            context: Context information
            
        Returns:
            Postprocessed data
        """
        context = context or {}
        postprocessors = self.get_plugins(PluginType.POSTPROCESSOR)
        
        processed_data = data
        
        for postprocessor in postprocessors:
            if postprocessor.can_handle(context):
                try:
                    processed_data = postprocessor.postprocess(processed_data, context)
                    LOG.debug(f"Postprocessor '{postprocessor.name}' completed")
                except Exception as e:
                    LOG.error(f"Postprocessor plugin '{postprocessor.name}' failed: {e}")
        
        return processed_data
    
    def list_plugins(self) -> Dict[str, Dict[str, Any]]:
        """
        List all available plugins.
        
        Returns:
            Dictionary of plugin information
        """
        return plugin_registry.list_all_plugins()
    
    def enable_plugin(self, plugin_name: str):
        """Enable a plugin."""
        plugin_registry.enable_plugin(plugin_name)
    
    def disable_plugin(self, plugin_name: str):
        """Disable a plugin."""
        plugin_registry.disable_plugin(plugin_name)
    
    def configure_plugin(self, plugin_name: str, config: Dict[str, Any]):
        """Configure a plugin."""
        plugin_registry.configure_plugin(plugin_name, config)
    
    def get_plugin_info(self, plugin_name: str) -> Optional[Dict[str, Any]]:
        """
        Get information about a specific plugin.
        
        Args:
            plugin_name: Name of plugin
            
        Returns:
            Plugin information dictionary or None if not found
        """
        all_plugins = self.list_plugins()
        return all_plugins.get(plugin_name)
    
    def reload_plugins(self):
        """Reload all plugins."""
        LOG.info("Reloading plugins...")
        plugin_registry.clear()
        self._initialized = False
        self.initialize()


# Global plugin manager instance
plugin_manager = PluginManager(auto_discover=False)  # Don't auto-discover to avoid import issues
