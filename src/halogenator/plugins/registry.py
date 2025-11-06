# -*- coding: ascii -*-
"""Plugin registry for managing plugin instances and classes."""

import logging
from typing import Dict, List, Type, Any, Optional
from collections import defaultdict

from .base import PluginBase, PluginType

LOG = logging.getLogger(__name__)


class PluginRegistry:
    """
    Central registry for managing plugins.
    
    Handles plugin registration, discovery, and lifecycle management.
    """
    
    def __init__(self):
        """Initialize empty registry."""
        self._plugin_classes: Dict[str, Type[PluginBase]] = {}
        self._plugin_instances: Dict[str, PluginBase] = {}
        self._plugins_by_type: Dict[PluginType, List[str]] = defaultdict(list)
        
    def register_class(self, plugin_class: Type[PluginBase], name: str = None):
        """
        Register a plugin class.
        
        Args:
            plugin_class: Plugin class to register
            name: Plugin name (defaults to class name)
        """
        if not issubclass(plugin_class, PluginBase):
            raise ValueError(f"Plugin class {plugin_class} must inherit from PluginBase")
        
        plugin_name = name or getattr(plugin_class, '_plugin_name', plugin_class.__name__)
        
        if plugin_name in self._plugin_classes:
            LOG.warning(f"Plugin '{plugin_name}' is already registered, overwriting")
        
        self._plugin_classes[plugin_name] = plugin_class
        LOG.info(f"Registered plugin class: {plugin_name}")
    
    def register_instance(self, plugin_instance: PluginBase):
        """
        Register a plugin instance.
        
        Args:
            plugin_instance: Plugin instance to register
        """
        if not isinstance(plugin_instance, PluginBase):
            raise ValueError(f"Plugin instance must inherit from PluginBase")
        
        plugin_name = plugin_instance.name
        
        if plugin_name in self._plugin_instances:
            LOG.warning(f"Plugin instance '{plugin_name}' is already registered, overwriting")
        
        self._plugin_instances[plugin_name] = plugin_instance
        
        # Add to type mapping
        plugin_type = plugin_instance.plugin_type
        if plugin_name not in self._plugins_by_type[plugin_type]:
            self._plugins_by_type[plugin_type].append(plugin_name)
        
        LOG.info(f"Registered plugin instance: {plugin_name} ({plugin_type.value})")
    
    def create_instance(self, plugin_name: str, *args, **kwargs) -> Optional[PluginBase]:
        """
        Create an instance of a registered plugin class.
        
        Args:
            plugin_name: Name of plugin to instantiate
            *args, **kwargs: Arguments to pass to plugin constructor
            
        Returns:
            Plugin instance or None if not found
        """
        if plugin_name not in self._plugin_classes:
            LOG.error(f"Plugin class '{plugin_name}' not found in registry")
            return None
        
        try:
            plugin_class = self._plugin_classes[plugin_name]
            instance = plugin_class(*args, **kwargs)
            self.register_instance(instance)
            return instance
        except Exception as e:
            LOG.error(f"Failed to create instance of plugin '{plugin_name}': {e}")
            return None
    
    def get_instance(self, plugin_name: str) -> Optional[PluginBase]:
        """
        Get a registered plugin instance.
        
        Args:
            plugin_name: Name of plugin
            
        Returns:
            Plugin instance or None if not found
        """
        return self._plugin_instances.get(plugin_name)
    
    def get_plugins_by_type(self, plugin_type: PluginType) -> List[PluginBase]:
        """
        Get all plugin instances of a specific type.
        
        Args:
            plugin_type: Type of plugins to retrieve
            
        Returns:
            List of plugin instances
        """
        plugin_names = self._plugins_by_type.get(plugin_type, [])
        return [self._plugin_instances[name] for name in plugin_names 
                if name in self._plugin_instances]
    
    def get_enabled_plugins_by_type(self, plugin_type: PluginType) -> List[PluginBase]:
        """
        Get all enabled plugin instances of a specific type.
        
        Args:
            plugin_type: Type of plugins to retrieve
            
        Returns:
            List of enabled plugin instances
        """
        return [plugin for plugin in self.get_plugins_by_type(plugin_type) 
                if plugin.enabled]
    
    def list_all_plugins(self) -> Dict[str, Dict[str, Any]]:
        """
        List all registered plugins with their metadata.
        
        Returns:
            Dictionary mapping plugin names to metadata
        """
        all_plugins = {}
        
        # Add instances
        for name, instance in self._plugin_instances.items():
            all_plugins[name] = instance.get_metadata()
            all_plugins[name]['status'] = 'instantiated'
        
        # Add classes that haven't been instantiated
        for name, plugin_class in self._plugin_classes.items():
            if name not in all_plugins:
                # Create temporary instance to get metadata
                try:
                    temp_instance = plugin_class(name)
                    metadata = temp_instance.get_metadata()
                    metadata['status'] = 'registered'
                    all_plugins[name] = metadata
                except Exception as e:
                    all_plugins[name] = {
                        'name': name,
                        'type': 'unknown',
                        'status': 'error',
                        'error': str(e)
                    }
        
        return all_plugins
    
    def unregister(self, plugin_name: str):
        """
        Unregister a plugin.
        
        Args:
            plugin_name: Name of plugin to unregister
        """
        # Remove instance
        if plugin_name in self._plugin_instances:
            instance = self._plugin_instances.pop(plugin_name)
            plugin_type = instance.plugin_type
            if plugin_name in self._plugins_by_type[plugin_type]:
                self._plugins_by_type[plugin_type].remove(plugin_name)
            LOG.info(f"Unregistered plugin instance: {plugin_name}")
        
        # Remove class
        if plugin_name in self._plugin_classes:
            self._plugin_classes.pop(plugin_name)
            LOG.info(f"Unregistered plugin class: {plugin_name}")
    
    def clear(self):
        """Clear all registered plugins."""
        self._plugin_classes.clear()
        self._plugin_instances.clear()
        self._plugins_by_type.clear()
        LOG.info("Cleared all plugins from registry")
    
    def enable_plugin(self, plugin_name: str):
        """
        Enable a plugin.
        
        Args:
            plugin_name: Name of plugin to enable
        """
        instance = self.get_instance(plugin_name)
        if instance:
            instance.enabled = True
            LOG.info(f"Enabled plugin: {plugin_name}")
        else:
            LOG.warning(f"Plugin not found: {plugin_name}")
    
    def disable_plugin(self, plugin_name: str):
        """
        Disable a plugin.
        
        Args:
            plugin_name: Name of plugin to disable
        """
        instance = self.get_instance(plugin_name)
        if instance:
            instance.enabled = False
            LOG.info(f"Disabled plugin: {plugin_name}")
        else:
            LOG.warning(f"Plugin not found: {plugin_name}")
    
    def configure_plugin(self, plugin_name: str, config: Dict[str, Any]):
        """
        Configure a plugin.
        
        Args:
            plugin_name: Name of plugin to configure
            config: Configuration dictionary
        """
        instance = self.get_instance(plugin_name)
        if instance:
            instance.configure(config)
            LOG.info(f"Configured plugin: {plugin_name}")
        else:
            LOG.warning(f"Plugin not found: {plugin_name}")


# Global plugin registry instance
plugin_registry = PluginRegistry()