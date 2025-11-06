# -*- coding: ascii -*-
"""Base classes and interfaces for halogenator plugins."""

import abc
from typing import Dict, Any, List, Optional, Callable, Tuple
from enum import Enum


class PluginType(Enum):
    """Plugin types for different extension points."""
    FILTER = "filter"                    # Product filtering plugins
    VALIDATOR = "validator"              # Validation plugins  
    EXPORTER = "exporter"               # Output format plugins
    ANALYZER = "analyzer"               # Analysis plugins
    CONSTRAINT = "constraint"           # Constraint plugins
    RULE = "rule"                       # Reaction rule plugins
    CLI_COMMAND = "cli_command"         # CLI command plugins
    PREPROCESSOR = "preprocessor"       # Data preprocessing plugins
    POSTPROCESSOR = "postprocessor"     # Data postprocessing plugins


class PluginBase(abc.ABC):
    """
    Base class for all halogenator plugins.
    
    All plugins must inherit from this class and implement the required methods.
    """
    
    def __init__(self, name: str, version: str = "1.0.0"):
        """
        Initialize plugin.
        
        Args:
            name: Plugin name (should be unique)
            version: Plugin version string
        """
        self.name = name
        self.version = version
        self._enabled = True
        self._config = {}
    
    @property
    @abc.abstractmethod
    def plugin_type(self) -> PluginType:
        """Return the type of this plugin."""
        pass
    
    @property
    @abc.abstractmethod
    def description(self) -> str:
        """Return a description of what this plugin does."""
        pass
    
    @property
    def enabled(self) -> bool:
        """Return whether this plugin is enabled."""
        return self._enabled
    
    @enabled.setter
    def enabled(self, value: bool):
        """Enable or disable this plugin."""
        self._enabled = value
    
    @property
    def config(self) -> Dict[str, Any]:
        """Return plugin configuration."""
        return self._config.copy()
    
    def configure(self, config: Dict[str, Any]):
        """
        Configure the plugin with settings.
        
        Args:
            config: Configuration dictionary
        """
        self._config.update(config)
        self.on_configure(config)
    
    def on_configure(self, config: Dict[str, Any]):
        """
        Called when plugin is configured. Override in subclasses.
        
        Args:
            config: Configuration dictionary
        """
        pass
    
    @abc.abstractmethod
    def can_handle(self, context: Dict[str, Any]) -> bool:
        """
        Check if this plugin can handle the given context.
        
        Args:
            context: Context information for the operation
            
        Returns:
            True if plugin can handle this context
        """
        pass
    
    def get_metadata(self) -> Dict[str, Any]:
        """
        Return plugin metadata.
        
        Returns:
            Dictionary with plugin metadata
        """
        return {
            'name': self.name,
            'version': self.version,
            'type': self.plugin_type.value,
            'description': self.description,
            'enabled': self.enabled,
            'config': self.config
        }


class FilterPlugin(PluginBase):
    """Base class for product filtering plugins."""
    
    @property
    def plugin_type(self) -> PluginType:
        return PluginType.FILTER
    
    @abc.abstractmethod
    def filter_products(self, products: List[Dict[str, Any]], 
                       context: Dict[str, Any]) -> List[Dict[str, Any]]:
        """
        Filter products based on plugin criteria.
        
        Args:
            products: List of product dictionaries
            context: Context information
            
        Returns:
            Filtered list of products
        """
        pass


class ValidatorPlugin(PluginBase):
    """Base class for validation plugins."""
    
    @property
    def plugin_type(self) -> PluginType:
        return PluginType.VALIDATOR
    
    @abc.abstractmethod
    def validate(self, data: Any, context: Dict[str, Any]) -> List[str]:
        """
        Validate data and return any validation errors.
        
        Args:
            data: Data to validate
            context: Context information
            
        Returns:
            List of validation error messages (empty if valid)
        """
        pass


class ExporterPlugin(PluginBase):
    """Base class for output format plugins."""
    
    @property
    def plugin_type(self) -> PluginType:
        return PluginType.EXPORTER
    
    @property
    @abc.abstractmethod
    def supported_formats(self) -> List[str]:
        """Return list of file formats this exporter supports."""
        pass
    
    @abc.abstractmethod
    def export(self, data: Any, output_path: str, 
              context: Dict[str, Any]) -> bool:
        """
        Export data to specified output path.
        
        Args:
            data: Data to export
            output_path: Path to write output file
            context: Context information
            
        Returns:
            True if export was successful
        """
        pass


class AnalyzerPlugin(PluginBase):
    """Base class for analysis plugins."""
    
    @property
    def plugin_type(self) -> PluginType:
        return PluginType.ANALYZER
    
    @abc.abstractmethod
    def analyze(self, data: Any, context: Dict[str, Any]) -> Dict[str, Any]:
        """
        Analyze data and return analysis results.
        
        Args:
            data: Data to analyze
            context: Context information
            
        Returns:
            Dictionary with analysis results
        """
        pass


class ConstraintPlugin(PluginBase):
    """Base class for constraint plugins."""
    
    @property
    def plugin_type(self) -> PluginType:
        return PluginType.CONSTRAINT
    
    @abc.abstractmethod
    def check_constraint(self, mol: Any, history: List[Dict[str, Any]], 
                        context: Dict[str, Any]) -> Tuple[bool, Dict[str, Any]]:
        """
        Check if molecule satisfies the constraint.
        
        Args:
            mol: Molecule object
            history: Substitution history
            context: Context information
            
        Returns:
            Tuple of (constraint_satisfied, violation_details)
        """
        pass


class CLICommandPlugin(PluginBase):
    """Base class for CLI command plugins."""
    
    @property
    def plugin_type(self) -> PluginType:
        return PluginType.CLI_COMMAND
    
    @property
    @abc.abstractmethod
    def command_name(self) -> str:
        """Return the name of the CLI command."""
        pass
    
    @abc.abstractmethod
    def add_arguments(self, parser):
        """
        Add arguments to the CLI parser.
        
        Args:
            parser: argparse subparser for this command
        """
        pass
    
    @abc.abstractmethod
    def execute(self, args: Any) -> int:
        """
        Execute the CLI command.
        
        Args:
            args: Parsed command line arguments
            
        Returns:
            Exit code (0 for success)
        """
        pass


class PreprocessorPlugin(PluginBase):
    """Base class for data preprocessing plugins."""
    
    @property
    def plugin_type(self) -> PluginType:
        return PluginType.PREPROCESSOR
    
    @abc.abstractmethod
    def preprocess(self, data: Any, context: Dict[str, Any]) -> Any:
        """
        Preprocess data before main processing.
        
        Args:
            data: Input data
            context: Context information
            
        Returns:
            Preprocessed data
        """
        pass


class PostprocessorPlugin(PluginBase):
    """Base class for data postprocessing plugins."""
    
    @property
    def plugin_type(self) -> PluginType:
        return PluginType.POSTPROCESSOR
    
    @abc.abstractmethod
    def postprocess(self, data: Any, context: Dict[str, Any]) -> Any:
        """
        Postprocess data after main processing.
        
        Args:
            data: Output data from main processing
            context: Context information
            
        Returns:
            Postprocessed data
        """
        pass


def plugin(plugin_type: PluginType, name: str = None, version: str = "1.0.0"):
    """
    Decorator to register a class as a plugin.
    
    Args:
        plugin_type: Type of plugin
        name: Plugin name (defaults to class name)
        version: Plugin version
        
    Returns:
        Decorated class
    """
    def decorator(cls):
        cls._plugin_type = plugin_type
        cls._plugin_name = name or cls.__name__
        cls._plugin_version = version
        
        # Auto-register with registry when class is defined
        from .registry import plugin_registry
        plugin_registry.register_class(cls)
        
        return cls
    
    return decorator
