# -*- coding: ascii -*-
"""Molecular weight filter plugin."""

import logging
from typing import Dict, Any, List

from src.halogenator.plugins.base import FilterPlugin, plugin, PluginType

LOG = logging.getLogger(__name__)


@plugin(PluginType.FILTER, "molecular_weight_filter", "1.0.0")
class MolecularWeightFilter(FilterPlugin):
    """Filter products based on molecular weight constraints."""
    
    @property
    def description(self) -> str:
        return "Filter products by molecular weight range"
    
    def can_handle(self, context: Dict[str, Any]) -> bool:
        # Can handle any context where molecular weight filtering might be useful
        return True
    
    def on_configure(self, config: Dict[str, Any]):
        """Configure molecular weight thresholds."""
        self.min_mw = config.get('min_molecular_weight', 0)
        self.max_mw = config.get('max_molecular_weight', 1000)
        LOG.info(f"Configured molecular weight filter: {self.min_mw} - {self.max_mw}")
    
    def filter_products(self, products: List[Dict[str, Any]], 
                       context: Dict[str, Any]) -> List[Dict[str, Any]]:
        """
        Filter products by molecular weight.
        
        Args:
            products: List of product dictionaries
            context: Context information
            
        Returns:
            Filtered list of products
        """
        if not products:
            return products
        
        # Get configured thresholds
        min_mw = self.config.get('min_molecular_weight', 0)
        max_mw = self.config.get('max_molecular_weight', 1000)
        
        filtered = []
        
        for product in products:
            # Try to get molecular weight from various possible fields
            mw = self._get_molecular_weight(product)
            
            if mw is not None and min_mw <= mw <= max_mw:
                filtered.append(product)
            elif mw is None:
                # If molecular weight is not available, include the product
                filtered.append(product)
        
        LOG.info(f"Molecular weight filter: {len(products)} -> {len(filtered)} products "
                   f"(range: {min_mw}-{max_mw})")
        
        return filtered
    
    def _get_molecular_weight(self, product: Dict[str, Any]) -> float:
        """
        Extract molecular weight from product data.
        
        Args:
            product: Product dictionary
            
        Returns:
            Molecular weight or None if not available
        """
        # Try various field names that might contain molecular weight
        mw_fields = [
            'molecular_weight', 'mw', 'mol_weight', 'molwt',
            'formula_weight', 'exact_mass', 'mass'
        ]
        
        for field in mw_fields:
            if field in product and product[field] is not None:
                try:
                    return float(product[field])
                except (ValueError, TypeError):
                    continue
        
        # Try to calculate from SMILES if RDKit is available
        smiles = product.get('smiles')
        if smiles:
            try:
                from rdkit import Chem
                from rdkit.Chem import Descriptors
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    return Descriptors.MolWt(mol)
            except ImportError:
                LOG.debug("RDKit not available for molecular weight calculation")
            except Exception as e:
                LOG.debug(f"Failed to calculate molecular weight for {smiles}: {e}")
        
        return None