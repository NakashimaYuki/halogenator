# -*- coding: ascii -*-
from contextlib import contextmanager
import os

def rdkit_guard_enabled() -> bool:
    return os.environ.get('HALO_RDKIT_GUARD', '0') == '1'

@contextmanager
def rdkit_guard(qa_dict: dict):
    """Catch RDKit-related exceptions and convert to QA counters."""
    try:
        yield
    except Exception as e:
        if rdkit_guard_enabled():
            # record as RDKit error path; do not crash
            qa_paths = qa_dict.setdefault('qa_paths', {})
            qa_paths['rdkit_error'] = qa_paths.get('rdkit_error', 0) + 1
            # also count as template_unsupported if applicable
            qa_dict['template_unsupported'] = qa_dict.get('template_unsupported', 0) + 1
        else:
            raise