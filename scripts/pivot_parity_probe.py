# scripts/pivot_parity_probe.py
import json
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))
from halogenator.enumerate_k import EnumConfig, enumerate_products, enumerate_with_stats

def run(smiles="c1ccccc1O"):
    cfg = EnumConfig(k_max=1, halogens=('F','Cl'))
    last = None
    prods = 0
    for prod, qa in enumerate_products(smiles, cfg, return_qa_stats=True):
        prods += 1
        last = qa
    _, snap = enumerate_with_stats(smiles, cfg)
    res = {
        "smiles": smiles,
        "n_streaming_products": prods,
        "streaming": {
            "qa_paths": last.get("qa_paths", {}),
            "dedup_top": {
                "dedup_hits_inchi": last.get("dedup_hits_inchi",0),
                "dedup_hits_statesig": last.get("dedup_hits_statesig",0),
            }
        },
        "snapshot": {
            "qa_paths": snap.get("qa_paths", {}),
            "dedup_top": {
                "dedup_hits_inchi": snap.get("dedup_hits_inchi",0),
                "dedup_hits_statesig": snap.get("dedup_hits_statesig",0),
            }
        }
    }
    print(json.dumps(res, indent=2, sort_keys=True))

if __name__ == "__main__":
    run()