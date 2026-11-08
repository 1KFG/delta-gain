#!/usr/bin/env python3
"""ESMFold (HF transformers) on the pilot representatives; writes per-protein mean pLDDT and PDBs."""
import os
import sys
import time

import torch
from transformers import AutoTokenizer, EsmForProteinFolding
from transformers.models.esm.openfold_utils.protein import to_pdb, Protein as OFProtein
from transformers.models.esm.openfold_utils.feats import atom14_to_atom37

fa, out = sys.argv[1], sys.argv[2]
os.makedirs(out + "/pdb", exist_ok=True)
seqs = []
name = None
for l in open(fa):
    if l.startswith(">"):
        name = l[1:].split()[0]
        seqs.append([name, ""])
    else:
        seqs[-1][1] += l.strip()
tok = AutoTokenizer.from_pretrained("facebook/esmfold_v1")
model = EsmForProteinFolding.from_pretrained("facebook/esmfold_v1", low_cpu_mem_usage=True).cuda().eval()
model.esm = model.esm.half()
model.trunk.set_chunk_size(64)
res = open(out + "/plddt.tsv", "w")
res.write("id\tlength\tmean_plddt\tseconds\n")
t0 = time.time()
for i, (n, s) in enumerate(seqs):
    t = time.time()
    try:
        ids = tok([s], return_tensors="pt", add_special_tokens=False)["input_ids"].cuda()
        with torch.no_grad():
            o = model(ids)
        pl = o["plddt"][0, :, 1].float().mean().item() if o["plddt"].dim() == 3 else float(o["mean_plddt"])
        if i < 60:
            fin = {k: v.cpu().numpy() for k, v in o.items() if hasattr(v, "cpu")}
            from transformers.models.esm.openfold_utils.protein import Protein
            prot = Protein(aatype=fin["aatype"][0], atom_positions=fin["positions"][-1][0], atom_mask=fin["atom37_atom_exists"][0],
                           residue_index=fin["residue_index"][0] + 1, b_factors=fin["plddt"][0], chain_index=None)
            open(f"{out}/pdb/{n}.pdb", "w").write(to_pdb(prot))
        res.write(f"{n}\t{len(s)}\t{pl:.2f}\t{time.time() - t:.1f}\n")
        res.flush()
    except Exception as e:
        res.write(f"{n}\t{len(s)}\tERR\t{repr(e)[:80]}\n")
        torch.cuda.empty_cache()
print("done", len(seqs), "in", time.time() - t0)
