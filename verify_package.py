#!/usr/bin/env python3
"""Verify every file in the release manifest using SHA-256 (standard library only)."""
from pathlib import Path
import csv,hashlib,json,sys
root=Path(__file__).resolve().parent
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()
with (root/'manifest.csv').open(encoding='utf-8-sig',newline='') as f:rows=list(csv.DictReader(f))
errors=[]
for row in rows:
    p=(root/row['path']).resolve()
    if root not in p.parents:errors.append({'path':row['path'],'problem':'outside package'});continue
    if not p.is_file():errors.append({'path':row['path'],'problem':'missing'});continue
    if p.stat().st_size!=int(row['bytes']) or sha(p)!=row['sha256']:errors.append({'path':row['path'],'problem':'size/hash mismatch'})
print(json.dumps({'passed':not errors,'files_checked':len(rows),'errors':errors},indent=2))
sys.exit(1 if errors else 0)
