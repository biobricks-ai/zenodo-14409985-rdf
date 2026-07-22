#!/usr/bin/env python3
"""Stream deepFPlearn labels and predictions into semantically separated RDF."""
from __future__ import annotations
import argparse, hashlib, json, os
from pathlib import Path
import pyarrow.parquet as pq
from rdkit import Chem

BASE = "https://biobricks.ai/zenodo-14409985-rdf/"
RDF = "http://www.w3.org/1999/02/22-rdf-syntax-ns#"
PROV = "http://www.w3.org/ns/prov#"
XSD = "http://www.w3.org/2001/XMLSchema#"
SKOS = "http://www.w3.org/2004/02/skos/core#"
PUBCHEM = "https://pubchem.ncbi.nlm.nih.gov/compound/"
DTXSID = "https://comptox.epa.gov/dashboard/chemical/details/"
INCHIKEY = "https://identifiers.org/inchikey/"
UNIPROT = "https://identifiers.org/uniprot/"
SOURCE = "https://doi.org/10.5281/zenodo.14409985"
MODEL = BASE + "model/deepFPlearn"
TARGETS = {"AR": UNIPROT+"P10275", "ER": UNIPROT+"P03372", "GR": UNIPROT+"P04150",
           "Aromatase": UNIPROT+"P11511", "PPARg": UNIPROT+"P37231",
           "TR": BASE+"target/thyroid-receptor-family", "ED": BASE+"endpoint/endocrine-disruption"}
TOX21 = ["NR-AR","NR-AR-LBD","NR-AhR","NR-Aromatase","NR-ER","NR-ER-LBD",
         "NR-PPAR-gamma","SR-ARE","SR-ATAD5","SR-HSE","SR-MMP","SR-p53"]

def iri(x): return f"<{x}>"
def lit(x, dtype=None):
    s=str(x).replace("\\","\\\\").replace('"','\\"').replace("\n","\\n")
    return f'"{s}"' + (f"^^<{dtype}>" if dtype else "")
def emit(f,s,p,o): f.write(f"{iri(s)} {iri(p)} {o if o.startswith(('<','\"')) else iri(o)} .\n")
def obs_id(kind, compound, endpoint): return BASE+"observation/"+hashlib.sha256(f"{kind}|{compound}|{endpoint}".encode()).hexdigest()[:32]
def smiles_iri(smiles):
    mol=Chem.MolFromSmiles(smiles or "")
    return INCHIKEY+Chem.MolToInchiKey(mol) if mol else None
def roots():
    root=Path(os.environ.get("ZENODO_14409985_BRICK", ".bb/dependencies/zenodo-14409985/brick"))
    if not root.exists(): root=Path("/mnt/raid2/biobricks/zenodo-14409985/brick")
    return root

def convert(root: Path, output: Path, report_path: Path, max_rows=None):
    c={"scanned_rows":0,"converted_records":0,"excluded_rows":0,"emitted_statements":0,
       "training_labels":0,"tox21_test_labels":0,"model_predictions":0,"identifier_rows":0,
       "source_nonempty_cells":0,"mapped_nonempty_cells":0,"n_hexane_records":0}
    output.parent.mkdir(parents=True,exist_ok=True)
    def stmt(f,s,p,o): emit(f,s,p,o); c["emitted_statements"]+=1
    with output.open("w") as f:
      stop=False
      for filename,kind,fields in [("train_fnn.parquet","BenchmarkLabel",list(TARGETS)),
                                   ("test_tox21.parquet","Tox21AssayLabel",TOX21),
                                   ("predictions_comptox.parquet","ModelPrediction",list(TARGETS))]:
        for batch in pq.ParquetFile(root/filename).iter_batches(batch_size=4096):
          for row in batch.to_pylist():
            if max_rows is not None and c["scanned_rows"]>=max_rows: stop=True; break
            c["scanned_rows"]+=1; c["source_nonempty_cells"]+=sum(v is not None for v in row.values())
            if kind=="BenchmarkLabel": compound=PUBCHEM+str(int(row["cid"])) if row.get("cid") is not None else None
            elif kind=="Tox21AssayLabel": compound=smiles_iri(row.get("smiles"))
            else:
              compound=DTXSID+str(row.get("toxid")) if row.get("toxid") else None
            if not compound: c["excluded_rows"]+=1; continue
            c["identifier_rows"]+=1
            if compound in (PUBCHEM+"8058", DTXSID+"DTXSID7020682"): c["n_hexane_records"]+=1
            # Preserve the source's standard structure identifier as an identity
            # bridge without treating raw InChI text as an RDF entity.
            key = row.get("key")
            if not key and kind=="BenchmarkLabel": key_iri=smiles_iri(row.get("smiles")); key=key_iri.removeprefix(INCHIKEY) if key_iri else None
            if key: stmt(f,compound,SKOS+"exactMatch",iri(INCHIKEY+str(key)))
            values=0
            for field in fields:
              value=row.get(field)
              if value is None: continue
              endpoint=TARGETS.get(field, BASE+"tox21-assay/"+field)
              obs=obs_id(kind,compound,endpoint); values+=1
              stmt(f,obs,RDF+"type",iri(BASE+kind)); stmt(f,obs,BASE+"compound",iri(compound))
              stmt(f,obs,BASE+"endpoint",iri(endpoint)); stmt(f,obs,PROV+"wasDerivedFrom",iri(SOURCE))
              pred=BASE+("probability" if kind=="ModelPrediction" else "activityLabel")
              dtype=XSD+("decimal" if kind=="ModelPrediction" else "integer")
              stmt(f,obs,pred,lit(value,dtype))
              if kind=="ModelPrediction": stmt(f,obs,PROV+"wasGeneratedBy",iri(MODEL)); c["model_predictions"]+=1
              elif kind=="BenchmarkLabel": c["training_labels"]+=1
              else: c["tox21_test_labels"]+=1
            if values: c["converted_records"]+=1; c["mapped_nonempty_cells"]+=values+1
          if stop: break
        if stop: break
    total=sum(pq.ParquetFile(root/f).metadata.num_rows for f in ("train_fnn.parquet","test_tox21.parquet","predictions_comptox.parquet"))
    c.update({"source_rows":total,"source_tables":7,"bounded_run":max_rows is not None,
      "scan_coverage":c["scanned_rows"]/total,"record_count_coverage":c["converted_records"]/(c["scanned_rows"]-c["excluded_rows"]) if c["scanned_rows"]>c["excluded_rows"] else 1,
      "identifier_row_coverage":c["identifier_rows"]/c["scanned_rows"] if c["scanned_rows"] else 1,
      "mapped_cell_coverage":c["mapped_nonempty_cells"]/c["source_nonempty_cells"] if c["source_nonempty_cells"] else 1,
      "explicit_exclusions":{"identity_only_train_ae":719997,"duplicate_comptox_identity_table":719997,"identity_only_ar_holdout":4082,"identity_only_er_holdout":4060},
      "license":"CC BY 4.0"})
    report_path.parent.mkdir(parents=True,exist_ok=True); report_path.write_text(json.dumps(c,indent=2,sort_keys=True)+"\n")
    return c

if __name__=="__main__":
 p=argparse.ArgumentParser(); p.add_argument("--source",type=Path); p.add_argument("--output",type=Path,default=Path("brick/zenodo-14409985-rdf.nt")); p.add_argument("--report",type=Path,default=Path("reports/source-coverage.json")); p.add_argument("--max-rows",type=int)
 a=p.parse_args(); print(json.dumps(convert(a.source or roots(),a.output,a.report,a.max_rows),indent=2))
