import importlib.util
from pathlib import Path
import pyarrow as pa, pyarrow.parquet as pq
from rdflib import Graph, URIRef
from rdflib.namespace import SKOS

ROOT=Path(__file__).parents[1]; spec=importlib.util.spec_from_file_location("convert",ROOT/"stages/convert.py"); C=importlib.util.module_from_spec(spec); spec.loader.exec_module(C)

def fixtures(root):
 pq.write_table(pa.Table.from_pylist([{"cid":8058.0,"smiles":"CCCCCC","AR":0.0,"ER":0.0,"GR":0.0,"Aromatase":0.0,"TR":None,"PPARg":0.0,"ED":0}]),root/"train_fnn.parquet")
 pq.write_table(pa.Table.from_pylist([{**{x:None for x in C.TOX21},"NR-AR":1.0,"smiles":"CCCCCC"}]),root/"test_tox21.parquet")
 pq.write_table(pa.Table.from_pylist([{"toxid":"DTXSID7020682","inchi":None,"key":"VLKZOEOYAKHREP-UHFFFAOYSA-N",**{x:0.25 for x in C.TARGETS}}]),root/"predictions_comptox.parquet")

def test_hexane_and_evidence_separation(tmp_path):
 fixtures(tmp_path); out=tmp_path/"x.nt"; report=C.convert(tmp_path,out,tmp_path/"r.json")
 g=Graph().parse(out,format="nt")
 assert report["training_labels"]==6 and report["tox21_test_labels"]==1 and report["model_predictions"]==7
 assert report["n_hexane_records"]==2
 assert (None,URIRef(C.BASE+"compound"),URIRef(C.PUBCHEM+"8058")) in g
 assert (None,URIRef(C.BASE+"compound"),URIRef(C.INCHIKEY+"VLKZOEOYAKHREP-UHFFFAOYSA-N")) in g
 assert (URIRef(C.PUBCHEM+"8058"),SKOS.exactMatch,URIRef(C.INCHIKEY+"VLKZOEOYAKHREP-UHFFFAOYSA-N")) in g
 assert (URIRef(C.DTXSID+"DTXSID7020682"),SKOS.exactMatch,URIRef(C.INCHIKEY+"VLKZOEOYAKHREP-UHFFFAOYSA-N")) in g
 assert len(list(g.subjects(URIRef(C.RDF+"type"),URIRef(C.BASE+"ModelPrediction"))))==7

def test_target_identifiers_and_invalid_smiles():
 assert C.TARGETS["AR"]=="https://identifiers.org/uniprot/P10275"
 assert C.TARGETS["Aromatase"]=="https://identifiers.org/uniprot/P11511"
 assert C.smiles_iri("CCCCCC")==C.INCHIKEY+"VLKZOEOYAKHREP-UHFFFAOYSA-N"
 assert C.smiles_iri("not smiles") is None

def test_bounded_real_source(tmp_path):
 source=Path("/mnt/raid2/biobricks/zenodo-14409985/brick")
 if not source.exists(): return
 r=C.convert(source,tmp_path/"sample.nt",tmp_path/"r.json",max_rows=100)
 assert r["scanned_rows"]==100 and r["record_count_coverage"]==1
