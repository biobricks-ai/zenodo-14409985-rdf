# zenodo-14409985-rdf

RDF conversion of the CC BY 4.0 deepFPlearn data release. The graph keeps three
scientifically different things separate:

- curated training labels are `BenchmarkLabel` records and are not promoted to
  primary experimental observations;
- Tox21 test labels are `Tox21AssayLabel` records tied to their named assay;
- CompTox predictions are `ModelPrediction` records with decimal probability and
  the deepFPlearn model family as provenance.

The identity-only autoencoder/CompTox tables and AR/ER holdout lists are audited
but do not create evidence. The duplicate `train_ae` and `comptox` tables are not
loaded twice. Compounds use PubChem CID, DSSTox substance identifiers, and
standard InChIKey IRIs. Tox21 SMILES are standardized with RDKit to InChIKey.
Targets use UniProt where the source label identifies a single protein; broad TR,
ED and assay labels remain explicit source endpoints rather than guessed proteins.

Run `uv run pytest` and then `dvc repro`. The source is Zenodo record 14409985,
licensed CC BY 4.0 and described by DOI 10.1093/bib/bbac257.
