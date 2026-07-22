import json
from pathlib import Path
def test_health_contract():
 root=Path(__file__).parents[1]
 for name in ("source-coverage.json","ontology-policy.json","build-status.json"):
  assert json.loads((root/"health"/name).read_text())
 assert (root/".bb/dependencies.txt").read_text().strip()=="zenodo-14409985"
