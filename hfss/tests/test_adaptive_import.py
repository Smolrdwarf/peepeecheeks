"""Checks that HFSS adaptive-mesh settings in an .aedt driven setup map onto FieldForge's adaptive passes.

Run:  python hfss/tests/test_adaptive_import.py   (no solver, no HFSS needed)
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from fieldforge.aedt import importer as I          # noqa: E402
from fieldforge.aedt.parser import parse            # noqa: E402
from fieldforge.model import Project, Setup         # noqa: E402

AEDT_SETUP = """$begin 'AnsoftProject'
$begin 'HFSSModel'
$begin 'AnalysisSetup'
$begin 'SolveSetups'
$begin 'Setup1'
SetupType='HfssDriven'
Frequency='10GHz'
MaxDeltaS=0.01
MaximumPasses=12
MinimumPasses=1
MinimumConvergedPasses=2
PercentRefinement=30
IsEnabled=true
$begin 'Sweeps'
$begin 'Sweep'
RangeType='LinearCount'
RangeStart='8GHz'
RangeEnd='12GHz'
RangeCount=201
$end 'Sweep'
$end 'Sweeps'
$end 'Setup1'
$end 'SolveSetups'
$end 'AnalysisSetup'
$end 'HFSSModel'
$end 'AnsoftProject'
"""


def importer_for(text):
    imp = I._Importer.__new__(I._Importer)   # skip geometry: only setup() is under test
    root = parse(text, is_text=True)
    imp.d = root.find("HFSSModel") if root.name != "HFSSModel" else root
    imp.report = I.ImportReport("test.aedt", "HFSSModel")
    imp.units, imp.var_values = "mm", {}
    return imp


def test_mapping():
    s = importer_for(AEDT_SETUP).setup()
    assert s.adaptive is True
    assert (s.max_passes, s.max_delta_s, s.min_converged_passes) == (12, 0.01, 2), s
    assert (float(s.f_start), float(s.f_stop), s.n_points) == (8.0, 12.0, 201), s


def test_missing_keys_keep_defaults():
    txt = "\n".join(l for l in AEDT_SETUP.splitlines() if not l.startswith(("MaxDeltaS", "MaximumPasses", "MinimumConv")))
    s = importer_for(txt).setup()
    d = Setup()
    assert s.adaptive and (s.max_passes, s.max_delta_s, s.min_converged_passes) == \
        (d.max_passes, d.max_delta_s, d.min_converged_passes)


def test_old_project_files_load():
    p = Project.from_dict({"setup": {"name": "Old", "f_start": "1", "f_stop": "2"}})
    assert p.setup.adaptive is False and p.setup.edge_refinement is True


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn()
            print("ok", name)
