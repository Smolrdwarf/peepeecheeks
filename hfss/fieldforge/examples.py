"""Built-in example projects (File > Examples)."""
from __future__ import annotations

from .model import LumpedPort, Primitive, Project, Setup


def dipole() -> Project:
    p = Project(name="HalfWaveDipole", units="mm")
    p.variables = {"L": "70", "gap": "3", "r": "0.5"}
    p.objects = [
        Primitive("ArmTop", "cylinder", "copper", axis="Z",
                  params={"x": "0", "y": "0", "z": "gap/2", "radius": "r", "height": "L/2 - gap/2"}),
        Primitive("ArmBottom", "cylinder", "copper", axis="Z",
                  params={"x": "0", "y": "0", "z": "-L/2", "radius": "r", "height": "L/2 - gap/2"}),
    ]
    p.ports = [LumpedPort("Feed", axis="Z", impedance="73",
                          params={"x1": "0", "y1": "0", "z1": "-gap/2",
                                  "x2": "0", "y2": "0", "z2": "gap/2"})]
    p.setup = Setup(f_start="1", f_stop="3", farfield=True, farfield_freq="2")
    return p


def patch_antenna() -> Project:
    p = Project(name="PatchAntenna", units="mm")
    p.variables = {
        "patch_W": "32", "patch_L": "40", "sub_W": "60", "sub_L": "60",
        "h": "1.524", "feed_x": "-6",
    }
    p.objects = [
        Primitive("Substrate", "box", "Rogers RO4003C", transparency=0.7,
                  params={"x": "-sub_W/2", "y": "-sub_L/2", "z": "0",
                          "dx": "sub_W", "dy": "sub_L", "dz": "h"}),
        Primitive("Ground", "rectangle", "pec", axis="Z",
                  params={"x": "-sub_W/2", "y": "-sub_L/2", "z": "0", "w": "sub_W", "l": "sub_L"}),
        Primitive("Patch", "rectangle", "pec", axis="Z",
                  params={"x": "-patch_W/2", "y": "-patch_L/2", "z": "h", "w": "patch_W", "l": "patch_L"}),
    ]
    p.ports = [LumpedPort("Probe", axis="Z", impedance="50",
                          params={"x1": "feed_x", "y1": "0", "z1": "0",
                                  "x2": "feed_x", "y2": "0", "z2": "h"})]
    p.setup = Setup(f_start="1", f_stop="3.5", farfield=True, farfield_freq="2.3",
                    end_criteria_db=-40)
    return p


EXAMPLES = {
    "Half-wave dipole (2 GHz)": dipole,
    "Probe-fed patch antenna (RO4003C)": patch_antenna,
}
