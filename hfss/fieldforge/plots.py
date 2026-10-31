"""Result viewers: rectangular S-parameter plots, Smith chart, impedance, far-field patterns."""
from __future__ import annotations

import numpy as np
import pyvista as pv
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.figure import Figure
from PySide6.QtWidgets import QVBoxLayout, QWidget
from pyvistaqt import QtInteractor

from .viewer import primitive_mesh

TRACE_COLORS = ["#d62728", "#1f77b4", "#2ca02c", "#9467bd", "#ff7f0e", "#8c564b"]


class Results:
    """Thin wrapper around a results.npz file."""

    def __init__(self, path: str):
        d = np.load(path, allow_pickle=False)
        self.path = path
        self.freq = d["freq"]
        self.S = d["S"]
        self.Zin = d["Zin"]
        self.excited = int(d["excited"])
        self.port_names = [str(n) for n in d["port_names"]]
        self.has_farfield = "ff_D" in d
        if self.has_farfield:
            self.ff_freq = float(d["ff_freq"])
            self.ff_theta = d["ff_theta"]
            self.ff_phi = d["ff_phi"]
            self.ff_D = d["ff_D"]
            self.ff_Dmax = float(d["ff_Dmax"])
            self.ff_Prad = float(d["ff_Prad"])
            self.ff_Pin = float(d["ff_Pin"])

    def s_label(self, i: int) -> str:
        return f"S({self.port_names[i]},{self.port_names[self.excited]})"


class _MplPanel(QWidget):
    def __init__(self, parent=None, polar=False):
        super().__init__(parent)
        self.figure = Figure(figsize=(6, 4), dpi=100, layout="constrained")
        self.canvas = FigureCanvasQTAgg(self.figure)
        self.ax = self.figure.add_subplot(111, projection="polar" if polar else None)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(NavigationToolbar2QT(self.canvas, self))
        layout.addWidget(self.canvas)

    def _style(self, title, xlabel=None, ylabel=None):
        self.ax.set_title(title, fontsize=11, fontweight="bold", loc="left")
        if xlabel:
            self.ax.set_xlabel(xlabel)
        if ylabel:
            self.ax.set_ylabel(ylabel)
        self.ax.grid(True, which="both", color="#c8c8c8", linewidth=0.6)
        self.canvas.draw_idle()


class SParameterPlot(_MplPanel):
    def __init__(self, res: Results, title="S Parameter Plot", parent=None):
        super().__init__(parent)
        f = res.freq / 1e9
        for i in range(res.S.shape[0]):
            db = 20 * np.log10(np.maximum(np.abs(res.S[i]), 1e-12))
            self.ax.plot(f, db, color=TRACE_COLORS[i % len(TRACE_COLORS)], lw=1.8,
                         label=f"dB({res.s_label(i)})")
            if i == res.excited:
                k = int(np.argmin(db))
                self.ax.plot(f[k], db[k], "v", color="black")
                self.ax.annotate(f"m1: {f[k]:.4f} GHz, {db[k]:.2f} dB", (f[k], db[k]),
                                 xytext=(12, 18), textcoords="offset points", fontsize=8,
                                 arrowprops=dict(arrowstyle="-", color="#666"),
                                 bbox=dict(boxstyle="round", fc="#ffffe0", ec="#888"))
        self.ax.set_xlim(f[0], f[-1])
        lo = min(line.get_ydata().min() for line in self.ax.get_lines())
        self.ax.set_ylim(min(lo - 3, -10), 1)
        self.ax.legend(loc="lower left", fontsize=8, framealpha=0.9)
        self._style(title, "Freq [GHz]", "dB")


class ImpedancePlot(_MplPanel):
    def __init__(self, res: Results, title="Input Impedance", parent=None):
        super().__init__(parent)
        f = res.freq / 1e9
        z = res.Zin[res.excited]
        name = res.port_names[res.excited]
        self.ax.plot(f, z.real, color=TRACE_COLORS[0], lw=1.8, label=f"re(Z({name}))")
        self.ax.plot(f, z.imag, color=TRACE_COLORS[1], lw=1.8, ls="--", label=f"im(Z({name}))")
        self.ax.axhline(0, color="black", lw=0.8)
        self.ax.set_xlim(f[0], f[-1])
        self.ax.legend(loc="best", fontsize=8)
        self._style(title, "Freq [GHz]", "Ohm")


class SmithChart(_MplPanel):
    def __init__(self, res: Results, title="Smith Chart", parent=None):
        super().__init__(parent)
        ax = self.ax
        ax.set_aspect("equal")
        ax.axis("off")
        t = np.linspace(0, 2 * np.pi, 400)
        ax.plot(np.cos(t), np.sin(t), color="#555", lw=1.2)
        grid = dict(color="#c0c0c0", lw=0.6)
        for r in (0.2, 0.5, 1, 2, 5):
            ax.plot(r / (1 + r) + np.cos(t) / (1 + r), np.sin(t) / (1 + r), **grid)
            ax.text(r / (1 + r) - 1 / (1 + r) + 0.01, 0.015, f"{r:g}", fontsize=7, color="#666")
        for x in (0.2, 0.5, 1, 2, 5):
            for sx in (x, -x):
                z = 1j * sx + np.logspace(-3, 3, 600)
                g = (z - 1) / (z + 1)
                ax.plot(g.real, g.imag, **grid)
        ax.plot([-1, 1], [0, 0], **grid)
        s = res.S[res.excited]
        ax.plot(s.real, s.imag, color=TRACE_COLORS[0], lw=1.8, label=res.s_label(res.excited))
        ax.plot(s.real[0], s.imag[0], "o", color=TRACE_COLORS[0], ms=4)
        ax.legend(loc="upper left", fontsize=8)
        ax.set_xlim(-1.08, 1.08)
        ax.set_ylim(-1.08, 1.08)
        self._style(title)


class FarField2D(_MplPanel):
    def __init__(self, res: Results, title="Radiation Pattern", parent=None):
        super().__init__(parent, polar=True)
        th = np.deg2rad(res.ff_theta)
        d_db = 10 * np.log10(np.maximum(res.ff_D, 1e-12))
        floor = d_db.max() - 30
        for k, phi_cut in enumerate((0.0, 90.0)):
            j = int(np.argmin(np.abs(res.ff_phi - phi_cut)))
            j2 = int(np.argmin(np.abs(res.ff_phi - (phi_cut + 180.0))))
            # full great-circle cut: phi and phi+180 halves
            ang = np.concatenate([th, 2 * np.pi - th[::-1]])
            val = np.concatenate([d_db[:, j], d_db[::-1, j2]])
            self.ax.plot(ang, np.maximum(val, floor), color=TRACE_COLORS[k], lw=1.8,
                         label=f"Phi = {phi_cut:g} deg")
        self.ax.set_theta_zero_location("N")
        self.ax.set_theta_direction(-1)
        self.ax.set_rlim(floor, d_db.max() + 1)
        self.ax.legend(loc="lower left", bbox_to_anchor=(-0.15, -0.12), fontsize=8)
        self._style(f"{title}  -  dB(Directivity) [dBi]  @ {res.ff_freq / 1e9:g} GHz")


class FarField3D(QWidget):
    def __init__(self, res: Results, project=None, parent=None):
        super().__init__(parent)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        self.plotter = QtInteractor(self)
        layout.addWidget(self.plotter)
        p = self.plotter
        p.set_background("#ffffff", top="#dfe7f2")
        p.add_axes(line_width=2, color="black")

        th, ph = np.meshgrid(np.deg2rad(res.ff_theta), np.deg2rad(res.ff_phi), indexing="ij")
        d_db = 10 * np.log10(np.maximum(res.ff_D, 1e-12))
        dyn = 30.0
        r = np.clip((d_db - (d_db.max() - dyn)) / dyn, 0, 1)

        scale = 1.0
        variables = None
        if project is not None and project.objects:
            try:
                variables = project.resolved_variables()
                from .solver import model_bounds
                lo, hi = model_bounds(project, variables)
                scale = float(np.linalg.norm(hi - lo)) * 0.9
            except Exception:  # noqa: BLE001
                variables = None
        x = scale * r * np.sin(th) * np.cos(ph)
        y = scale * r * np.sin(th) * np.sin(ph)
        z = scale * r * np.cos(th)
        grid = pv.StructuredGrid(x, y, z)
        grid["dB(DirTotal)"] = d_db.ravel(order="F")
        p.add_mesh(grid, scalars="dB(DirTotal)", cmap="jet", opacity=0.85, smooth_shading=True,
                   clim=(d_db.max() - dyn, d_db.max()),
                   scalar_bar_args=dict(title="dB(DirTotal) [dBi]", vertical=True, fmt="%.1f",
                                        position_x=0.86, position_y=0.25, height=0.55, width=0.06,
                                        color="black"))
        if variables is not None:
            for obj in project.objects:
                if obj.visible:
                    try:
                        p.add_mesh(primitive_mesh(obj, variables), color=project.material(obj.material).color,
                                   opacity=0.5)
                    except Exception:  # noqa: BLE001
                        pass
        eff = 100 * res.ff_Prad / res.ff_Pin if res.ff_Pin > 0 else float("nan")
        p.add_text(f"3D Radiation Pattern @ {res.ff_freq / 1e9:g} GHz\n"
                   f"Peak directivity: {10 * np.log10(res.ff_Dmax):.2f} dBi\n"
                   f"Radiation efficiency: {eff:.1f} %",
                   position="upper_left", font_size=8, color="black")
        p.view_isometric()
        p.reset_camera()
