"""Main window: project manager, properties, 3D modeler, message manager and solver progress."""
from __future__ import annotations

import copy
import os
import re
import sys
import tempfile
import time

import numpy as np
from PySide6.QtCore import QProcess, QProcessEnvironment, Qt
from PySide6.QtGui import QAction, QColor, QKeySequence
from PySide6.QtWidgets import (
    QAbstractItemView, QApplication, QCheckBox, QColorDialog, QComboBox, QDockWidget, QFileDialog,
    QHBoxLayout, QHeaderView, QInputDialog, QLabel, QListWidget, QListWidgetItem, QMainWindow,
    QMenu, QMessageBox, QPlainTextEdit, QProgressBar, QPushButton, QSplitter, QStyle, QTableWidget,
    QTableWidgetItem, QTabWidget, QTextBrowser, QTreeWidget, QTreeWidgetItem, QVBoxLayout, QWidget,
    QDialog, QDialogButtonBox,
)

from . import examples
from .icons import color_icon, icon
from .model import (
    AXES, AXIS_KINDS, BOUNDARY_TYPES, FACE_NAMES, PRIMITIVE_PARAMS, UNITS, LumpedPort, Material,
    Primitive, Project, evaluate,
)
from . import reference as refmod
from .plots import ComparisonPlot, FarField2D, FarField3D, ImpedancePlot, Results, SParameterPlot, SmithChart
from .viewer import ModelerView

APP_NAME = "FieldForge"
FILE_FILTER = "FieldForge project (*.ffproj);;JSON (*.json)"
AEDT_FILTER = "Ansys Electronics Desktop project (*.aedt *.aedtz)"
REF_FILTER = "Reference results (*.s1p *.s2p *.s3p *.s4p *.snp *.csv);;Touchstone (*.s*p);;HFSS report CSV (*.csv)"
HFSS_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
IS_MAC = sys.platform == "darwin"

ROLE = Qt.UserRole

TS_RE = re.compile(r"Timestep:\s*(\d+)")
ENERGY_RE = re.compile(r"Energy:.*\(\s*([-\d.]+)\s*dB\)")


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.project = Project()
        self.path: str | None = None
        self.dirty = False
        self.results: Results | None = None
        self.undo_stack: list[dict] = []
        self.redo_stack: list[dict] = []
        self.process: QProcess | None = None
        self.references: list = []          # loaded reference results (refmod.Reference)
        self.import_report = None
        self.selected: tuple | None = None
        self._loading_props = False

        self.resize(1500, 920)
        self._build_central()
        self._build_docks()
        self._build_actions()
        self._build_menus()
        self._build_toolbars()
        self.statusBar().showMessage("Ready")
        self._units_label = QLabel()
        self.statusBar().addPermanentWidget(self._units_label)

        self.load_project(examples.patch_antenna(), None)
        self.message("info", f"{APP_NAME} started. openEMS FDTD solver backend.")

    # ================================================================== layout
    def _build_central(self):
        self.tabs = QTabWidget()
        self.tabs.setTabsClosable(True)
        self.tabs.setDocumentMode(True)
        self.tabs.tabCloseRequested.connect(self._close_tab)
        self.view = ModelerView(self.tabs)
        self.view.picked.connect(self._on_pick)
        self.tabs.addTab(self.view, icon("box"), "3D Modeler")
        self.tabs.tabBar().setTabButton(0, self.tabs.tabBar().ButtonPosition.RightSide, None)
        self.setCentralWidget(self.tabs)

    def _build_docks(self):
        # Project manager + properties
        self.tree = QTreeWidget()
        self.tree.setHeaderHidden(True)
        self.tree.setSelectionMode(QAbstractItemView.SingleSelection)
        self.tree.itemSelectionChanged.connect(self._on_tree_select)
        self.tree.itemDoubleClicked.connect(self._on_tree_double)
        self.tree.setContextMenuPolicy(Qt.CustomContextMenu)
        self.tree.customContextMenuRequested.connect(self._tree_menu)

        self.props = QTableWidget(0, 2)
        self.props.setHorizontalHeaderLabels(["Name", "Value"])
        self.props.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeToContents)
        self.props.horizontalHeader().setSectionResizeMode(1, QHeaderView.Stretch)
        self.props.verticalHeader().setVisible(False)
        self.props.setAlternatingRowColors(True)
        self.props.itemChanged.connect(self._on_prop_edited)
        self.props_title = QLabel(" Properties")
        self.props_title.setStyleSheet("font-weight: bold; padding: 3px;")
        pw = QWidget()
        pl = QVBoxLayout(pw)
        pl.setContentsMargins(0, 0, 0, 0)
        pl.setSpacing(0)
        pl.addWidget(self.props_title)
        pl.addWidget(self.props)

        split = QSplitter(Qt.Vertical)
        split.addWidget(self.tree)
        split.addWidget(pw)
        split.setSizes([480, 360])
        dock = QDockWidget("Project Manager", self)
        dock.setObjectName("project_manager")
        dock.setWidget(split)
        dock.setMinimumWidth(330)
        self.addDockWidget(Qt.LeftDockWidgetArea, dock)
        self.dock_project = dock

        # Message manager
        self.messages = QListWidget()
        self.messages.setWordWrap(True)
        dm = QDockWidget("Message Manager", self)
        dm.setObjectName("messages")
        dm.setWidget(self.messages)
        self.addDockWidget(Qt.BottomDockWidgetArea, dm)
        self.dock_messages = dm

        # Progress
        w = QWidget()
        lay = QVBoxLayout(w)
        lay.setContentsMargins(4, 4, 4, 4)
        row = QHBoxLayout()
        self.progress_label = QLabel("Idle")
        self.progress = QProgressBar()
        self.progress.setRange(0, 1000)
        self.progress.setTextVisible(False)
        self.abort_btn = QPushButton(icon("abort"), "Abort")
        self.abort_btn.setEnabled(False)
        self.abort_btn.clicked.connect(self.abort)
        row.addWidget(self.progress_label, 1)
        row.addWidget(self.progress, 2)
        row.addWidget(self.abort_btn)
        lay.addLayout(row)
        self.log = QPlainTextEdit()
        self.log.setReadOnly(True)
        self.log.setMaximumBlockCount(4000)
        self.log.setStyleSheet(f"font-family: {'Menlo' if IS_MAC else 'Consolas'}, monospace; font-size: 9pt;")
        lay.addWidget(self.log)
        dp = QDockWidget("Progress", self)
        dp.setObjectName("progress")
        dp.setWidget(w)
        self.addDockWidget(Qt.BottomDockWidgetArea, dp)
        self.splitDockWidget(dm, dp, Qt.Horizontal)
        self.resizeDocks([dm], [230], Qt.Vertical)
        self.dock_progress = dp

    def _act(self, text, slot, shortcut=None, ico=None, tip=None):
        a = QAction(text, self)
        if ico is not None:
            a.setIcon(icon(ico) if isinstance(ico, str) else ico)
        if shortcut:
            a.setShortcut(QKeySequence(shortcut))
        a.setStatusTip(tip or text.replace("&", ""))
        a.triggered.connect(slot)
        return a

    def _build_actions(self):
        st = self.style()
        A = self._act
        self.a_new = A("&New", self.new_project, "Ctrl+N", st.standardIcon(QStyle.SP_FileIcon))
        self.a_open = A("&Open...", self.open_project, "Ctrl+O", st.standardIcon(QStyle.SP_DialogOpenButton))
        self.a_save = A("&Save", self.save_project, "Ctrl+S", st.standardIcon(QStyle.SP_DialogSaveButton))
        self.a_saveas = A("Save &As...", self.save_project_as, "Ctrl+Shift+S")
        self.a_exit = A("E&xit", self.close, QKeySequence.Quit if IS_MAC else "Alt+F4")
        self.a_import_aedt = A("&Import Ansys HFSS Design (.aedt)...", self.import_aedt_dialog, "Ctrl+I", "import",
                               "Import geometry, materials, ports, boundaries and setup from an HFSS project")
        self.a_report = A("Import &Report", self.show_import_report)
        self.a_ref = A("Load &Reference Results (Touchstone / HFSS CSV)...", self.load_reference, None, "reference",
                       "Load S-parameters exported from HFSS to compare against")
        self.a_compare = A("&Compare with Reference", lambda: self.open_result("compare"), None, "compare")
        self.a_undo = A("&Undo", self.undo, "Ctrl+Z", st.standardIcon(QStyle.SP_ArrowBack))
        self.a_redo = A("&Redo", self.redo, QKeySequence.Redo, st.standardIcon(QStyle.SP_ArrowForward))
        # The key labelled "delete" on a Mac keyboard is Backspace to Qt.
        self.a_delete = A("&Delete", self.delete_selected, "Backspace" if IS_MAC else "Delete",
                          st.standardIcon(QStyle.SP_TrashIcon))
        self.a_dup = A("D&uplicate", self.duplicate_selected, "Ctrl+D")
        self.a_box = A("&Box", lambda: self.draw("box"), None, "box", "Draw a box")
        self.a_cyl = A("&Cylinder", lambda: self.draw("cylinder"), None, "cylinder", "Draw a cylinder")
        self.a_sph = A("&Sphere", lambda: self.draw("sphere"), None, "sphere", "Draw a sphere")
        self.a_rect = A("&Rectangle", lambda: self.draw("rectangle"), None, "rectangle", "Draw a 2D sheet")
        self.a_port = A("&Lumped Port", self.add_port, None, "port", "Assign a lumped port excitation")
        self.a_region = A("&Region / Boundaries", lambda: self.select_key(("region",)), None, "region")
        self.a_setup = A("Analysis &Setup", lambda: self.select_key(("setup",)), None, "setup")
        self.a_validate = A("&Validation Check", self.validate, None, "validate")
        self.a_analyze = A("&Analyze All", self.analyze, "F5", "analyze")
        self.a_abort = A("A&bort", self.abort, None, "abort")
        self.a_var = A("Design &Variables...", self.add_variable, None, "variable")
        self.a_units = A("&Units...", self.set_units)
        self.a_mat = A("New &Material...", self.add_material)
        self.a_s = A("&S Parameter Plot", lambda: self.open_result("sparam"), None, "plot")
        self.a_smith = A("S&mith Chart", lambda: self.open_result("smith"), None, "smith")
        self.a_z = A("&Input Impedance", lambda: self.open_result("zin"), None, "plot")
        self.a_ff3 = A("&3D Far Field", lambda: self.open_result("ff3d"), None, "farfield")
        self.a_ff2 = A("&2D Radiation Pattern", lambda: self.open_result("ff2d"), None, "farfield")
        self.a_touch = A("Export &Touchstone...", self.export_touchstone)
        self.a_fit = A("&Fit All", self.view.fit_all, "Ctrl+Shift+F", "fit")
        self.a_iso = A("&Isometric", lambda: (self.view.view_isometric(), self.view.reset_camera()), "Ctrl+0")
        self.a_top = A("&Top (XY)", self.view.view_xy, "Ctrl+1")
        self.a_front = A("F&ront (XZ)", self.view.view_xz, "Ctrl+2")
        self.a_right = A("&Right (YZ)", self.view.view_yz, "Ctrl+3")
        self.a_grid = A("Show &Grid", lambda c: (self.view.set_grid_visible(c), self.refresh_view()))
        self.a_grid.setCheckable(True)
        self.a_grid.setChecked(True)
        self.a_about = A("&About", self.about)

    def _build_menus(self):
        mb = self.menuBar()
        m = mb.addMenu("&File")
        for a in (self.a_new, self.a_open, None, self.a_import_aedt):
            m.addSeparator() if a is None else m.addAction(a)
        ex = m.addMenu("&Examples")
        for title, factory in examples.EXAMPLES.items():
            ex.addAction(title, lambda f=factory: self._open_example(f))
        m.addSeparator()
        for a in (self.a_save, self.a_saveas):
            m.addAction(a)
        m.addSeparator()
        m.addAction(self.a_exit)
        m = mb.addMenu("&Edit")
        for a in (self.a_undo, self.a_redo, None, self.a_delete, self.a_dup):
            m.addSeparator() if a is None else m.addAction(a)
        m = mb.addMenu("&View")
        for a in (self.a_fit, self.a_iso, self.a_top, self.a_front, self.a_right, None, self.a_grid):
            m.addSeparator() if a is None else m.addAction(a)
        m.addSeparator()
        for d in (self.dock_project, self.dock_messages, self.dock_progress):
            m.addAction(d.toggleViewAction())
        m = mb.addMenu("&Draw")
        for a in (self.a_box, self.a_cyl, self.a_sph, None, self.a_rect):
            m.addSeparator() if a is None else m.addAction(a)
        m = mb.addMenu("&Modeler")
        for a in (self.a_units, self.a_var, self.a_mat):
            m.addAction(a)
        m = mb.addMenu("&Simulation")
        sub = m.addMenu("&Excitations")
        sub.addAction(self.a_port)
        sub = m.addMenu("&Boundaries")
        sub.addAction(self.a_region)
        m.addAction(self.a_setup)
        m.addSeparator()
        for a in (self.a_validate, self.a_analyze, self.a_abort):
            m.addAction(a)
        m = mb.addMenu("&Results")
        for a in (self.a_s, self.a_smith, self.a_z, None, self.a_ff3, self.a_ff2, None, self.a_ref,
                  self.a_compare, None, self.a_touch):
            m.addSeparator() if a is None else m.addAction(a)
        m = mb.addMenu("&Help")
        m.addAction(self.a_about)

    def _build_toolbars(self):
        def tb(name, actions):
            t = self.addToolBar(name)
            t.setObjectName(name)
            t.setIconSize(t.iconSize() * 0.9)
            for a in actions:
                t.addSeparator() if a is None else t.addAction(a)
            return t
        tb("File", [self.a_new, self.a_open, self.a_save, self.a_import_aedt, None, self.a_undo, self.a_redo,
                    self.a_delete])
        tb("Draw", [self.a_box, self.a_cyl, self.a_sph, self.a_rect])
        tb("View", [self.a_fit])
        tb("Simulation", [self.a_port, self.a_region, self.a_setup, self.a_var, None,
                          self.a_validate, self.a_analyze, self.a_abort])
        tb("Results", [self.a_s, self.a_smith, self.a_ff3, None, self.a_ref, self.a_compare])

    # ================================================================== messages
    def message(self, level: str, text: str):
        st = self.style()
        ico = {"info": QStyle.SP_MessageBoxInformation, "warn": QStyle.SP_MessageBoxWarning,
               "error": QStyle.SP_MessageBoxCritical}[level]
        item = QListWidgetItem(st.standardIcon(ico), f"[{time.strftime('%H:%M:%S')}] {text}")
        self.messages.addItem(item)
        self.messages.scrollToBottom()

    # ================================================================== project state
    def load_project(self, project: Project, path: str | None):
        self.project = project
        self.path = path
        self.dirty = False
        self.undo_stack.clear()
        self.redo_stack.clear()
        self.results = None
        self.selected = None
        for i in range(self.tabs.count() - 1, 0, -1):
            self.tabs.removeTab(i)
        self._try_load_results()
        self._load_references()
        self.refresh_all(reset_camera=True)

    def _load_references(self):
        self.references = []
        for r in self.project.references:
            try:
                self.references.append(refmod.load(r["path"], r.get("label")))
            except Exception as exc:  # noqa: BLE001
                self.message("warn", f"Reference '{r.get('path')}' could not be loaded: {exc}")

    def _results_dir(self) -> str:
        if self.path:
            base = os.path.splitext(self.path)[0] + ".ffresults"
        else:
            base = os.path.join(tempfile.gettempdir(), "fieldforge", self.project.name)
        return os.path.join(base, self.project.setup.name)

    def _try_load_results(self):
        f = os.path.join(self._results_dir(), "results.npz")
        if os.path.exists(f):
            try:
                self.results = Results(f)
                self.message("info", f"Loaded existing solution from {f}")
            except Exception as exc:  # noqa: BLE001
                self.message("warn", f"Could not load results: {exc}")

    def _snapshot(self):
        self.undo_stack.append(copy.deepcopy(self.project.to_dict()))
        if len(self.undo_stack) > 200:
            self.undo_stack.pop(0)
        self.redo_stack.clear()
        self.dirty = True

    def undo(self):
        if not self.undo_stack:
            return
        self.redo_stack.append(self.project.to_dict())
        self.project = Project.from_dict(self.undo_stack.pop())
        self.dirty = True
        self._load_references()
        self.refresh_all()

    def redo(self):
        if not self.redo_stack:
            return
        self.undo_stack.append(self.project.to_dict())
        self.project = Project.from_dict(self.redo_stack.pop())
        self.dirty = True
        self._load_references()
        self.refresh_all()

    def refresh_all(self, reset_camera=False):
        self.refresh_tree()
        self.refresh_view(reset_camera)
        self.show_properties(self.selected)
        self._update_title()

    def refresh_view(self, reset_camera=False):
        errors: list[str] = []
        self.view.rebuild(self.project, errors, keep_camera=not reset_camera)
        for e in errors:
            self.message("error", e)
        if self.selected and self.selected[0] in ("object", "port"):
            self.view.highlight(self.selected[1])

    def _update_title(self):
        name = os.path.basename(self.path) if self.path else self.project.name
        self.setWindowTitle(f"{APP_NAME} - {name}{' *' if self.dirty else ''}")
        self._units_label.setText(f"  Units: {self.project.units}  ")

    # ================================================================== project tree
    def _item(self, parent, text, key, ico=None, bold=False):
        it = QTreeWidgetItem(parent, [text])
        it.setData(0, ROLE, key)
        if ico is not None:
            it.setIcon(0, icon(ico) if isinstance(ico, str) else ico)
        if bold:
            f = it.font(0)
            f.setBold(True)
            it.setFont(0, f)
        return it

    def refresh_tree(self):
        p = self.project
        self.tree.blockSignals(True)
        expanded = {self._key_str(self.tree.topLevelItem(0))} if self.tree.topLevelItemCount() else set()
        expanded |= {self._key_str(it) for it in self._iter_items() if it.isExpanded()}
        self.tree.clear()
        st = self.style()
        folder = st.standardIcon(QStyle.SP_DirIcon)
        root = self._item(self.tree, p.name, ("project",), st.standardIcon(QStyle.SP_DriveHDIcon), True)
        design = self._item(root, f"{p.name}Design (DrivenModal, FDTD)", ("design",),
                            st.standardIcon(QStyle.SP_ComputerIcon))
        model = self._item(design, "Model", ("model",), icon("box"))
        by_mat: dict[str, list] = {}
        for o in p.objects:
            by_mat.setdefault(o.material, []).append(o)
        for mat_name, objs in by_mat.items():
            mat = p.material(mat_name)
            grp = self._item(model, mat_name, ("matgroup", mat_name), color_icon(mat.color))
            for o in objs:
                self._item(grp, o.name, ("object", o.name), "box" if o.kind == "cad" else o.kind)
        b = self._item(design, "Boundaries", ("boundaries",), folder)
        bset = sorted(set(p.region.boundaries))
        self._item(b, f"Region ({', '.join(bset)})", ("region",), "region")
        e = self._item(design, "Excitations", ("excitations",), folder)
        for port in p.ports:
            self._item(e, port.name + ("" if port.excited else " (passive)"), ("port", port.name), "port")
        a = self._item(design, "Analysis", ("analysis",), folder)
        self._item(a, p.setup.name, ("setup",), "setup")
        r = self._item(design, "Results", ("results",), folder)
        if self.results is not None:
            self._item(r, "S Parameter Plot", ("result", "sparam"), "plot")
            self._item(r, "Smith Chart", ("result", "smith"), "smith")
            self._item(r, "Input Impedance", ("result", "zin"), "plot")
            if self.results.has_farfield:
                self._item(r, "3D Far Field", ("result", "ff3d"), "farfield")
                self._item(r, "2D Radiation Pattern", ("result", "ff2d"), "farfield")
            if self.references:
                self._item(r, "Comparison with Reference", ("result", "compare"), "compare")
        else:
            it = self._item(r, "(no solution data)", ("noresult",))
            it.setForeground(0, QColor("#888"))
            if self.references:
                self._item(r, "S Parameter Plot (reference only)", ("result", "sparam"), "plot")
        for k, ref in enumerate(self.references):
            self._item(r, f"Reference: {ref.label}", ("reference", k), "reference")
        if p.source:
            self._item(root, f"Source: {os.path.basename(p.source.get('path', ''))} / {p.source.get('design')}",
                       ("source",), "import")
        d = self._item(root, "Definitions", ("definitions",), folder)
        self._item(d, f"Design Variables ({len(p.variables)})", ("variables",), "variable")
        mats = self._item(d, "Materials", ("materials",), folder)
        for name, mat in p.all_materials().items():
            label = name + ("  [project]" if name in p.materials else "")
            self._item(mats, label, ("material", name), color_icon(mat.color))

        for it in self._iter_items():
            k = self._key_str(it)
            if k in expanded or it.data(0, ROLE)[0] in ("project", "design", "model", "matgroup",
                                                        "boundaries", "excitations", "analysis", "results"):
                it.setExpanded(True)
        self.tree.blockSignals(False)
        if self.selected:
            self._select_tree_key(self.selected)

    def _iter_items(self):
        stack = [self.tree.topLevelItem(i) for i in range(self.tree.topLevelItemCount())]
        while stack:
            it = stack.pop()
            yield it
            stack.extend(it.child(i) for i in range(it.childCount()))

    @staticmethod
    def _key_str(it):
        return repr(it.data(0, ROLE)) if it is not None else ""

    def _select_tree_key(self, key):
        self.tree.blockSignals(True)
        for it in self._iter_items():
            if it.data(0, ROLE) == key:
                self.tree.setCurrentItem(it)
                break
        else:
            self.tree.clearSelection()
        self.tree.blockSignals(False)

    def select_key(self, key):
        self.selected = key
        self._select_tree_key(key)
        self.show_properties(key)
        self.view.highlight(key[1] if key and key[0] in ("object", "port") else None)

    def _on_tree_select(self):
        it = self.tree.currentItem()
        key = it.data(0, ROLE) if it is not None else None
        self.selected = key
        self.show_properties(key)
        self.view.highlight(key[1] if key and key[0] in ("object", "port") else None)

    def _on_tree_double(self, it, _col):
        key = it.data(0, ROLE)
        if key and key[0] == "result":
            self.open_result(key[1])
        elif key and key[0] == "reference":
            self.open_result("compare" if self.results is not None else "sparam")
        elif key and key[0] == "source":
            self.show_import_report()

    def _on_pick(self, name):
        if not name:
            self.select_key(None)
            return
        kind = "port" if any(p.name == name for p in self.project.ports) else "object"
        self.select_key((kind, name))

    def _tree_menu(self, pos):
        it = self.tree.itemAt(pos)
        key = it.data(0, ROLE) if it else None
        m = QMenu(self)
        kind = key[0] if key else None
        if kind in ("object", "port"):
            m.addAction(self.a_dup)
            m.addAction(self.a_delete)
        elif kind == "model":
            for a in (self.a_box, self.a_cyl, self.a_sph, self.a_rect):
                m.addAction(a)
        elif kind == "excitations":
            m.addAction(self.a_port)
        elif kind in ("analysis", "setup"):
            m.addAction(self.a_validate)
            m.addAction(self.a_analyze)
        elif kind == "variables":
            m.addAction("Add Variable...", self.add_variable)
            if self.project.variables:
                m.addAction("Delete Variable...", self.delete_variable)
        elif kind in ("materials", "material"):
            m.addAction(self.a_mat)
            if kind == "material" and key[1] in self.project.materials:
                m.addAction(self.a_delete)
        elif kind == "results" or kind == "result":
            for a in (self.a_s, self.a_smith, self.a_z, self.a_ff3, self.a_ff2, self.a_touch, None, self.a_ref,
                      self.a_compare):
                m.addSeparator() if a is None else m.addAction(a)
        elif kind == "reference":
            m.addAction(self.a_compare)
            m.addAction("Remove Reference", lambda k=key[1]: self.remove_reference(k))
        elif kind in ("source", "project") and self.project.source:
            m.addAction(self.a_report)
        if not m.isEmpty():
            m.exec(self.tree.viewport().mapToGlobal(pos))

    # ================================================================== properties
    def _find(self, kind, name):
        seq = self.project.objects if kind == "object" else self.project.ports
        return next((o for o in seq if o.name == name), None)

    def show_properties(self, key):
        self._loading_props = True
        self.props.clearContents()
        self.props.setRowCount(0)
        self._prop_setters = []
        p = self.project
        title = "Properties"
        rows = []  # (label, value, kind, setter, options)

        def expr_setter(target: dict, k):
            def s(text):
                evaluate(text, p.resolved_variables())  # validate before committing
                target[k] = text
            return s

        def attr_setter(obj, attr, conv=str):
            def s(val):
                setattr(obj, attr, conv(val))
            return s

        kind = key[0] if key else None
        if kind == "object":
            o = self._find("object", key[1])
            if o is not None:
                title = f"Properties: {o.name} ({o.kind.capitalize()})"
                rows.append(("Name", o.name, "text", self._renamer(o), None))
                rows.append(("Material", o.material, "combo", attr_setter(o, "material"),
                             list(p.all_materials())))
                if o.kind in AXIS_KINDS:
                    rows.append(("Axis" if o.kind == "cylinder" else "Normal", o.axis, "combo",
                                 attr_setter(o, "axis"), list(AXES)))
                for k, label, _ in PRIMITIVE_PARAMS[o.kind]:
                    rows.append((f"{label} ({p.units})", o.params[k], "text", expr_setter(o.params, k), None))
                rows.append(("Transparency", f"{o.transparency:g}", "text",
                             attr_setter(o, "transparency", lambda v: min(max(float(v), 0.0), 0.95)), None))
                rows.append(("Visible", o.visible, "check", attr_setter(o, "visible", bool), None))
                if o.kind == "cad":
                    title = f"Properties: {o.name} (imported)"
                    rows.append(("Source part", o.part, "readonly", None, None))
                    try:
                        from . import cad
                        rows.append(("Geometry", cad.kind(p, o), "readonly", None, None))
                    except Exception as exc:  # noqa: BLE001
                        rows.append(("Geometry", f"unavailable: {exc}", "readonly", None, None))
                    rows.append(("Dimensions", "driven by design variables (Definitions)", "readonly", None, None))
        elif kind == "port":
            port = self._find("port", key[1])
            if port is not None:
                title = f"Properties: {port.name} (Lumped Port)"
                rows.append(("Name", port.name, "text", self._renamer(port), None))
                rows.append(("Excited", port.excited, "check", attr_setter(port, "excited", bool), None))
                rows.append(("Impedance (Ohm)", port.impedance, "text",
                             expr_setter(port.__dict__, "impedance"), None))
                rows.append(("Current direction", port.axis, "combo", attr_setter(port, "axis"), list(AXES)))
                for k in ("x1", "y1", "z1", "x2", "y2", "z2"):
                    label = {"1": "Start", "2": "End"}[k[1]] + f" {k[0].upper()} ({p.units})"
                    rows.append((label, port.params[k], "text", expr_setter(port.params, k), None))
        elif kind == "region":
            title = "Properties: Region / Boundaries"
            rows.append((f"Padding ({p.units}, or 'auto')", p.region.padding, "text",
                         self._padding_setter(), None))
            for i, face in enumerate(FACE_NAMES):
                def s(v, i=i):
                    p.region.boundaries[i] = v
                rows.append((f"Boundary {face}", p.region.boundaries[i], "combo", s, BOUNDARY_TYPES))
        elif kind == "setup":
            s = p.setup
            title = f"Properties: {s.name} (Driven, broadband FDTD)"

            def num(attr, conv=float):
                def f(v):
                    conv(v)
                    setattr(s, attr, v if attr in ("f_start", "f_stop", "farfield_freq") else conv(v))
                return f
            rows += [
                ("Start frequency (GHz)", s.f_start, "text", num("f_start"), None),
                ("Stop frequency (GHz)", s.f_stop, "text", num("f_stop"), None),
                ("Frequency points", str(s.n_points), "text", num("n_points", int), None),
                ("Mesh: cells per wavelength", str(s.cells_per_wavelength), "text",
                 num("cells_per_wavelength", int), None),
                ("Convergence: energy decay (dB)", f"{s.end_criteria_db:g}", "text", num("end_criteria_db"), None),
                ("Max time steps", str(s.max_timesteps), "text", num("max_timesteps", int), None),
                ("Compute far field", s.farfield, "check", attr_setter(s, "farfield", bool), None),
                ("Far-field frequency (GHz)", s.farfield_freq, "text", num("farfield_freq"), None),
                ("Far-field angular step (deg)", f"{s.farfield_step_deg:g}", "text",
                 num("farfield_step_deg"), None),
            ]
        elif kind == "variables":
            title = "Design Variables"
            vals = {}
            try:
                vals = p.resolved_variables()
            except Exception:  # noqa: BLE001
                pass
            for name, expr in p.variables.items():
                def s(v, name=name):
                    trial = dict(p.variables)
                    trial[name] = v
                    Project(variables=trial).resolved_variables()
                    p.variables[name] = v
                shown = f"{expr}" if name not in vals or str(vals[name]) == expr else expr
                rows.append((f"{name}  (= {vals.get(name, float('nan')):g})", shown, "text", s, None))
        elif kind == "material":
            mat = p.all_materials()[key[1]]
            editable = key[1] in p.materials
            title = f"Material: {mat.name}" + ("" if editable else "  (library, read-only)")
            ro = "text" if editable else "readonly"
            rows += [
                ("Relative permittivity", f"{mat.epsr:g}", ro, attr_setter(mat, "epsr", float), None),
                ("Relative permeability", f"{mat.mur:g}", ro, attr_setter(mat, "mur", float), None),
                ("Dielectric loss tangent", f"{mat.tand:g}", ro, attr_setter(mat, "tand", float), None),
                ("Perfect conductor", mat.conductor, "check" if editable else "readonly",
                 attr_setter(mat, "conductor", bool), None),
            ]
            if editable:
                rows.append(("Color", "Pick...", "color", self._color_setter(mat), None))
        elif kind == "project":
            title = "Project"
            rows.append(("Name", p.name, "text", attr_setter(p, "name"), None))
            rows.append(("Units", p.units, "combo", self._units_setter(), list(UNITS)))
            rows.append(("File", self.path or "(unsaved)", "readonly", None, None))
            if p.source:
                rows.append(("Imported from", p.source.get("path", ""), "readonly", None, None))
                rows.append(("HFSS design", p.source.get("design", ""), "readonly", None, None))
        elif kind == "source":
            title = "Imported design"
            rows.append(("File", p.source.get("path", ""), "readonly", None, None))
            rows.append(("Design", p.source.get("design", ""), "readonly", None, None))
            rows.append(("Model units", p.source.get("model_units", ""), "readonly", None, None))
            rows.append(("Geometry", "rebuilt from the HFSS modeler history; edit design variables to change it",
                         "readonly", None, None))
        elif kind == "reference" and key[1] < len(self.references):
            ref = self.references[key[1]]
            title = f"Reference: {ref.label}"
            rows.append(("Label", ref.label, "text", self._ref_label_setter(key[1]), None))
            rows.append(("File", ref.path, "readonly", None, None))
            rows.append(("Traces", ", ".join(t.label for t in ref.traces[:8]), "readonly", None, None))
            tr = ref.reflection()
            if tr is not None:
                m = refmod.metrics(tr.freq, tr.db)
                rows.append(("Resonance", f"{m.f_res / 1e9:.4f} GHz, {m.depth_db:.2f} dB", "readonly", None, None))
        elif kind == "result" and self.results is not None:
            r = self.results
            title = "Solution data"
            rows.append(("File", r.path, "readonly", None, None))
            rows.append(("Ports", ", ".join(r.port_names), "readonly", None, None))
            rows.append(("Frequency points", str(len(r.freq)), "readonly", None, None))
            if r.has_farfield:
                rows.append(("Peak directivity (dBi)", f"{10 * np.log10(r.ff_Dmax):.2f}", "readonly", None, None))

        self.props_title.setText(" " + title)
        self.props.setRowCount(len(rows))
        for i, (label, value, typ, setter, options) in enumerate(rows):
            li = QTableWidgetItem(label)
            li.setFlags(Qt.ItemIsEnabled)
            self.props.setItem(i, 0, li)
            self._prop_setters.append(setter)
            if typ == "text":
                it = QTableWidgetItem(str(value))
                it.setData(ROLE, str(value))
                self.props.setItem(i, 1, it)
            elif typ == "readonly":
                it = QTableWidgetItem(str(value))
                it.setFlags(Qt.ItemIsEnabled | Qt.ItemIsSelectable)
                it.setForeground(QColor("#666"))
                self.props.setItem(i, 1, it)
            elif typ == "combo":
                cb = QComboBox()
                cb.addItems(options)
                cb.setCurrentText(str(value))
                cb.currentTextChanged.connect(lambda v, s=setter: self._commit(s, v))
                self.props.setCellWidget(i, 1, cb)
            elif typ == "check":
                chk = QCheckBox()
                chk.setChecked(bool(value))
                chk.toggled.connect(lambda v, s=setter: self._commit(s, v))
                self.props.setCellWidget(i, 1, chk)
            elif typ == "color":
                btn = QPushButton(value)
                btn.clicked.connect(lambda _=False, s=setter: s())
                self.props.setCellWidget(i, 1, btn)
        self._loading_props = False

    def _commit(self, setter, value) -> bool:
        before = copy.deepcopy(self.project.to_dict())
        try:
            setter(value)
        except Exception as exc:  # noqa: BLE001
            self.project = Project.from_dict(before)
            self.message("error", f"Invalid value '{value}': {exc}")
            self.show_properties(self.selected)
            return False
        self.undo_stack.append(before)
        self.redo_stack.clear()
        self.dirty = True
        self.refresh_tree()
        self.refresh_view()
        self._update_title()
        return True

    def _on_prop_edited(self, item: QTableWidgetItem):
        if self._loading_props or item.column() != 1:
            return
        setter = self._prop_setters[item.row()]
        if setter is None:
            return
        new = item.text().strip()
        if new == item.data(ROLE):
            return
        if self._commit(setter, new):
            self.show_properties(self.selected)

    def _renamer(self, obj):
        def s(new):
            new = new.strip()
            if not re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*", new):
                raise ValueError("names must be identifiers (letters, digits, underscore)")
            taken = {o.name for o in self.project.objects} | {p.name for p in self.project.ports}
            if new != obj.name and new in taken:
                raise ValueError(f"'{new}' already exists")
            obj.name = new
            self.selected = (self.selected[0], new)
        return s

    def _padding_setter(self):
        def s(v):
            if v.strip().lower() != "auto":
                evaluate(v, self.project.resolved_variables())
            self.project.region.padding = v.strip()
        return s

    def _units_setter(self):
        def s(v):
            self.project.units = v
        return s

    def _color_setter(self, mat: Material):
        def s():
            c = QColorDialog.getColor(QColor.fromRgbF(*mat.color), self, "Material color")
            if c.isValid():
                self._commit(lambda _: setattr(mat, "color", (c.redF(), c.greenF(), c.blueF())), None)
        return s

    # ================================================================== editing commands
    def draw(self, kind):
        self._snapshot()
        name = self.project.unique_name(kind.capitalize())
        mat = "pec" if kind == "rectangle" else "vacuum"
        obj = Primitive(name, kind, mat)
        self.project.objects.append(obj)
        self.selected = ("object", name)
        self.refresh_all()
        self.message("info", f"Created {kind} '{name}'. Edit its dimensions in the Properties panel.")

    def add_port(self):
        self._snapshot()
        name = self.project.unique_name("Port")
        self.project.ports.append(LumpedPort(name, excited=not any(p.excited for p in self.project.ports)))
        self.selected = ("port", name)
        self.refresh_all()
        self.message("info", f"Created lumped port '{name}'. Set its start/end points in Properties.")

    def delete_selected(self):
        key = self.selected
        if not key:
            return
        if key[0] == "object":
            self._snapshot()
            self.project.objects = [o for o in self.project.objects if o.name != key[1]]
        elif key[0] == "port":
            self._snapshot()
            self.project.ports = [o for o in self.project.ports if o.name != key[1]]
        elif key[0] == "material" and key[1] in self.project.materials:
            if any(o.material == key[1] for o in self.project.objects):
                self.message("error", f"Material '{key[1]}' is in use and cannot be deleted.")
                return
            self._snapshot()
            del self.project.materials[key[1]]
        else:
            return
        self.message("info", f"Deleted '{key[1]}'.")
        self.selected = None
        self.refresh_all()

    def duplicate_selected(self):
        key = self.selected
        if not key or key[0] not in ("object", "port"):
            return
        src = self._find(key[0], key[1])
        self._snapshot()
        dup = copy.deepcopy(src)
        base = re.sub(r"\d+$", "", src.name) or src.name
        dup.name = self.project.unique_name(base)
        (self.project.objects if key[0] == "object" else self.project.ports).append(dup)
        self.selected = (key[0], dup.name)
        self.refresh_all()

    def add_variable(self):
        name, ok = QInputDialog.getText(self, "Add Design Variable", "Name:")
        if not ok or not name.strip():
            return
        name = name.strip()
        if not re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*", name):
            self.message("error", f"'{name}' is not a valid variable name.")
            return
        value, ok = QInputDialog.getText(self, "Add Design Variable", f"Value of {name} ({self.project.units}):",
                                         text="1")
        if not ok:
            return
        try:
            evaluate(value, self.project.resolved_variables())
        except Exception as exc:  # noqa: BLE001
            self.message("error", str(exc))
            return
        self._snapshot()
        self.project.variables[name] = value
        self.selected = ("variables",)
        self.refresh_all()

    def delete_variable(self):
        names = list(self.project.variables)
        name, ok = QInputDialog.getItem(self, "Delete Design Variable", "Variable:", names, 0, False)
        if not ok:
            return
        self._snapshot()
        del self.project.variables[name]
        try:
            self.project.resolved_variables()
            for o in self.project.objects:
                o.values(self.project.resolved_variables())
        except Exception as exc:  # noqa: BLE001
            self.message("warn", f"'{name}' is still referenced: {exc}")
        self.refresh_all()

    def add_material(self):
        name, ok = QInputDialog.getText(self, "New Material", "Material name:")
        if not ok or not name.strip():
            return
        name = name.strip()
        if name in self.project.all_materials():
            self.message("error", f"Material '{name}' already exists.")
            return
        eps, ok = QInputDialog.getDouble(self, "New Material", "Relative permittivity:", 2.0, 1.0, 1000.0, 3)
        if not ok:
            return
        self._snapshot()
        rng = np.random.default_rng(abs(hash(name)) % 2**32)
        self.project.materials[name] = Material(name, eps, color=tuple(rng.uniform(0.3, 0.9, 3).round(3)))
        self.selected = ("material", name)
        self.refresh_all()

    def set_units(self):
        u, ok = QInputDialog.getItem(self, "Units", "Model units:", list(UNITS),
                                     list(UNITS).index(self.project.units), False)
        if ok and u != self.project.units:
            self._snapshot()
            self.project.units = u
            self.refresh_all()

    # ================================================================== files
    def _confirm_discard(self) -> bool:
        if not self.dirty:
            return True
        r = QMessageBox.question(self, APP_NAME, "Save changes to the current project?",
                                 QMessageBox.Save | QMessageBox.Discard | QMessageBox.Cancel)
        if r == QMessageBox.Save:
            return self.save_project()
        return r == QMessageBox.Discard

    def new_project(self):
        if self._confirm_discard():
            self.load_project(Project(), None)
            self.message("info", "New project.")

    def _open_example(self, factory):
        if self._confirm_discard():
            self.load_project(factory(), None)
            self.message("info", f"Opened example '{self.project.name}'.")

    def open_project(self):
        if not self._confirm_discard():
            return
        path, _ = QFileDialog.getOpenFileName(self, "Open Project", "", FILE_FILTER)
        if path:
            try:
                self.load_project(Project.load(path), path)
                self.message("info", f"Opened {path}")
            except Exception as exc:  # noqa: BLE001
                self.message("error", f"Could not open {path}: {exc}")

    def save_project(self) -> bool:
        if not self.path:
            return self.save_project_as()
        self.project.save(self.path)
        self.dirty = False
        self._update_title()
        self.message("info", f"Saved {self.path}")
        return True

    def save_project_as(self) -> bool:
        path, _ = QFileDialog.getSaveFileName(self, "Save Project As", f"{self.project.name}.ffproj", FILE_FILTER)
        if not path:
            return False
        self.path = path
        return self.save_project()

    def closeEvent(self, event):
        if self.process is not None:
            self.abort()
        if self._confirm_discard():
            self.view.close()
            event.accept()
        else:
            event.ignore()

    # ================================================================== simulation
    def validate(self) -> bool:
        p = self.project
        ok = True
        self.message("info", "Validation check:")
        try:
            variables = p.resolved_variables()
        except Exception as exc:  # noqa: BLE001
            self.message("error", f"  Design variables: {exc}")
            return False
        self.message("info", f"  Design variables: {len(variables)} OK")
        if p.source:
            try:
                from . import cad
                res = cad.replay(p)
                bad = [i for i in res.issues if i.level == "error"]
                self.message("warn" if bad else "info",
                             f"  Imported geometry: {len(res.shapes)} objects rebuilt from "
                             f"{os.path.basename(p.source['path'])}" + (f", {len(bad)} error(s)" if bad else ""))
            except Exception as exc:  # noqa: BLE001
                self.message("error", f"  Imported geometry: {exc}")
                return False
        for o in p.objects:
            try:
                o.values(variables)
            except Exception as exc:  # noqa: BLE001
                self.message("error", f"  {o.name}: {exc}")
                ok = False
        if not p.ports:
            self.message("error", "  Excitations: no ports defined")
            ok = False
        elif not any(port.excited for port in p.ports):
            self.message("error", "  Excitations: no port is excited")
            ok = False
        else:
            for port in p.ports:
                try:
                    a, b = port.corners(variables)
                    evaluate(port.impedance, variables)
                    ax = AXES.index(port.axis)
                    if a[ax] == b[ax]:
                        self.message("error", f"  {port.name}: zero length along current direction {port.axis}")
                        ok = False
                except Exception as exc:  # noqa: BLE001
                    self.message("error", f"  {port.name}: {exc}")
                    ok = False
            if ok:
                self.message("info", f"  Excitations: {len(p.ports)} port(s) OK")
        s = p.setup
        try:
            if not 0 <= float(s.f_start) < float(s.f_stop):
                raise ValueError("start frequency must be below stop frequency")
            if s.farfield and not float(s.f_start) <= float(s.farfield_freq) <= float(s.f_stop):
                raise ValueError("far-field frequency must lie inside the sweep")
            self.message("info", f"  Analysis setup: {s.f_start} - {s.f_stop} GHz OK")
        except Exception as exc:  # noqa: BLE001
            self.message("error", f"  Analysis setup: {exc}")
            ok = False
        if ok:
            try:
                from .solver import _import_openems, build_mesh, region_bounds
                _import_openems()
                lo, hi = region_bounds(p, variables)
                mesh, max_res = build_mesh(p, variables, lo, hi)
                cells = int(np.prod([len(m) for m in mesh]))
                self.message("info", f"  Mesh: {' x '.join(str(len(m)) for m in mesh)} = {cells:,} cells, "
                                     f"max cell {max_res:.3g} {p.units}, ~{cells * 200 / 1e6:.0f} MB")
            except Exception as exc:  # noqa: BLE001
                self.message("error", f"  Mesh: {exc}")
                ok = False
        self.message("info" if ok else "error", "Validation " + ("passed." if ok else "FAILED."))
        return ok

    def analyze(self):
        if self.process is not None:
            self.message("warn", "A simulation is already running.")
            return
        if not self.validate():
            return
        out = self._results_dir()
        os.makedirs(out, exist_ok=True)
        proj_file = os.path.join(out, "project_snapshot.json")
        self.project.save(proj_file)

        self.log.clear()
        self.progress.setValue(0)
        self.progress_label.setText(f"{self.project.setup.name}: starting solver...")
        self.abort_btn.setEnabled(True)
        self._t0 = time.time()
        self._last_ts, self._last_db = 0, 0.0

        proc = QProcess(self)
        env = QProcessEnvironment.systemEnvironment()
        env.insert("PYTHONUNBUFFERED", "1")
        proc.setProcessEnvironment(env)
        proc.setWorkingDirectory(HFSS_DIR)
        proc.setProcessChannelMode(QProcess.MergedChannels)
        proc.readyReadStandardOutput.connect(self._on_solver_output)
        proc.finished.connect(self._on_solver_finished)
        self.process = proc
        self._out_dir = out
        proc.start(sys.executable, ["-u", "-m", "fieldforge.solver", proj_file, out])
        self.message("info", f"Analyzing {self.project.name} : {self.project.setup.name} ...")

    def _on_solver_output(self):
        data = bytes(self.process.readAllStandardOutput()).decode("utf-8", "replace")
        for line in data.splitlines():
            if not line.strip():
                continue
            # openEMS sometimes wraps the energy figure onto its own line, so track both separately
            m_ts, m_en = TS_RE.search(line), ENERGY_RE.search(line)
            if m_ts:
                self._last_ts = int(m_ts.group(1))
            if m_en:
                self._last_db = float(m_en.group(1))
            if m_ts or m_en:
                ts, db = self._last_ts, self._last_db
                s = self.project.setup
                frac = max(ts / max(s.max_timesteps, 1), min(db / s.end_criteria_db, 1.0) if db < 0 else 0)
                self.progress.setValue(int(1000 * min(frac, 0.99)))
                self.progress_label.setText(
                    f"{s.name}: FDTD time step {ts:,}  |  energy {db:.1f} dB  /  {s.end_criteria_db:g} dB")
            elif line.startswith("[FieldForge]"):
                msg = line[len("[FieldForge]"):].strip()
                self.progress_label.setText(f"{self.project.setup.name}: {msg}")
                self.message("error" if msg.startswith("ERROR") else "info", msg)
            self.log.appendPlainText(line)

    def _on_solver_finished(self, code, _status):
        self.abort_btn.setEnabled(False)
        self.process = None
        dt = time.time() - self._t0
        if code == 0:
            self.progress.setValue(1000)
            self.progress_label.setText(f"Simulation completed in {dt:.0f} s")
            self.message("info", f"Simulation completed in {dt:.1f} s.")
            try:
                self.results = Results(os.path.join(self._out_dir, "results.npz"))
            except Exception as exc:  # noqa: BLE001
                self.message("error", f"Could not read results: {exc}")
                return
            for i in range(self.tabs.count() - 1, 0, -1):
                self.tabs.removeTab(i)
            self.refresh_tree()
            self.open_result("sparam")
            if self.results.has_farfield:
                self.open_result("ff3d", focus=False)
        else:
            self.progress.setValue(0)
            self.progress_label.setText("Simulation failed or aborted")
            self.message("error", f"Simulation stopped (exit code {code}). See the Progress log.")

    def abort(self):
        if self.process is not None:
            self.process.kill()
            self.message("warn", "Simulation aborted by user.")

    # ================================================================== results
    def open_result(self, which, focus=True):
        if self.results is None and not (which in ("sparam", "smith") and self.references):
            self.message("warn", "No solution data. Run Simulation > Analyze All first.")
            return
        if which == "compare":
            if not self.references:
                self.message("warn", "No reference loaded. Use Results > Load Reference Results first.")
                return
            if self.results is None:
                self.message("warn", "Run the simulation first (Simulation > Analyze All), then compare.")
                return
        r = self.results
        if which in ("ff3d", "ff2d") and (r is None or not r.has_farfield):
            self.message("warn", "No far-field data: enable 'Compute far field' in the analysis setup.")
            return
        titles = {"sparam": "S Parameter Plot", "smith": "Smith Chart", "zin": "Input Impedance",
                  "ff3d": "3D Far Field", "ff2d": "2D Radiation Pattern", "compare": "Comparison"}
        title = titles[which]
        for i in range(self.tabs.count()):
            if self.tabs.tabText(i) == title:
                self.tabs.setCurrentIndex(i)
                return
        if which == "sparam":
            w, ico = SParameterPlot(r, title, references=self.references), "plot"
        elif which == "smith":
            w, ico = SmithChart(r, title, references=self.references), "smith"
        elif which == "compare":
            w, ico = ComparisonPlot(r, self.references[0]), "compare"
            c = w.metrics
            self.message("info", f"Comparison with {self.references[0].label}: resonance "
                                 f"{c['a'].f_res / 1e9:.4f} vs {c['b'].f_res / 1e9:.4f} GHz "
                                 f"({c['df_res_pct']:+.2f} %), depth {c['a'].depth_db:.1f} vs "
                                 f"{c['b'].depth_db:.1f} dB"
                                 + (f", RMS difference {c['rms_db']:.2f} dB" if c['rms_db'] is not None else ""))
        elif which == "zin":
            w, ico = ImpedancePlot(r, title), "plot"
        elif which == "ff3d":
            w, ico = FarField3D(r, self.project), "farfield"
        else:
            w, ico = FarField2D(r, title), "farfield"
        idx = self.tabs.addTab(w, icon(ico), title)
        if focus:
            self.tabs.setCurrentIndex(idx)

    def _close_tab(self, idx):
        if idx == 0:
            return
        w = self.tabs.widget(idx)
        self.tabs.removeTab(idx)
        if isinstance(w, FarField3D):
            w.plotter.close()
        w.deleteLater()

    # ================================================================== Ansys import / references
    def import_aedt_dialog(self):
        if not self._confirm_discard():
            return
        path, _ = QFileDialog.getOpenFileName(self, "Import Ansys HFSS Design", "", AEDT_FILTER)
        if not path:
            return
        from .aedt.importer import import_aedt, list_designs
        QApplication.setOverrideCursor(Qt.WaitCursor)
        try:
            designs = list_designs(path)
        except Exception as exc:  # noqa: BLE001
            QApplication.restoreOverrideCursor()
            self.message("error", f"Could not read {path}: {exc}")
            return
        QApplication.restoreOverrideCursor()
        if not designs:
            self.message("error", f"{os.path.basename(path)} contains no HFSS designs.")
            return
        name = designs[0][0]
        if len(designs) > 1:
            labels = [f"{n}    [{t}]" for n, t in designs]
            choice, ok = QInputDialog.getItem(self, "Import Ansys HFSS Design",
                                              f"{os.path.basename(path)} contains {len(designs)} designs.\n"
                                              "Choose one to import:", labels, 0, False)
            if not ok:
                return
            name = designs[labels.index(choice)][0]
        self.import_aedt_file(path, name)

    def import_aedt_file(self, path, design=None, show_report=True):
        from .aedt.importer import import_aedt
        QApplication.setOverrideCursor(Qt.WaitCursor)
        t0 = time.time()
        try:
            project, report = import_aedt(path, design)
        except Exception as exc:  # noqa: BLE001
            QApplication.restoreOverrideCursor()
            self.message("error", f"Import failed: {exc}")
            return
        finally:
            QApplication.restoreOverrideCursor()
        self.load_project(project, None)
        self.import_report = report
        self.message("info", report.summary() + f" in {time.time() - t0:.1f} s")
        for i in report.errors:
            self.message("error", f"  {i}")
        for n in report.notes:
            self.message("warn", f"  {n}")
        if show_report:
            self.show_import_report()

    def show_import_report(self):
        if self.import_report is None:
            self.message("info", "No import report for this session (re-import the .aedt to see it).")
            return
        dlg = QDialog(self)
        dlg.setWindowTitle("HFSS Import Report")
        dlg.resize(760, 680)
        lay = QVBoxLayout(dlg)
        view = QTextBrowser()
        view.setHtml(self.import_report.to_html())
        lay.addWidget(view)
        bb = QDialogButtonBox(QDialogButtonBox.Ok)
        bb.accepted.connect(dlg.accept)
        lay.addWidget(bb)
        dlg.exec()

    def load_reference(self):
        start = os.path.dirname(self.project.source.get("path", "")) if self.project.source else ""
        path, _ = QFileDialog.getOpenFileName(self, "Load Reference Results", start, REF_FILTER)
        if path:
            self.add_reference(path)

    def add_reference(self, path, label=None):
        try:
            ref = refmod.load(path, label or ("HFSS" if not self.references else None))
        except Exception as exc:  # noqa: BLE001
            self.message("error", f"Could not load reference {path}: {exc}")
            return
        self._snapshot()
        self.project.references.append({"path": os.path.abspath(path), "label": ref.label})
        self.references.append(ref)
        tr = ref.reflection()
        msg = f"Loaded reference '{ref.label}' ({len(ref.traces)} trace(s)) from {os.path.basename(path)}"
        if tr is not None:
            m = refmod.metrics(tr.freq, tr.db)
            msg += f"; {tr.label} minimum {m.depth_db:.2f} dB at {m.f_res / 1e9:.4f} GHz"
        self.message("info", msg)
        self._close_result_tabs()
        self.refresh_tree()
        self._update_title()
        self.open_result("compare" if self.results is not None else "sparam")

    def remove_reference(self, k):
        self._snapshot()
        del self.project.references[k]
        del self.references[k]
        self._close_result_tabs()
        self.refresh_all()

    def _ref_label_setter(self, k):
        def s(v):
            self.project.references[k]["label"] = v
            self.references[k].label = v
            self._close_result_tabs()
        return s

    def _close_result_tabs(self):
        for i in range(self.tabs.count() - 1, 0, -1):
            if self.tabs.tabText(i) in ("S Parameter Plot", "Smith Chart", "Comparison"):
                self._close_tab(i)

    def export_touchstone(self):
        if self.results is None:
            self.message("warn", "No solution data to export.")
            return
        r = self.results
        n = len(r.port_names)
        path, _ = QFileDialog.getSaveFileName(self, "Export Touchstone", f"{self.project.name}.s{n}p",
                                              f"Touchstone (*.s{n}p)")
        if not path:
            return
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(f"! {APP_NAME} export - column {r.port_names[r.excited]} of the S-matrix only\n")
            fh.write("# Hz S RI R 50\n")
            for k, f in enumerate(r.freq):
                vals = " ".join(f"{r.S[i][k].real:.8e} {r.S[i][k].imag:.8e}" for i in range(n))
                fh.write(f"{f:.6e} {vals}\n")
        self.message("info", f"Exported {path}")

    def about(self):
        QMessageBox.about(self, f"About {APP_NAME}",
                          f"<b>{APP_NAME}</b><br>A 3D electromagnetic simulation environment.<br><br>"
                          "Solver: openEMS (EC-FDTD), CSXCAD geometry kernel.<br>"
                          "GUI: PySide6 + PyVista/VTK, matplotlib.")


def apply_light_palette(app):
    """Force a light engineering-tool look regardless of the OS dark-mode setting."""
    from PySide6.QtGui import QPalette
    app.styleHints().setColorScheme(Qt.ColorScheme.Light)
    pal = QPalette()
    for role, color in [
        (QPalette.Window, "#f0f0f0"), (QPalette.WindowText, "#000000"), (QPalette.Base, "#ffffff"),
        (QPalette.AlternateBase, "#f5f8fc"), (QPalette.Text, "#000000"), (QPalette.Button, "#f0f0f0"),
        (QPalette.ButtonText, "#000000"), (QPalette.Highlight, "#3d7bd9"),
        (QPalette.HighlightedText, "#ffffff"), (QPalette.ToolTipBase, "#ffffe0"),
        (QPalette.ToolTipText, "#000000"), (QPalette.PlaceholderText, "#808080"),
    ]:
        pal.setColor(role, QColor(color))
    pal.setColor(QPalette.Disabled, QPalette.Text, QColor("#a0a0a0"))
    pal.setColor(QPalette.Disabled, QPalette.ButtonText, QColor("#a0a0a0"))
    app.setPalette(pal)
    app.setStyleSheet("QDockWidget::title { background: #dde5f0; padding: 4px; }"
                      "QTabBar::tab:selected { background: #ffffff; }")


def main():
    from .solver import _import_openems
    try:
        _import_openems()
    except Exception as exc:  # noqa: BLE001
        print(f"warning: openEMS not importable ({exc}); simulation will be unavailable", file=sys.stderr)
    app = QApplication.instance() or QApplication(sys.argv)
    app.setApplicationName(APP_NAME)
    app.setStyle("Fusion")
    apply_light_palette(app)
    win = MainWindow()
    win.show()
    if len(sys.argv) > 1 and os.path.exists(sys.argv[1]):
        win.load_project(Project.load(sys.argv[1]), os.path.abspath(sys.argv[1]))
    return app.exec()


if __name__ == "__main__":
    sys.exit(main())
