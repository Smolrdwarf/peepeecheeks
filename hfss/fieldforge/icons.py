"""Small vector icons drawn with QPainter so the app needs no image assets."""
from __future__ import annotations

from PySide6.QtCore import QPointF, QRectF, Qt
from PySide6.QtGui import QBrush, QColor, QIcon, QPainter, QPen, QPixmap, QPolygonF

SIZE = 32


def _canvas():
    pm = QPixmap(SIZE, SIZE)
    pm.fill(Qt.transparent)
    p = QPainter(pm)
    p.setRenderHint(QPainter.Antialiasing)
    p.setPen(QPen(QColor("#303030"), 1.6))
    return pm, p


def _poly(p, pts, color):
    p.setBrush(QBrush(QColor(color)))
    p.drawPolygon(QPolygonF([QPointF(x, y) for x, y in pts]))


def icon(name: str) -> QIcon:
    pm, p = _canvas()
    if name == "box":
        _poly(p, [(6, 11), (18, 11), (18, 27), (6, 27)], "#7fa7d9")
        _poly(p, [(6, 11), (13, 5), (25, 5), (18, 11)], "#a9c6ea")
        _poly(p, [(18, 11), (25, 5), (25, 21), (18, 27)], "#5b86bd")
    elif name == "cylinder":
        p.setBrush(QBrush(QColor("#7fa7d9")))
        p.drawRect(QRectF(7, 8, 18, 17))
        p.drawEllipse(QRectF(7, 21, 18, 7))
        p.setPen(Qt.NoPen)
        p.drawRect(QRectF(8, 9, 16, 15))
        p.setPen(QPen(QColor("#303030"), 1.6))
        p.setBrush(QBrush(QColor("#a9c6ea")))
        p.drawEllipse(QRectF(7, 4, 18, 8))
    elif name == "sphere":
        from PySide6.QtGui import QRadialGradient
        g = QRadialGradient(12, 12, 16)
        g.setColorAt(0, QColor("#d8e6f7"))
        g.setColorAt(1, QColor("#4f7db5"))
        p.setBrush(QBrush(g))
        p.drawEllipse(QRectF(4, 4, 24, 24))
    elif name == "rectangle":
        _poly(p, [(3, 22), (13, 10), (29, 10), (19, 22)], "#e3c26a")
    elif name == "port":
        _poly(p, [(10, 4), (22, 4), (22, 28), (10, 28)], "#e05050")
        p.setPen(QPen(QColor("#ffffff"), 2.2))
        p.drawLine(16, 24, 16, 9)
        p.drawLine(16, 8, 12, 13)
        p.drawLine(16, 8, 20, 13)
    elif name == "analyze":
        p.setBrush(QBrush(QColor("#2fa84f")))
        p.drawEllipse(QRectF(3, 3, 26, 26))
        p.setPen(Qt.NoPen)
        p.setBrush(QBrush(QColor("white")))
        p.drawPolygon(QPolygonF([QPointF(13, 9), QPointF(23, 16), QPointF(13, 23)]))
    elif name == "validate":
        p.setBrush(QBrush(QColor("#ffffff")))
        p.drawRect(QRectF(5, 4, 20, 24))
        p.setPen(QPen(QColor("#2fa84f"), 3))
        p.drawLine(9, 17, 14, 22)
        p.drawLine(14, 22, 23, 10)
    elif name == "abort":
        p.setBrush(QBrush(QColor("#d03030")))
        p.drawRect(QRectF(7, 7, 18, 18))
    elif name == "setup":
        p.setBrush(QBrush(QColor("#c0c0c0")))
        p.drawEllipse(QRectF(5, 5, 22, 22))
        p.setBrush(QBrush(QColor("#ffffff")))
        p.drawEllipse(QRectF(12, 12, 8, 8))
    elif name == "region":
        p.setBrush(Qt.NoBrush)
        p.setPen(QPen(QColor("#5070a0"), 1.6, Qt.DashLine))
        p.drawRect(QRectF(4, 8, 20, 20))
        p.drawLine(4, 8, 10, 3)
        p.drawLine(24, 8, 29, 3)
        p.drawLine(10, 3, 29, 3)
        p.drawLine(29, 3, 29, 22)
        p.drawLine(24, 28, 29, 22)
    elif name == "plot":
        p.setBrush(QBrush(QColor("#ffffff")))
        p.drawRect(QRectF(3, 4, 26, 24))
        p.setPen(QPen(QColor("#d62728"), 2))
        p.drawPolyline(QPolygonF([QPointF(5, 9), QPointF(11, 11), QPointF(16, 25),
                                  QPointF(20, 12), QPointF(27, 9)]))
    elif name == "farfield":
        p.setBrush(QBrush(QColor("#f0a030")))
        p.drawEllipse(QRectF(9, 2, 14, 15))
        p.setBrush(QBrush(QColor("#30a0f0")))
        p.drawEllipse(QRectF(9, 15, 14, 15))
    elif name == "smith":
        p.setBrush(QBrush(QColor("#ffffff")))
        p.drawEllipse(QRectF(3, 3, 26, 26))
        p.drawEllipse(QRectF(13, 9, 15, 14))
        p.drawLine(3, 16, 29, 16)
    elif name == "variable":
        p.setFont(p.font())
        f = p.font()
        f.setBold(True)
        f.setPointSize(13)
        p.setFont(f)
        p.setPen(QColor("#6a3fa0"))
        p.drawText(QRectF(0, 0, SIZE, SIZE), Qt.AlignCenter, "x=")
    elif name == "fit":
        p.setBrush(Qt.NoBrush)
        for x1, y1, x2, y2 in [(4, 11, 4, 4), (4, 4, 11, 4), (21, 4, 28, 4), (28, 4, 28, 11),
                               (4, 21, 4, 28), (4, 28, 11, 28), (21, 28, 28, 28), (28, 21, 28, 28)]:
            p.drawLine(x1, y1, x2, y2)
        p.drawRect(QRectF(11, 11, 10, 10))
    else:
        p.setBrush(QBrush(QColor("#c0c0c0")))
        p.drawRect(QRectF(6, 6, 20, 20))
    p.end()
    return QIcon(pm)


def color_icon(rgb) -> QIcon:
    pm, p = _canvas()
    p.setBrush(QBrush(QColor.fromRgbF(*rgb)))
    p.drawRect(QRectF(6, 6, 20, 20))
    p.end()
    return QIcon(pm)
