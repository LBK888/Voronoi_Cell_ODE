import json
import math
from PyQt6.QtWidgets import (
    QGraphicsScene, QGraphicsView, QGraphicsItem, QGraphicsPathItem,
    QMenu, QGraphicsSceneMouseEvent, QGraphicsProxyWidget
)
from PyQt6.QtCore import Qt, QPointF, QRectF, pyqtSignal, QObject
from PyQt6.QtGui import QPen, QBrush, QColor, QPainterPath, QKeySequence, QAction

class NodeSocket(QGraphicsItem):
    def __init__(self, node, socket_type="input", name="", value_type="Any"):
        super().__init__(node)
        self.node = node
        self.socket_type = socket_type  # "input" or "output"
        self.name = name
        self.value_type = value_type # "Any", "Float", "Array"
        self.radius = 6
        self.edges = []
        
        if self.socket_type in ["seq_in", "seq_out"]:
            self.color = QColor(50, 150, 255) # Blue
        elif self.value_type == "Array":
            self.color = QColor(255, 200, 50)  # Yellow
        elif self.value_type == "Float":
            self.color = QColor(100, 200, 255) # Light Blue
        else:
            self.color = QColor(200, 200, 200) # Gray
            
        self.setAcceptHoverEvents(True)
        self.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsSelectable, False)
        
    def can_connect(self, other_socket):
        if self.node == other_socket.node:
            return False
        if self.socket_type == other_socket.socket_type:
            return False
            
        is_seq_self = self.socket_type in ["seq_in", "seq_out"]
        is_seq_other = other_socket.socket_type in ["seq_in", "seq_out"]
        
        if is_seq_self or is_seq_other:
            return is_seq_self and is_seq_other
            
        if self.value_type == "Any" or other_socket.value_type == "Any":
            return True
        if self.value_type == other_socket.value_type:
            return True
        return False
        
    def boundingRect(self):
        return QRectF(-self.radius-2, -self.radius-2, 2*self.radius+4, 2*self.radius+4)

    def paint(self, painter, option, widget=None):
        painter.setBrush(QBrush(self.color))
        painter.setPen(QPen(Qt.GlobalColor.black, 1))
        painter.drawEllipse(QPointF(0, 0), self.radius, self.radius)

    def add_edge(self, edge):
        # seq_out only allows 1 edge max
        if self.socket_type == "seq_out":
            for e in list(self.edges):
                e.remove()
        if edge not in self.edges:
            self.edges.append(edge)

    def remove_edge(self, edge):
        if edge in self.edges:
            self.edges.remove(edge)

    def scene_pos(self):
        return self.scenePos()

    def mousePressEvent(self, event):
        if event.button() == Qt.MouseButton.LeftButton:
            if self.socket_type in ["input", "seq_in"]:
                if len(self.edges) > 0:
                    for e in list(self.edges):
                        e.remove()
                    return
            self.scene().start_drag_edge(self)
        else:
            super().mousePressEvent(event)


class NodeEdge(QGraphicsPathItem):
    def __init__(self, source_socket, dest_socket=None, edge_type="Data"):
        super().__init__()
        self.source_socket = source_socket
        self.dest_socket = dest_socket
        self.edge_type = edge_type  # Default to "Data" now
        self.draft_pos = None

        self.setZValue(-1)
        self.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsSelectable, True)
        self.update_path()
        self.update_style()

    def update_style(self):
        if self.isSelected():
            color = QColor(255, 255, 255)
            width = 4
        else:
            if self.edge_type == "Sequence":
                color = QColor(50, 150, 255)
                width = 4
            else:
                color = QColor(255, 255, 255)
                width = 3
        
        pen = QPen(color, width)
        self.setPen(pen)

    def shape(self):
        from PyQt6.QtGui import QPainterPathStroker
        stroker = QPainterPathStroker()
        stroker.setWidth(10)
        return stroker.createStroke(self.path())

    def update_path(self):
        if not self.source_socket:
            return

        source_pos = self.source_socket.scene_pos()
        dest_pos = self.draft_pos if self.dest_socket is None else self.dest_socket.scene_pos()

        if dest_pos is None:
            return

        routing_style = "Spline"
        if self.scene() and hasattr(self.scene(), "routing_style"):
            routing_style = self.scene().routing_style

        path = QPainterPath(source_pos)
        
        if self.edge_type == "Sequence":
            if routing_style == "Spline":
                dist_y = abs(dest_pos.y() - source_pos.y()) * 0.5
                dist_y = max(dist_y, 40)
                if dest_pos.y() < source_pos.y():
                    dist_y += (source_pos.y() - dest_pos.y()) * 0.5
                ctrl1 = QPointF(source_pos.x(), source_pos.y() + dist_y if self.source_socket.socket_type == "seq_out" else source_pos.y() - dist_y)
                ctrl2 = QPointF(dest_pos.x(), dest_pos.y() - dist_y if self.source_socket.socket_type == "seq_out" else dest_pos.y() + dist_y)
                path.cubicTo(ctrl1, ctrl2, dest_pos)
            elif routing_style == "Orthogonal":
                mid_y = (source_pos.y() + dest_pos.y()) / 2
                if dest_pos.y() > source_pos.y() + 20:
                    path.lineTo(source_pos.x(), mid_y)
                    path.lineTo(dest_pos.x(), mid_y)
                    path.lineTo(dest_pos.x(), dest_pos.y())
                else:
                    path.lineTo(source_pos.x(), source_pos.y() + 30)
                    mid_x = (source_pos.x() + dest_pos.x()) / 2
                    path.lineTo(mid_x, source_pos.y() + 30)
                    path.lineTo(mid_x, dest_pos.y() - 30)
                    path.lineTo(dest_pos.x(), dest_pos.y() - 30)
                    path.lineTo(dest_pos.x(), dest_pos.y())
            else:
                path.lineTo(dest_pos)
        else:
            if routing_style == "Spline":
                dist_x = abs(dest_pos.x() - source_pos.x()) * 0.5
                dist_x = max(dist_x, 50)
                if dest_pos.x() < source_pos.x():
                    dist_x += (source_pos.x() - dest_pos.x()) * 0.5
                ctrl1 = QPointF(source_pos.x() + dist_x, source_pos.y())
                ctrl2 = QPointF(dest_pos.x() - dist_x, dest_pos.y())
                path.cubicTo(ctrl1, ctrl2, dest_pos)
            elif routing_style == "Orthogonal":
                mid_x = (source_pos.x() + dest_pos.x()) / 2
                if dest_pos.x() > source_pos.x() + 20: 
                    path.lineTo(mid_x, source_pos.y())
                    path.lineTo(mid_x, dest_pos.y())
                    path.lineTo(dest_pos.x(), dest_pos.y())
                else:
                    path.lineTo(source_pos.x() + 30, source_pos.y())
                    mid_y = (source_pos.y() + dest_pos.y()) / 2
                    path.lineTo(source_pos.x() + 30, mid_y)
                    path.lineTo(dest_pos.x() - 30, mid_y)
                    path.lineTo(dest_pos.x() - 30, dest_pos.y())
                    path.lineTo(dest_pos.x(), dest_pos.y())
            else:
                path.lineTo(dest_pos)

        self.setPath(path)

    def paint(self, painter, option, widget=None):
        self.update_style()
        super().paint(painter, option, widget)

    def remove(self):
        if self.source_socket:
            self.source_socket.remove_edge(self)
        if self.dest_socket:
            self.dest_socket.remove_edge(self)
        if self.scene():
            self.scene().removeItem(self)

    def toggle_type(self):
        self.edge_type = "Repression" if self.edge_type == "Activation" else "Activation"
        self.update_style()
        
    def mouseDoubleClickEvent(self, event):
        if event.button() == Qt.MouseButton.LeftButton:
            self.toggle_type()
            event.accept()
        else:
            super().mouseDoubleClickEvent(event)


class NodeScene(QGraphicsScene):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.drag_edge = None
        self.drag_start_socket = None

    def start_drag_edge(self, socket):
        self.drag_start_socket = socket
        edge_type = "Data"
        if hasattr(socket, 'socket_type') and socket.socket_type in ["seq_in", "seq_out"]:
            edge_type = "Sequence"
        elif hasattr(socket, 'name'):
            if socket.name == "Inhibit (Out)":
                edge_type = "Repression"
            elif socket.name == "Active (Out)":
                edge_type = "Activation"
                
        self.drag_edge = NodeEdge(socket, None, edge_type)
        self.drag_edge.draft_pos = socket.scene_pos()
        self.addItem(self.drag_edge)

    def mouseMoveEvent(self, event):
        if self.drag_edge:
            self.drag_edge.draft_pos = event.scenePos()
            self.drag_edge.update_path()
        else:
            super().mouseMoveEvent(event)

    def mouseReleaseEvent(self, event):
        if self.drag_edge:
            # Find socket under cursor
            item = self.itemAt(event.scenePos(), self.views()[0].transform())
            if isinstance(item, NodeSocket) and item != self.drag_start_socket:
                if hasattr(self.drag_start_socket, 'can_connect') and self.drag_start_socket.can_connect(item):
                    # Successfully connected
                    is_seq = self.drag_start_socket.socket_type in ["seq_in", "seq_out"]
                    
                    if (self.drag_start_socket.socket_type == "output" or self.drag_start_socket.socket_type == "seq_out"):
                        self.drag_edge.source_socket = self.drag_start_socket
                        self.drag_edge.dest_socket = item
                    else:
                        self.drag_edge.source_socket = item
                        self.drag_edge.dest_socket = self.drag_start_socket
                    
                    self.drag_start_socket.add_edge(self.drag_edge)
                    item.add_edge(self.drag_edge)
                    self.drag_edge.update_path()
                    self.drag_edge = None
                    self.drag_start_socket = None
                    return
            
            # Connection failed or dropped
            self.removeItem(self.drag_edge)
            self.drag_edge = None
            self.drag_start_socket = None
        else:
            super().mouseReleaseEvent(event)

    def arrange_nodes(self):
        nodes = [item for item in self.items() if hasattr(item, 'inputs') and hasattr(item, 'outputs')]
        if not nodes: return
        
        layers = {}
        visiting = set()
        
        def get_layer(node):
            if node in layers:
                return layers[node]
            if node in visiting:
                return 0 # cyclic fallback
            visiting.add(node)
                
            max_p = -1
            for sock in getattr(node, 'inputs', []):
                for edge in sock.edges:
                    p_sock = edge.source_socket if hasattr(edge, 'source_socket') and edge.source_socket and edge.source_socket.socket_type == 'output' else edge.dest_socket
                    if p_sock and p_sock.node != node:
                        max_p = max(max_p, get_layer(p_sock.node))
                        
            visiting.remove(node)
            layers[node] = max_p + 1
            return layers[node]

        for n in nodes:
            get_layer(n)
            
        layer_groups = {}
        for n, l in layers.items():
            layer_groups.setdefault(l, []).append(n)
            
        col_spacing = 350
        padding = 50
        
        for l in sorted(layer_groups.keys()):
            x_pos = l * col_spacing
            y_pos = 100
            layer_nodes = sorted(layer_groups[l], key=lambda n: n.pos().y())
            for n in layer_nodes:
                n.setPos(x_pos, y_pos)
                if hasattr(n, 'boundingRect'):
                    y_pos += n.boundingRect().height() + padding
                
        for n in nodes:
            for sock in getattr(n, 'inputs', []) + getattr(n, 'outputs', []):
                for edge in sock.edges:
                    if hasattr(edge, 'update_path'):
                        edge.update_path()



class NodeView(QGraphicsView):
    def __init__(self, scene, parent=None):
        super().__init__(scene, parent)
        self.setRenderHint(self.renderHints() | self.renderHints().Antialiasing)
        self.setBackgroundBrush(QBrush(QColor(30, 30, 30)))
        
        self.setDragMode(QGraphicsView.DragMode.NoDrag)
        self.setViewportUpdateMode(QGraphicsView.ViewportUpdateMode.FullViewportUpdate)
        
        self._is_panning = False
        self._pan_start = QPointF()

    def mousePressEvent(self, event):
        if event.button() == Qt.MouseButton.LeftButton:
            if event.modifiers() & Qt.KeyboardModifier.ControlModifier:
                self.setDragMode(QGraphicsView.DragMode.RubberBandDrag)
                super().mousePressEvent(event)
            else:
                item = self.itemAt(event.pos())
                if not item:
                    self._is_panning = True
                    self._pan_start = event.position()
                    self.setCursor(Qt.CursorShape.ClosedHandCursor)
                else:
                    self.setDragMode(QGraphicsView.DragMode.NoDrag)
                    super().mousePressEvent(event)
        elif event.button() == Qt.MouseButton.MiddleButton:
            self._is_panning = True
            self._pan_start = event.position()
            self.setCursor(Qt.CursorShape.ClosedHandCursor)
        else:
            super().mousePressEvent(event)

    def mouseMoveEvent(self, event):
        if self._is_panning:
            delta = event.position() - self._pan_start
            self._pan_start = event.position()
            
            h_bar = self.horizontalScrollBar()
            v_bar = self.verticalScrollBar()
            
            old_h, old_v = h_bar.value(), v_bar.value()
            new_h, new_v = old_h - int(delta.x()), old_v - int(delta.y())
            
            h_bar.setValue(new_h)
            v_bar.setValue(new_v)
            
            if h_bar.value() != new_h or v_bar.value() != new_v:
                scene = self.scene()
                if scene:
                    vr = self.mapToScene(self.viewport().rect()).boundingRect()
                    pad_x, pad_y = vr.width() * 2, vr.height() * 2
                    scene.setSceneRect(scene.sceneRect().united(vr.adjusted(-pad_x, -pad_y, pad_x, pad_y)))
                    h_bar.setValue(new_h)
                    v_bar.setValue(new_v)
        else:
            super().mouseMoveEvent(event)

    def mouseReleaseEvent(self, event):
        if self._is_panning and (event.button() == Qt.MouseButton.LeftButton or event.button() == Qt.MouseButton.MiddleButton):
            self._is_panning = False
            self.setCursor(Qt.CursorShape.ArrowCursor)
        else:
            super().mouseReleaseEvent(event)
        
        if not self._is_panning:
            self.setDragMode(QGraphicsView.DragMode.NoDrag)

    def wheelEvent(self, event):
        # Zoom support
        zoom_in_factor = 1.25
        zoom_out_factor = 1 / zoom_in_factor

        self.setTransformationAnchor(QGraphicsView.ViewportAnchor.AnchorUnderMouse)
        if event.angleDelta().y() > 0:
            zoom_factor = zoom_in_factor
        else:
            zoom_factor = zoom_out_factor

        self.scale(zoom_factor, zoom_factor)
        
        scene = self.scene()
        if scene:
            vr = self.mapToScene(self.viewport().rect()).boundingRect()
            pad_x, pad_y = vr.width() * 2, vr.height() * 2
            scene.setSceneRect(scene.sceneRect().united(vr.adjusted(-pad_x, -pad_y, pad_x, pad_y)))
