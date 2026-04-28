import json
from PyQt6.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QPushButton, QMenu, QFileDialog, QMessageBox, QInputDialog,
    QWidget, QListWidget, QListWidgetItem, QSplitter, QLabel
)
from PyQt6.QtCore import Qt
from PyQt6.QtGui import QAction

from .node_core import NodeScene, NodeView, NodeEdge
from .node_types import GeneNode, CellNode, MathNode, TextNode, VariableReadNode, ParameterNode, NamedConstantNode, StateNode
from .node_parser import PythonToNodeParser, NodeToPythonGenerator
from PyQt6.QtGui import QPainter, QColor, QPen, QBrush
from PyQt6.QtCore import QRectF, QPointF, QTimer

class MinimapWidget(QWidget):
    def __init__(self, view, parent=None):
        super().__init__(parent)
        self.view = view
        self.setFixedSize(200, 150)
        self.setStyleSheet("background-color: rgba(20, 20, 20, 200); border: 1px solid #555;")
        self.setCursor(Qt.CursorShape.PointingHandCursor)
        self.timer = QTimer(self)
        self.timer.timeout.connect(self.update)
        self.timer.start(50)

    def get_items_rect(self):
        scene = self.view.scene()
        from .node_types import BaseNode
        min_x, min_y, max_x, max_y = float('inf'), float('inf'), float('-inf'), float('-inf')
        has_nodes = False
        for item in scene.items():
            if isinstance(item, BaseNode):
                has_nodes = True
                rect = item.sceneBoundingRect()
                min_x = min(min_x, rect.left())
                min_y = min(min_y, rect.top())
                max_x = max(max_x, rect.right())
                max_y = max(max_y, rect.bottom())

        if not has_nodes:
            items_rect = QRectF(-1000, -1000, 2000, 2000)
        else:
            items_rect = QRectF(QPointF(min_x, min_y), QPointF(max_x, max_y))
            
        items_rect.adjust(-items_rect.width()*0.1 - 100, -items_rect.height()*0.1 - 100, 
                          items_rect.width()*0.1 + 100, items_rect.height()*0.1 + 100)
        if items_rect.width() == 0: items_rect.setWidth(1)
        if items_rect.height() == 0: items_rect.setHeight(1)
        return items_rect

    def paintEvent(self, event):
        super().paintEvent(event)
        scene = self.view.scene()
        if not scene: return
        
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        
        items_rect = self.get_items_rect()
        
        scale_x = self.width() / items_rect.width()
        scale_y = self.height() / items_rect.height()
        scale = min(scale_x, scale_y)
        
        offset_x = (self.width() - items_rect.width() * scale) / 2
        offset_y = (self.height() - items_rect.height() * scale) / 2
        
        def map_to_minimap(scene_pos):
            x = (scene_pos.x() - items_rect.left()) * scale + offset_x
            y = (scene_pos.y() - items_rect.top()) * scale + offset_y
            return QPointF(x, y)
            
        from .node_types import BaseNode
        from .node_core import NodeEdge
        
        painter.setPen(Qt.PenStyle.NoPen)
        for item in scene.items():
            if isinstance(item, BaseNode):
                pos = map_to_minimap(item.scenePos())
                w = item.width * scale
                h = item.height * scale
                painter.setBrush(QBrush(item.title_color if item.title_color.alpha() > 0 else item.bg_color))
                painter.drawRect(QRectF(pos.x(), pos.y(), w, h))
                
        painter.setPen(QPen(QColor(150, 150, 150, 100), 1))
        for item in scene.items():
            if isinstance(item, NodeEdge) and item.source_socket and item.dest_socket:
                p1 = map_to_minimap(item.source_socket.scenePos())
                p2 = map_to_minimap(item.dest_socket.scenePos())
                painter.drawLine(p1, p2)
                
        viewport_rect = self.view.mapToScene(self.view.viewport().rect()).boundingRect()
        vp_pos = map_to_minimap(viewport_rect.topLeft())
        vp_w = viewport_rect.width() * scale
        vp_h = viewport_rect.height() * scale
        
        painter.setPen(QPen(QColor(255, 255, 255, 200), 1))
        painter.setBrush(QBrush(QColor(255, 255, 255, 30)))
        painter.drawRect(QRectF(vp_pos.x(), vp_pos.y(), vp_w, vp_h))

    def move_view(self, event):
        scene = self.view.scene()
        if not scene: return
        
        items_rect = self.get_items_rect()
        
        scale_x = self.width() / items_rect.width()
        scale_y = self.height() / items_rect.height()
        scale = min(scale_x, scale_y)
        
        offset_x = (self.width() - items_rect.width() * scale) / 2
        offset_y = (self.height() - items_rect.height() * scale) / 2
        
        scene_x = (event.pos().x() - offset_x) / scale + items_rect.left()
        scene_y = (event.pos().y() - offset_y) / scale + items_rect.top()
        
        self.view.centerOn(scene_x, scene_y)

    def mousePressEvent(self, event):
        if event.button() == Qt.MouseButton.LeftButton:
            self.move_view(event)
            
    def mouseMoveEvent(self, event):
        if event.buttons() & Qt.MouseButton.LeftButton:
            self.move_view(event)


class VariableSidebar(QWidget):
    def __init__(self, scene, parent=None):
        super().__init__(parent)
        self.scene = scene
        self.layout = QVBoxLayout(self)
        self.layout.setContentsMargins(0, 0, 0, 0)
        
        self.label = QLabel("Variables & Constants")
        self.label.setStyleSheet("font-weight: bold; color: #ddd; padding: 5px;")
        self.layout.addWidget(self.label)
        
        self.list_widget = QListWidget()
        self.list_widget.setDragEnabled(True)
        self.layout.addWidget(self.list_widget)
        
        self.timer = QTimer(self)
        self.timer.timeout.connect(self.refresh_list)
        self.timer.start(1000)
        
        self.current_items = set()

    def refresh_list(self):
        new_items = set()
        
        from .node_types import StateNode, ParameterNode, NamedConstantNode
        
        for item in self.scene.items():
            if type(item).__name__ == "StateNode":
                var_name = getattr(item, "var_name", "")
                if var_name:
                    new_items.add(("Variable", var_name))
                else:
                    new_items.add(("Variable", f"Y[:, {item.index_edit.value()}]"))
            elif type(item).__name__ == "ParameterNode":
                p_name = item.param_name_edit.text()
                if p_name:
                    new_items.add(("Parameter", p_name))
            elif type(item).__name__ == "NamedConstantNode":
                c_name = item.name_edit.text()
                if c_name:
                    new_items.add(("Constant", c_name))
                    
        if new_items != self.current_items:
            self.current_items = new_items
            self.list_widget.clear()
            for t, name in sorted(list(self.current_items)):
                item = QListWidgetItem(f"[{t}] {name}")
                item.setData(Qt.ItemDataRole.UserRole, (t, name))
                self.list_widget.addItem(item)

class NodeEditorDialog(QDialog):

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Node-based ODE Editor")
        self.resize(1300, 800)
        
        main_layout = QVBoxLayout(self)
        
        # Toolbar
        toolbar = QHBoxLayout()
        
        self.add_node_btn = QPushButton("Add Node...")
        self.apply_ode_btn = QPushButton("Apply ODE Text to Nodes")
        self.load_config_btn = QPushButton("Load Config Template")
        self.export_btn = QPushButton("Export JSON")
        self.import_btn = QPushButton("Import JSON")
        self.generate_code_btn = QPushButton("Generate ODE Code & Apply")
        
        self.routing_btn = QPushButton("Style: Spline")
        self.scene = NodeScene(self)
        self.scene.routing_style = "Spline"
        
        toolbar.addWidget(self.add_node_btn)
        toolbar.addWidget(self.apply_ode_btn)
        toolbar.addWidget(self.load_config_btn)
        toolbar.addWidget(self.export_btn)
        toolbar.addWidget(self.import_btn)
        toolbar.addWidget(self.generate_code_btn)
        toolbar.addWidget(self.routing_btn)
        toolbar.addStretch()
        main_layout.addLayout(toolbar)
        
        # Splitter
        splitter = QSplitter(Qt.Orientation.Horizontal)
        main_layout.addWidget(splitter)
        
        # Scene & View
        self.scene.setSceneRect(-2000, -2000, 4000, 4000)
        self.view = NodeView(self.scene, self)
        
        self.view.setAcceptDrops(True)
        self.view.dragEnterEvent = self.view_dragEnterEvent
        self.view.dragMoveEvent = self.view_dragMoveEvent
        self.view.dropEvent = self.view_dropEvent
        
        splitter.addWidget(self.view)
        
        self.sidebar = VariableSidebar(self.scene, self)
        splitter.addWidget(self.sidebar)
        splitter.setSizes([1000, 250])
        
        # Minimap
        self.minimap = MinimapWidget(self.view, self.view)
        
        # Connections
        self.add_node_btn.clicked.connect(self.show_add_menu)
        self.apply_ode_btn.clicked.connect(self.apply_ode_text)
        self.load_config_btn.clicked.connect(self.load_config_template)
        self.export_btn.clicked.connect(self.export_json)
        self.import_btn.clicked.connect(self.import_json)
        self.generate_code_btn.clicked.connect(self.generate_code)
        self.routing_btn.clicked.connect(self.toggle_routing)

        # Setup Parser/Generator
        self.parser = PythonToNodeParser(self.scene)
        self.generator = NodeToPythonGenerator(self.scene)
        
        # Override view double click and resize
        self.view.mouseDoubleClickEvent = self.view_double_click
        self.view.resizeEvent = self.view_resizeEvent
        
        self.keyPressEvent = self.handle_shortcuts
        

    def view_dragEnterEvent(self, event):
        if event.mimeData().hasFormat("application/x-qabstractitemmodeldatalist"):
            event.acceptProposedAction()

    def view_dragMoveEvent(self, event):
        if event.mimeData().hasFormat("application/x-qabstractitemmodeldatalist"):
            event.acceptProposedAction()

    def view_dropEvent(self, event):
        if event.mimeData().hasFormat("application/x-qabstractitemmodeldatalist"):
            item = self.sidebar.list_widget.currentItem()
            if not item: return
            item_data = item.data(Qt.ItemDataRole.UserRole)
            if item_data:
                node_type, name = item_data
                if hasattr(event, 'position'):
                    pos = event.position().toPoint()
                else:
                    pos = event.pos()
                scene_pos = self.view.mapToScene(pos)
                
                from .node_types import VariableReadNode, ParameterNode
                
                if node_type == "Variable" or node_type == "Constant":
                    node = VariableReadNode()
                    node.var_combo.setCurrentText(name)
                elif node_type == "Parameter":
                    node = ParameterNode()
                    node.param_name_edit.setText(name)
                    node.access_style = "get"
                else:
                    return
                    
                node.setPos(scene_pos.x() - node.width/2, scene_pos.y() - node.height/2)
                self.scene.addItem(node)
                event.acceptProposedAction()

    def load_config_template(self):
        filename, _ = QFileDialog.getOpenFileName(self, "Load Config JSON", "", "JSON Files (*.json)")
        if filename:
            try:
                import json
                with open(filename, 'r', encoding='utf-8') as f:
                    data = json.load(f)
                    
                ode_str = data.get('ode', '')
                if not ode_str:
                    ode_str = data.get('ode_code', '')
                params_dict = data.get('params', {})
                # Format params exactly as the editor does (as a stringified dict or similar)
                import pprint
                params_str = pprint.pformat(params_dict) if isinstance(params_dict, dict) else str(params_dict)
                color_func_str = data.get('color_func', '')
                
                self.parser.parse_ode_string(ode_str, params_str, color_func_str)
                QMessageBox.information(self, "Success", "Loaded config and parsed to nodes successfully.")
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Failed to load config: {str(e)}")

    def toggle_routing(self):
        if self.scene.routing_style == "Spline":
            self.scene.routing_style = "Orthogonal"
            self.routing_btn.setText("Style: Orthogonal")
        else:
            self.scene.routing_style = "Spline"
            self.routing_btn.setText("Style: Spline")
            
        from .node_core import NodeEdge
        for item in self.scene.items():
            if isinstance(item, NodeEdge):
                item.update_path()

    def view_resizeEvent(self, event):
        super(NodeView, self.view).resizeEvent(event)
        # Position minimap at bottom right of the view
        if hasattr(self, 'minimap') and hasattr(self, 'view'):
            margin = 20
            # map from view's local coordinates to find bottom right
            self.minimap.move(self.view.width() - self.minimap.width() - margin, 
                              self.view.height() - self.minimap.height() - margin)

    def view_double_click(self, event):
        if event.button() == Qt.MouseButton.LeftButton:
            self.show_add_menu(pos=self.view.mapToGlobal(event.pos()))
        else:
            super(NodeView, self.view).mouseDoubleClickEvent(event)

    def show_add_menu(self, pos=None):
        if isinstance(pos, bool):
            pos = None
        menu = QMenu(self)
        menu.setStyleSheet("QMenu { background-color: #333; color: white; padding: 5px; }")
        
        from .node_types import GeneNode, CellNode, MathNode, TextNode, ParameterNode, StateNode, DerivativeNode, PythonCodeNode
        nodes = {
            "Gene (Legacy)": GeneNode,
            "Parameter": ParameterNode,
            "State (Y)": StateNode,
            "Derivative (dY)": DerivativeNode,
            "Math": MathNode,
            "Python Code": PythonCodeNode,
            "Cell": CellNode,
            "Text Note": TextNode
        }
        
        for name, cls in nodes.items():
            action = menu.addAction(name)
            action.triggered.connect(lambda checked, c=cls: self.add_node(c))
            
        if pos is None:
            pos = self.add_node_btn.mapToGlobal(self.add_node_btn.rect().bottomLeft())
        menu.exec(pos)

    def add_node(self, node_class):
        node = node_class()
        # Default spawn in center of view
        view_center = self.view.mapToScene(self.view.viewport().rect().center())
        node.setPos(view_center.x() - node.width/2, view_center.y() - node.height/2)
        self.scene.addItem(node)

    def apply_ode_text(self):
        reply = QMessageBox.question(
            self, "Apply ODE Text", 
            "This will replace your current node graph with parsed text. Continue?",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No
        )
        if reply == QMessageBox.StandardButton.Yes:
            if self.parent():
                ode_str = self.parent().ode_edit.toPlainText() if hasattr(self.parent(), 'ode_edit') else ""
                params_str = self.parent().params_edit.toPlainText() if hasattr(self.parent(), 'params_edit') else ""
                color_func_str = self.parent().color_func_edit.toPlainText() if hasattr(self.parent(), 'color_func_edit') else ""
                
                self.parser.parse_ode_string(ode_str, params_str, color_func_str)
                
                # Verification
                generated_ode = self.generator.generate_ode().strip()
                import re
                def normalize(s):
                    return re.sub(r'\s+', '', s) if s else ""
                
                gen_params = ""
                gen_color = ""
                for item in self.scene.items():
                    if type(item).__name__ == "OdeBeginNode":
                        gen_params = getattr(item, "params_str", "")
                        gen_color = getattr(item, "color_func_str", "")
                        break
                        
                mismatch = []
                if normalize(generated_ode) != normalize(ode_str):
                    mismatch.append("ODE Text")
                if normalize(gen_params) != normalize(params_str):
                    mismatch.append("Parameters")
                if normalize(gen_color) != normalize(color_func_str):
                    mismatch.append("Color Function")
                    
                if mismatch:
                    fields = ", ".join(mismatch)
                    QMessageBox.warning(self, "Verification Warning", f"The generated graph differs from the original text in: {fields}!\nSome structures might be missing or altered.")
                else:
                    if hasattr(self.parent(), 'status_bar'):
                        self.parent().status_bar.showMessage("Parsed ODE to Nodes (Verified matching).")

    def export_json(self):
        filename, _ = QFileDialog.getSaveFileName(self, "Export Node Graph", "", "JSON Files (*.json)")
        if filename:
            try:
                data = self.generator.serialize_scene()
                with open(filename, 'w') as f:
                    f.write(data)
                QMessageBox.information(self, "Success", "Node graph exported successfully.")
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Failed to export: {str(e)}")

    def import_json(self):
        filename, _ = QFileDialog.getOpenFileName(self, "Import Node Graph", "", "JSON Files (*.json)")
        if filename:
            try:
                with open(filename, 'r') as f:
                    data = f.read()
                self.generator.deserialize_scene(data)
                QMessageBox.information(self, "Success", "Node graph imported successfully.")
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Failed to import: {str(e)}")

    def generate_code(self):
        code = self.generator.generate_ode()
        if self.parent():
            self.parent().ode_edit.setPlainText(code)
            self.parent().status_bar.showMessage("ODE code updated from Node Editor")
        self.accept()

    def handle_shortcuts(self, event):
        if event.modifiers() & Qt.KeyboardModifier.ControlModifier:
            if event.key() == Qt.Key.Key_C:
                pass # TODO Copy
            elif event.key() == Qt.Key.Key_V:
                pass # TODO Paste
            elif event.key() == Qt.Key.Key_X:
                pass # TODO Cut
        elif event.key() == Qt.Key.Key_Delete:
            nodes_to_delete = [item for item in self.scene.selectedItems() if type(item).__name__ != "NodeEdge"]
            edges_to_delete = [item for item in self.scene.selectedItems() if type(item).__name__ == "NodeEdge"]
            
            for edge in edges_to_delete:
                if hasattr(edge, 'remove'):
                    edge.remove()
                
            for node in nodes_to_delete:
                seq_in_sources = []
                seq_out_dests = []
                connected_inputs = []
                connected_outputs = []
                
                if hasattr(node, "seq_in") and node.seq_in:
                    for edge in node.seq_in.edges:
                        if edge.source_socket: seq_in_sources.append(edge.source_socket)
                if hasattr(node, "seq_out") and node.seq_out:
                    for edge in node.seq_out.edges:
                        if edge.dest_socket: seq_out_dests.append(edge.dest_socket)
                        
                for sock in getattr(node, 'inputs', []):
                    for edge in sock.edges:
                        if edge.source_socket and edge.source_socket.node != node:
                            connected_inputs.append(edge.source_socket)
                for sock in getattr(node, 'outputs', []):
                    for edge in sock.edges:
                        if edge.dest_socket and edge.dest_socket.node != node:
                            connected_outputs.append(edge.dest_socket)

                for socket in getattr(node, 'inputs', []) + getattr(node, 'outputs', []):
                    for edge in list(getattr(socket, 'edges', [])):
                        if hasattr(edge, 'remove'): edge.remove()
                        
                if hasattr(node, "seq_in") and node.seq_in:
                    for edge in list(node.seq_in.edges):
                        if hasattr(edge, 'remove'): edge.remove()
                if hasattr(node, "seq_out") and node.seq_out:
                    for edge in list(node.seq_out.edges):
                        if hasattr(edge, 'remove'): edge.remove()
                        
                self.scene.removeItem(node)
                
                from .node_core import NodeEdge
                if seq_in_sources and seq_out_dests:
                    source = seq_in_sources[0]
                    dest = seq_out_dests[0]
                    new_edge = NodeEdge(source, dest, "Sequence")
                    source.add_edge(new_edge)
                    dest.add_edge(new_edge)
                    self.scene.addItem(new_edge)
                    new_edge.update_path()
                    
                if len(connected_inputs) == 1 and len(connected_outputs) == 1:
                    source = connected_inputs[0]
                    dest = connected_outputs[0]
                    if source.socket_type != "seq_in" and dest.socket_type != "seq_out":
                        new_edge = NodeEdge(source, dest, "Data")
                        source.add_edge(new_edge)
                        dest.add_edge(new_edge)
                        self.scene.addItem(new_edge)
                        new_edge.update_path()

        else:
            super().keyPressEvent(event)

    def closeEvent(self, event):
        # Auto-save layout mapping logic here (if integrating with config.json directly)
        super().closeEvent(event)
