import ast
import json
from dataclasses import dataclass
from PyQt6.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QPushButton, QMenu, QFileDialog, QMessageBox, QInputDialog,
    QWidget, QListWidget, QListWidgetItem, QSplitter, QLabel, QFrame
)
from PyQt6.QtCore import Qt
from PyQt6.QtGui import QAction

from .node_core import NodeScene, NodeView, NodeEdge
from .node_types import GeneNode, CellNode, MathNode, TextNode, VariableReadNode, ParameterNode, NamedConstantNode, StateNode
from .node_parser import PythonToNodeParser, NodeToPythonGenerator
from PyQt6.QtGui import QPainter, QColor, QPen, QBrush
from PyQt6.QtCore import QRectF, QPointF, QTimer


@dataclass(frozen=True)
class NodeEditorDocument:
    """The complete model payload exchanged with the visual node editor."""

    ode_source: str = ""
    params_source: str = ""
    color_func_source: str = ""

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
        self.label.setObjectName("sidebarTitle")
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

    def __init__(self, parent=None, document=None):
        super().__init__(parent)
        self._source_document = document or NodeEditorDocument()
        self._result_document = None
        self.setWindowTitle("Visual Node Editor")
        self.resize(1300, 800)
        self.setMinimumSize(960, 640)
        
        main_layout = QVBoxLayout(self)
        main_layout.setContentsMargins(14, 12, 14, 12)
        main_layout.setSpacing(10)

        header = QHBoxLayout()
        header_copy = QVBoxLayout()
        title = QLabel("Visual Node Editor")
        title.setObjectName("editorTitle")
        subtitle = QLabel(
            "Build the ODE as a graph. Parameters and color mapping travel with the model."
        )
        subtitle.setObjectName("editorSubtitle")
        subtitle.setWordWrap(True)
        header_copy.addWidget(title)
        header_copy.addWidget(subtitle)
        header.addLayout(header_copy)
        header.addStretch(1)
        model_badge = QLabel("CURRENT MODEL")
        model_badge.setObjectName("modelBadge")
        header.addWidget(model_badge)
        main_layout.addLayout(header)
        
        toolbar_frame = QFrame()
        toolbar_frame.setObjectName("toolbarFrame")
        toolbar = QHBoxLayout(toolbar_frame)
        toolbar.setContentsMargins(8, 6, 8, 6)
        
        self.add_node_btn = QPushButton("Add Node...")
        self.add_node_btn.setObjectName("accentButton")
        self.apply_ode_btn = QPushButton("Reload Current Model")
        self.load_config_btn = QPushButton("Open Config...")
        self.export_btn = QPushButton("Export Graph")
        self.import_btn = QPushButton("Import Graph")
        
        self.routing_btn = QPushButton("Style: Spline")
        self.scene = NodeScene(self)
        self.scene.routing_style = "Spline"
        
        toolbar.addWidget(self.add_node_btn)
        toolbar.addWidget(self.apply_ode_btn)
        toolbar.addStretch(1)
        toolbar.addWidget(self.load_config_btn)
        toolbar.addWidget(self.import_btn)
        toolbar.addWidget(self.export_btn)
        toolbar.addWidget(self.routing_btn)
        main_layout.addWidget(toolbar_frame)
        
        # Splitter
        splitter = QSplitter(Qt.Orientation.Horizontal)
        splitter.setChildrenCollapsible(False)
        main_layout.addWidget(splitter, 1)
        
        # Scene & View
        self.scene.setSceneRect(-2000, -2000, 4000, 4000)
        self.view = NodeView(self.scene, self)
        
        self.view.setAcceptDrops(True)
        self.view.dragEnterEvent = self.view_dragEnterEvent
        self.view.dragMoveEvent = self.view_dragMoveEvent
        self.view.dropEvent = self.view_dropEvent
        
        splitter.addWidget(self.view)
        
        self.sidebar = VariableSidebar(self.scene, self)
        self.sidebar.setObjectName("sidePanel")
        self.sidebar.setMinimumWidth(210)
        splitter.addWidget(self.sidebar)
        splitter.setSizes([1000, 250])
        splitter.setStretchFactor(0, 1)
        splitter.setStretchFactor(1, 0)
        
        # Minimap
        self.minimap = MinimapWidget(self.view, self.view)

        footer = QHBoxLayout()
        self.status_label = QLabel()
        self.status_label.setObjectName("editorStatus")
        self.close_btn = QPushButton("Close")
        self.generate_code_btn = QPushButton("Apply to Model")
        self.generate_code_btn.setObjectName("primaryButton")
        footer.addWidget(self.status_label, 1)
        footer.addWidget(self.close_btn)
        footer.addWidget(self.generate_code_btn)
        main_layout.addLayout(footer)
        
        # Connections
        self.add_node_btn.clicked.connect(self.show_add_menu)
        self.apply_ode_btn.clicked.connect(self.apply_ode_text)
        self.load_config_btn.clicked.connect(self.load_config_template)
        self.export_btn.clicked.connect(self.export_json)
        self.import_btn.clicked.connect(self.import_json)
        self.generate_code_btn.clicked.connect(self.generate_code)
        self.routing_btn.clicked.connect(self.toggle_routing)
        self.close_btn.clicked.connect(self.reject)

        # Setup Parser/Generator
        self.parser = PythonToNodeParser(self.scene)
        self.generator = NodeToPythonGenerator(self.scene)
        
        # Override view double click and resize
        self.view.mouseDoubleClickEvent = self.view_double_click
        self.view.resizeEvent = self.view_resizeEvent
        
        self.keyPressEvent = self.handle_shortcuts

        self.add_node_btn.setToolTip("Add a computational or control-flow node")
        self.apply_ode_btn.setToolTip(
            "Discard graph edits and rebuild from the model that opened this window"
        )
        self.load_config_btn.setToolTip("Load ODE, parameters and color mapping from a config file")
        self.import_btn.setToolTip("Import a previously exported node graph")
        self.export_btn.setToolTip("Save this node graph without applying it to the model")
        self.routing_btn.setToolTip("Switch connection paths between spline and orthogonal")
        self.generate_code_btn.setToolTip(
            "Generate Python from this graph and return all model fields to the main window"
        )

        self._apply_styles()
        if self._source_document.ode_source.strip():
            self.load_document(self._source_document)
        else:
            self._set_status("Empty graph - add a node or open a config.", "info")

    def _apply_styles(self):
        self.setStyleSheet("""
            QDialog { background: #f3f5f7; color: #374151; }
            QLabel#editorTitle { font-size: 20px; font-weight: 650; color: #25313c; }
            QLabel#editorSubtitle { color: #52606d; }
            QLabel#modelBadge {
                color: #6D3BBF; background: #efe8fb; border: 1px solid #d9c9f2;
                border-radius: 9px; padding: 3px 8px; font-size: 9px; font-weight: 700;
            }
            QFrame#toolbarFrame {
                background: white; border: 1px solid #d8dee5; border-radius: 6px;
            }
            QPushButton {
                min-height: 28px; padding: 2px 11px; color: #374151;
                background: white; border: 1px solid #cbd5df; border-radius: 4px;
            }
            QPushButton:hover { background: #eef3f7; border-color: #9aa8b6; }
            QPushButton#accentButton {
                color: #1A6FB5; border-color: #1A6FB5; font-weight: 600;
            }
            QPushButton#accentButton:hover { background: #e8f2fb; }
            QPushButton#primaryButton {
                min-width: 145px; color: white; background: #1A6FB5;
                border-color: #1A6FB5; font-weight: 650;
            }
            QPushButton#primaryButton:hover { background: #145b96; }
            QWidget#sidePanel {
                background: white; border: 1px solid #d8dee5; border-radius: 5px;
            }
            QLabel#sidebarTitle { color: #25313c; font-weight: 650; padding: 8px 8px 3px 8px; }
            QListWidget {
                color: #374151; background: white; border: none; padding: 4px;
            }
            QListWidget::item { padding: 5px; border-radius: 3px; }
            QListWidget::item:selected { color: white; background: #1A6FB5; }
            QLabel#editorStatus { color: #52606d; padding-left: 2px; }
            QLabel#editorStatus[statusKind="success"] { color: #0D8F6B; }
            QLabel#editorStatus[statusKind="warning"] { color: #B45309; }
            QLabel#editorStatus[statusKind="error"] { color: #C4211C; }
            QMenu {
                color: #374151; background: white; border: 1px solid #cbd5df; padding: 4px;
            }
            QMenu::item { padding: 6px 22px 6px 10px; }
            QMenu::item:selected { color: white; background: #1A6FB5; }
        """)

    def _set_status(self, message, kind="info"):
        self.status_label.setText(message)
        self.status_label.setProperty("statusKind", kind)
        self.status_label.style().unpolish(self.status_label)
        self.status_label.style().polish(self.status_label)

    def result_document(self):
        return self._result_document

    def load_document(self, document):
        self.parser.parse_ode_string(
            document.ode_source,
            document.params_source,
            document.color_func_source,
        )
        self.scene.arrange_nodes()
        self._source_document = document
        node_count = sum(
            1 for item in self.scene.items() if hasattr(item, "serialize")
        )
        self._set_status(
            f"Loaded current model into {node_count} nodes. Graph edits are not applied yet.",
            "success",
        )

    def build_document(self):
        ode_source = self.generator.generate_ode()
        params_source = self._source_document.params_source
        color_func_source = self._source_document.color_func_source
        for item in self.scene.items():
            if type(item).__name__ == "OdeBeginNode":
                params_source = getattr(item, "params_str", params_source)
                color_func_source = getattr(item, "color_func_str", color_func_source)
                break
        return NodeEditorDocument(
            ode_source=ode_source,
            params_source=params_source,
            color_func_source=color_func_source,
        )

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
                with open(filename, 'r', encoding='utf-8') as f:
                    data = json.load(f)
                    
                ode_str = data.get('ode', '')
                if not ode_str:
                    ode_str = data.get('ode_code', '')
                params_value = data.get('params', '')
                import pprint
                params_str = "params = " + pprint.pformat(params_value, width=60) if isinstance(params_value, dict) else str(params_value)
                color_func_str = data.get('color_func', '')

                document = NodeEditorDocument(
                    ode_source=ode_str,
                    params_source=params_str,
                    color_func_source=color_func_str,
                )
                self.load_document(document)
                self._set_status("Config loaded. Review the graph, then apply it to the model.", "success")
            except Exception as e:
                self._set_status(f"Could not load config: {e}", "error")
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
        self._set_status(f"Connection style: {self.scene.routing_style}.", "info")

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
        
        from .node_types import (GeneNode, CellNode, MathNode, TextNode, ParameterNode, StateNode,
                                 DerivativeNode, PythonCodeNode, AssignmentNode, ForLoopNode)
        nodes = {
            "Gene (Legacy)": GeneNode,
            "Parameter": ParameterNode,
            "State (Y)": StateNode,
            "Derivative (dY)": DerivativeNode,
            "Assign": AssignmentNode,
            "For Loop": ForLoopNode,
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
        self._set_status(f"Added {type(node).__name__}. Apply when the graph is ready.", "info")

    def apply_ode_text(self):
        reply = QMessageBox.question(
            self, "Reload Current Model",
            "Discard graph edits and rebuild this graph from the model that opened the editor?",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No
        )
        if reply == QMessageBox.StandardButton.Yes:
            self.load_document(self._source_document)

    def export_json(self):
        filename, _ = QFileDialog.getSaveFileName(self, "Export Node Graph", "", "JSON Files (*.json)")
        if filename:
            try:
                data = self.generator.serialize_scene()
                with open(filename, 'w', encoding='utf-8') as f:
                    f.write(data)
                self._set_status("Node graph exported.", "success")
            except Exception as e:
                self._set_status(f"Could not export graph: {e}", "error")
                QMessageBox.critical(self, "Error", f"Failed to export: {str(e)}")

    def import_json(self):
        filename, _ = QFileDialog.getOpenFileName(self, "Import Node Graph", "", "JSON Files (*.json)")
        if filename:
            try:
                with open(filename, 'r', encoding='utf-8') as f:
                    data = f.read()
                self.generator.deserialize_scene(data)
                node_count = sum(
                    1 for item in self.scene.items() if hasattr(item, "serialize")
                )
                self._set_status(
                    f"Imported {node_count} nodes. Review the graph before applying.",
                    "success",
                )
            except Exception as e:
                self._set_status(f"Could not import graph: {e}", "error")
                QMessageBox.critical(self, "Error", f"Failed to import: {str(e)}")

    def generate_code(self):
        try:
            result = self.build_document()
            ast.parse(result.ode_source)
        except Exception as e:
            self._set_status(f"Generated ODE is invalid: {e}", "error")
            QMessageBox.critical(
                self,
                "Cannot Apply Model",
                f"The generated ODE is not valid Python:\n\n{e}",
            )
            return

        warnings = getattr(self.generator, "warnings", [])
        if warnings:
            preview = "\n".join(warnings[:5])
            if len(warnings) > 5:
                preview += f"\n...and {len(warnings) - 5} more."
            reply = QMessageBox.warning(
                self,
                "Generation Warnings",
                f"The graph generated code with warnings:\n\n{preview}\n\nApply anyway?",
                QMessageBox.StandardButton.Apply | QMessageBox.StandardButton.Cancel,
                QMessageBox.StandardButton.Cancel,
            )
            if reply != QMessageBox.StandardButton.Apply:
                self._set_status("Apply cancelled; resolve graph warnings first.", "warning")
                return

        self._result_document = result
        self._set_status("Model generated and ready for the main window.", "success")
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
