from PyQt6.QtWidgets import (
    QGraphicsItem, QGraphicsRectItem, QGraphicsTextItem, QGraphicsProxyWidget,
    QWidget, QVBoxLayout, QLineEdit, QComboBox, QMenu, QFormLayout, QSpinBox, QDoubleSpinBox
)
from PyQt6.QtCore import Qt, QRectF, QPointF
from PyQt6.QtGui import QPen, QBrush, QColor, QFont, QPixmap, QIcon
from .node_core import NodeSocket

NODE_COLORS = [
    QColor(0, 0, 0, 0), # None
    QColor(255, 100, 100), QColor(100, 255, 100), QColor(100, 100, 255),
    QColor(255, 255, 100), QColor(255, 100, 255), QColor(100, 255, 255),
    QColor(200, 100, 100), QColor(100, 200, 100), QColor(100, 100, 200),
    QColor(200, 200, 100), QColor(200, 100, 200), QColor(100, 200, 200),
    QColor(150, 150, 150), QColor(80, 80, 80), QColor(220, 150, 50),
    QColor(50, 150, 220)
]

class BaseNode(QGraphicsItem):
    def __init__(self, title="Node"):
        super().__init__()
        self.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsMovable)
        self.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsSelectable)
        self.setFlag(QGraphicsItem.GraphicsItemFlag.ItemSendsGeometryChanges)
        
        self.title = title
        self.width = 180
        self.height = 100
        self.bg_color = QColor(40, 40, 40, 220)
        self.title_color = NODE_COLORS[13]
        self.bypassed = False
        
        self._is_resizing = False
        self._resize_start_pos = QPointF()
        self._start_width = 180
        self._start_height = 100
        
        self.inputs = []   # list of NodeSocket
        self.outputs = []  # list of NodeSocket
        self.seq_in = None
        self.seq_out = None
        
        self.title_item = QGraphicsTextItem(self.title, self)
        self.title_item.setDefaultTextColor(Qt.GlobalColor.white)
        font = QFont("Arial", 10, QFont.Weight.Bold)
        self.title_item.setFont(font)
        self.title_item.setPos(5, 5)

        # Widget container
        self.proxy = QGraphicsProxyWidget(self)
        self.widget_container = QWidget()
        self.widget_container.setAttribute(Qt.WidgetAttribute.WA_TranslucentBackground)
        self.widget_container.setStyleSheet("""
            QWidget { color: white; font-size: 10pt; }
            QLineEdit, QComboBox, QSpinBox, QDoubleSpinBox {
                background-color: #555;
                color: white;
                border: 1px solid #777;
                border-radius: 3px;
                padding: 2px;
            }
            QComboBox::drop-down { border: 0px; }
            QComboBox QAbstractItemView {
                background-color: #444;
                color: white;
                selection-background-color: #666;
            }
        """)
        self.layout = QVBoxLayout(self.widget_container)
        self.layout.setContentsMargins(5, 5, 5, 5)
        self.proxy.setWidget(self.widget_container)
        self.proxy.setPos(10, 30)
        
    def boundingRect(self):
        return QRectF(0, 0, self.width, self.height)

    def paint(self, painter, option, widget=None):
        rect = self.boundingRect()
        
        # Resize handle hit area
        resize_rect = QRectF(self.width - 20, self.height - 20, 20, 20)
        
        # Draw background
        painter.setBrush(QBrush(self.bg_color))
        pen = QPen(Qt.GlobalColor.white if self.isSelected() else Qt.GlobalColor.black, 2)
        if self.bypassed:
            pen.setStyle(Qt.PenStyle.DashLine)
            painter.setOpacity(0.5)
        painter.setPen(pen)
        painter.drawRoundedRect(rect, 5, 5)
        
        # Draw resize handle
        painter.setBrush(QBrush(QColor(100, 100, 100, 200)))
        painter.setPen(Qt.PenStyle.NoPen)
        # Draw a tiny triangle block
        pts = [QPointF(self.width, self.height - 15), QPointF(self.width - 15, self.height), QPointF(self.width, self.height)]
        from PyQt6.QtGui import QPolygonF
        painter.drawPolygon(QPolygonF(pts))
        
        # Draw title rect
        title_rect = QRectF(0, 0, self.width, 25)
        painter.setBrush(QBrush(self.title_color))
        painter.setPen(Qt.PenStyle.NoPen)
        # Top-rounded only, simplified
        painter.drawRoundedRect(title_rect, 5, 5)
        
        # Dim if bypassed
        if self.bypassed:
            painter.setBrush(QBrush(QColor(0, 0, 0, 100)))
            painter.drawRoundedRect(rect, 5, 5)
            painter.setOpacity(1.0)

    def add_input(self, name, value_type="Any"):
        socket = NodeSocket(self, "input", name, value_type)
        self.inputs.append(socket)
        self.update_sockets()
        return socket

    def add_output(self, name, value_type="Any"):
        socket = NodeSocket(self, "output", name, value_type)
        self.outputs.append(socket)
        self.update_sockets()
        return socket

    def compute_size(self):
        from PyQt6.QtGui import QFontMetrics
        fm = QFontMetrics(self.title_item.font())
        title_width = fm.horizontalAdvance(self.title) + 40
        
        num_inputs = len(self.inputs)
        num_outputs = len(self.outputs)
        slots_height = max(num_inputs, num_outputs) * 20
        
        widget_hint = self.widget_container.sizeHint()
        
        min_width = max(title_width, widget_hint.width() + 20, 120)
        min_height = 40 + slots_height + widget_hint.height() + 10
        
        return min_width, min_height

    def update_sockets(self):
        y_offset = 40
        num_inputs = len(self.inputs)
        num_outputs = len(self.outputs)
        spacing = 20
        
        for i, socket in enumerate(self.inputs):
            socket.setPos(0, y_offset + i * spacing)
            
        for i, socket in enumerate(self.outputs):
            socket.setPos(self.width, y_offset + i * spacing)
            
        if hasattr(self, 'seq_in') and self.seq_in:
            self.seq_in.setPos(self.width / 2, 0)
        if hasattr(self, 'seq_out') and self.seq_out:
            self.seq_out.setPos(self.width / 2, self.height)
            
        min_w, min_h = self.compute_size()
        
        if not self._is_resizing:
            if self.width < min_w:
                self.width = min_w
            if self.height < min_h:
                self.height = min_h
                
        self.proxy.resize(max(self.width - 20, 10), max(self.height - 40, 10))
        self.update()
        
        socks = self.inputs + self.outputs
        if hasattr(self, 'seq_in') and self.seq_in: socks.append(self.seq_in)
        if hasattr(self, 'seq_out') and self.seq_out: socks.append(self.seq_out)
        for socket in socks:
            for edge in getattr(socket, 'edges', []):
                if hasattr(edge, 'update_path'):
                    edge.update_path()

    def itemChange(self, change, value):
        if change == QGraphicsItem.GraphicsItemChange.ItemPositionHasChanged:
            socks = self.inputs + self.outputs
            if hasattr(self, 'seq_in') and self.seq_in: socks.append(self.seq_in)
            if hasattr(self, 'seq_out') and self.seq_out: socks.append(self.seq_out)
            for socket in socks:
                for edge in socket.edges:
                    edge.update_path()
        return super().itemChange(change, value)

    def mousePressEvent(self, event):
        resize_rect = QRectF(self.width - 20, self.height - 20, 20, 20)
        if resize_rect.contains(event.pos()) and event.button() == Qt.MouseButton.LeftButton:
            self._is_resizing = True
            self._resize_start_pos = event.scenePos()
            self._start_width = self.width
            self._start_height = self.height
            event.accept()
        else:
            super().mousePressEvent(event)

    def mouseMoveEvent(self, event):
        if self._is_resizing:
            delta = event.scenePos() - self._resize_start_pos
            min_w, min_h = self.compute_size()
            new_w = max(min_w, self._start_width + delta.x())
            new_h = max(min_h, self._start_height + delta.y())
            self.prepareGeometryChange()
            self.width = new_w
            self.height = new_h
            self.proxy.resize(max(self.width - 20, 10), max(self.height - 40, 10))
            self.update_sockets()
            self.update()
            event.accept()
        else:
            super().mouseMoveEvent(event)

    def mouseReleaseEvent(self, event):
        if self._is_resizing:
            self._is_resizing = False
            event.accept()
        else:
            super().mouseReleaseEvent(event)

    def contextMenuEvent(self, event):
        menu = QMenu()
        menu.setStyleSheet("QMenu { background-color: #333; color: white; }")
        
        bypass_action = menu.addAction("Bypass" if not self.bypassed else "Enable")
        menu.addSeparator()
        
        color_menu = menu.addMenu("Set Color")
        color_menu.setStyleSheet("QMenu { background-color: #333; color: white; }")
        for i, color in enumerate(NODE_COLORS):
            pix = QPixmap(16, 16)
            if color.alpha() == 0:
                pix.fill(Qt.GlobalColor.lightGray)
            else:
                pix.fill(color)
            name = "None" if color.alpha() == 0 else f"Color {i}"
            action = color_menu.addAction(QIcon(pix), name)
            action.triggered.connect(lambda checked, c=color: self.set_title_color(c))

        action = menu.exec(event.screenPos())
        if action == bypass_action:
            self.bypassed = not self.bypassed
            self.update()

    def set_title_color(self, color):
        self.title_color = color
        self.update()

    def serialize(self):
        return {
            "type": self.__class__.__name__,
            "x": self.pos().x(),
            "y": self.pos().y(),
            "width": self.width,
            "height": self.height,
            "bypassed": self.bypassed,
            "title_color": self.title_color.name()
        }

    def deserialize(self, data):
        self.setPos(data.get("x", 0), data.get("y", 0))
        if "width" in data and "height" in data:
            self.prepareGeometryChange()
            self.width = data["width"]
            self.height = data["height"]
            self.update_sockets()
        self.bypassed = data.get("bypassed", False)
        self.set_title_color(QColor(data.get("title_color", self.title_color.name())))


class SeqNode(BaseNode):
    def __init__(self, title="Seq Node"):
        super().__init__(title)
        self.seq_in = NodeSocket(self, "seq_in", "Sequence In", "Sequence")
        self.seq_out = NodeSocket(self, "seq_out", "Sequence Out", "Sequence")

class DataNode(BaseNode):
    def __init__(self, title="Data Node"):
        super().__init__(title)

class OdeBeginNode(SeqNode):
    def __init__(self):
        super().__init__("ODE begin")
        self.title_color = QColor(150, 50, 200) # Purple for special nodes
        self.width = 160
        self.height = 60
        # Only has seq_out, so we hide seq_in
        if self.seq_in:
            if self.scene(): self.scene().removeItem(self.seq_in)
            self.seq_in = None
        
        self.update_sockets()

class ForEachCellNode(SeqNode):
    def __init__(self):
        super().__init__("For each cell")
        self.title_color = QColor(150, 50, 200)
        self.width = 160
        self.height = 60
        self.add_output("i (cell index)", "Int")
        self.add_output("Cell Position", "Array")
        self.add_output("x_norm[i]", "Float")
        self.update_sockets()

class ForEachNeighborNode(SeqNode):
    def __init__(self):
        super().__init__("For each neighbor cell")
        self.title_color = QColor(150, 50, 200)
        self.width = 200
        self.height = 60
        self.update_sockets()

class EndOfForLoopNode(SeqNode):
    def __init__(self):
        super().__init__("End of For loop")
        self.title_color = QColor(150, 50, 200)
        self.width = 160
        self.height = 60
        self.update_sockets()

class ElseNode(SeqNode):
    def __init__(self):
        super().__init__("Else")
        self.title_color = QColor(150, 50, 200)
        self.width = 160
        self.height = 60
        self.update_sockets()

class GeneNode(DataNode):
    def __init__(self):
        super().__init__("Gene / Entity")
        self.title_color = NODE_COLORS[0]
        self.width = 260
        
        self.add_input("Modulators")
        self.add_output("Active (Out)")
        self.add_output("Inhibit (Out)")
        
        # Form Layout inside proxy
        form = QFormLayout()
        form.setContentsMargins(0,0,0,0)
        
        self.name_edit = QLineEdit("Gene_A")
        self.conc_edit = QLineEdit("0.0")
        self.dyn_edit = QLineEdit("1.0")
        self.deg_edit = QLineEdit("0.1")
        self.color_edit = QLineEdit("red")
        
        form.addRow("Name:", self.name_edit)
        form.addRow("Init Conc:", self.conc_edit)
        form.addRow("Produce:", self.dyn_edit)
        form.addRow("Degrade:", self.deg_edit)
        form.addRow("Color:", self.color_edit)
        
        self.layout.addLayout(form)

        self.update_sockets()

    def serialize(self):
        data = super().serialize()
        data["name"] = self.name_edit.text()
        data["conc"] = self.conc_edit.text()
        data["dyn"] = self.dyn_edit.text()
        data["deg"] = self.deg_edit.text()
        data["color"] = self.color_edit.text()
        return data

    def deserialize(self, data):
        super().deserialize(data)
        self.name_edit.setText(data.get("name", "Gene_A"))
        self.conc_edit.setText(data.get("conc", "0.0"))
        self.dyn_edit.setText(data.get("dyn", "1.0"))
        self.deg_edit.setText(data.get("deg", "0.1"))
        self.color_edit.setText(data.get("color", "red"))


class CellNode(DataNode):
    def __init__(self):
        super().__init__("Cell Node")
        self.title_color = NODE_COLORS[2]
        self.width = 160
        
        self.add_input("Death Signals")
        self.add_input("Div Signals")
        
        self.add_output("State")
        
        self.type_combo = QComboBox()
        self.type_combo.addItems(["Local (本體)", "Adjacent (相鄰)", "Range (範圍)"])
        self.layout.addWidget(self.type_combo)
        

        self.update_sockets()

    def serialize(self):
        data = super().serialize()
        data["cell_type_index"] = self.type_combo.currentIndex()
        return data

    def deserialize(self, data):
        super().deserialize(data)
        self.type_combo.setCurrentIndex(data.get("cell_type_index", 0))


class MathNode(DataNode):
    def __init__(self):
        super().__init__("Math")
        self.title_color = NODE_COLORS[5]
        self.width = 150
        
        self.add_input("A", "Any")
        self.add_input("B", "Any")
        self.add_output("Result", "Any")
        
        self.op_combo = QComboBox()
        self.op_combo.addItems([
            "Add (+)", "Subtract (-)", "Multiply (*)", "Divide (/)", "Power (**)", 
            "Sine (sin)", "Cosine (cos)", "Tangent (tan)", 
            "Exp (exp)", "Log (log)", "Square Root (sqrt)", "Absolute (abs)"
        ])
        self.layout.addWidget(self.op_combo)
        
        self.op_combo.currentIndexChanged.connect(self.on_op_changed)
        self.update_sockets()

    def on_op_changed(self, index):
        is_unary = index >= 5
        if is_unary:
            while len(self.inputs) > 1:
                sock = self.inputs.pop()
                for edge in list(sock.edges): edge.remove()
                if self.scene(): self.scene().removeItem(sock)
        else:
            while len(self.inputs) < 2:
                self.add_input(f"{'A' if len(self.inputs)==0 else 'B'}", "Any")
        self.update_sockets()

    def serialize(self):
        data = super().serialize()
        data["op_index"] = self.op_combo.currentIndex()
        return data

    def deserialize(self, data):
        super().deserialize(data)
        self.op_combo.setCurrentIndex(data.get("op_index", 0))

class ConstantNode(DataNode):
    def __init__(self):
        super().__init__("Constant")
        self.title_color = NODE_COLORS[6]
        self.width = 140
        self.add_output("Value", "Float")
        
        form = QFormLayout()
        form.setContentsMargins(0,0,0,0)
        from PyQt6.QtGui import QDoubleValidator
        self.value_edit = QLineEdit("0.0")
        self.value_edit.setValidator(QDoubleValidator())
        form.addRow("Val:", self.value_edit)
        self.layout.addLayout(form)

        self.update_sockets()

    def serialize(self):
        data = super().serialize()
        data["value"] = self.value_edit.text() or "0.0"
        return data
        
    def deserialize(self, data):
        super().deserialize(data)
        self.value_edit.setText(str(data.get("value", 0.0)))

class ParameterNode(DataNode):
    def __init__(self):
        super().__init__("Parameter")
        self.title_color = NODE_COLORS[1]
        self.width = 180
        self.add_output("Value", "Float")
        
        form = QFormLayout()
        form.setContentsMargins(0,0,0,0)
        self.param_name_edit = QLineEdit("betaD")
        self.default_value_edit = QLineEdit("0.0")
        form.addRow("Name:", self.param_name_edit)
        form.addRow("Default:", self.default_value_edit)
        self.layout.addLayout(form)
        self.access_style = "get"

        self.update_sockets()

    def serialize(self):
        data = super().serialize()
        data["param_name"] = self.param_name_edit.text()
        data["default_value"] = self.default_value_edit.text()
        data["access_style"] = getattr(self, "access_style", "get")
        return data
        
    def deserialize(self, data):
        super().deserialize(data)
        self.param_name_edit.setText(data.get("param_name", "betaD"))
        self.default_value_edit.setText(data.get("default_value", "0.0"))
        self.access_style = data.get("access_style", "get")

class StateNode(DataNode):
    def __init__(self):
        super().__init__("State (Y)")
        self.title_color = NODE_COLORS[2]
        self.width = 150
        self.add_output("Concentration", "Array")
        self.color_socket = self.add_output("→ Color Func", "ColorMapping")
        self.param_socket = self.add_output("→ Params Label", "ParamLabel")
        
        form = QFormLayout()
        form.setContentsMargins(0,0,0,0)
        self.index_edit = QSpinBox()
        self.index_edit.setMinimum(0)
        self.index_edit.setMaximum(99)
        form.addRow("Y Index:", self.index_edit)
        self.layout.addLayout(form)
        self.var_name = ""

        self.update_sockets()

    def serialize(self):
        data = super().serialize()
        data["state_idx"] = self.index_edit.value()
        data["var_name"] = getattr(self, "var_name", "")
        return data

    def deserialize(self, data):
        super().deserialize(data)
        self.index_edit.setValue(data.get("state_idx", 0))
        self.var_name = data.get("var_name", "")

class DerivativeNode(SeqNode):
    def __init__(self):
        super().__init__("Derivative (dY)")
        self.title_color = NODE_COLORS[3]
        self.bg_color = QColor(60, 20, 20, 220)
        self.width = 150
        self.add_input("dVar", "Array")
        
        form = QFormLayout()
        form.setContentsMargins(0,0,0,0)
        self.index_edit = QSpinBox()
        self.index_edit.setMinimum(0)
        self.index_edit.setMaximum(99)
        form.addRow("dY Index:", self.index_edit)
        self.layout.addLayout(form)

        self.update_sockets()

    def serialize(self):
        data = super().serialize()
        data["state_idx"] = self.index_edit.value()
        return data

    def deserialize(self, data):
        super().deserialize(data)
        self.index_edit.setValue(data.get("state_idx", 0))

class AssignmentNode(SeqNode):
    def __init__(self):
        super().__init__("Assign")
        self.title_color = NODE_COLORS[4]
        self.width = 180
        self.add_input("Value", "Any")
        self.add_output("Variable", "Any")

        form = QFormLayout()
        form.setContentsMargins(0,0,0,0)
        self.var_name_edit = QLineEdit("var")
        form.addRow("Name:", self.var_name_edit)
        self.layout.addLayout(form)
        self.value_mode = "scalar"

        self.update_sockets()

    def serialize(self):
        data = super().serialize()
        data["var_name"] = self.var_name_edit.text()
        data["value_mode"] = getattr(self, "value_mode", "scalar")
        return data

    def deserialize(self, data):
        super().deserialize(data)
        self.var_name_edit.setText(data.get("var_name", "var"))
        self.value_mode = data.get("value_mode", "scalar")

class PythonCodeNode(SeqNode):
    def __init__(self):
        super().__init__("Python Code")
        self.title_color = NODE_COLORS[8]
        self.width = 250
        self.height = 180
        
        from PyQt6.QtWidgets import QTextEdit
        
        self.code_edit = QTextEdit()
        self.code_edit.setPlaceholderText("out_0 = in_0 + in_1")
        self.code_edit.setStyleSheet("background-color: #2b2b2b; color: #a9b7c6; font-family: Consolas, monospace;")
        self.code_edit.setLineWrapMode(QTextEdit.LineWrapMode.NoWrap)
        self.layout.addWidget(self.code_edit)
        
        self.add_output("out_0", "Any")
        self.update_sockets()

    def compute_size(self):
        min_width, min_height = super().compute_size()
        doc_size = self.code_edit.document().size()
        return max(int(min_width), int(doc_size.width()) + 50), max(int(min_height), int(doc_size.height()) + 80)

    def add_dynamic_input(self):
        idx = len(self.inputs)
        self.add_input(f"in_{idx}", "Any")
        
    def rem_dynamic_input(self):
        if len(self.inputs) > 0:
            socket = self.inputs.pop()
            for edge in list(socket.edges):
                edge.remove()
            if self.scene():
                self.scene().removeItem(socket)
    
        self.update_sockets()

    def serialize(self):
        data = super().serialize()
        data["code"] = self.code_edit.toPlainText()
        data["num_inputs"] = len(self.inputs)
        data["output_names"] = [s.name for s in self.outputs]
        return data

    def deserialize(self, data):
        super().deserialize(data)
        self.code_edit.setPlainText(data.get("code", ""))
        num_in = data.get("num_inputs", 0)
        while len(self.inputs) < num_in:
            self.add_dynamic_input()
        out_names = data.get("output_names", ["out_0"])
        if len(out_names) > 0:
            self.outputs[0].name = out_names[0]
            for i in range(1, len(out_names)):
                self.add_output(out_names[i], "Any")
        self.update_sockets()


class TextNode(DataNode):
    def __init__(self):
        super().__init__("Text Note")
        self.title_color = QColor(0,0,0,0) # Invisible title
        self.bg_color = QColor(100, 100, 50, 150)
        self.width = 200
        self.height = 150
        self.setZValue(-10) # Always at bottom
        
        self.title_item.hide()
        
        from PyQt6.QtWidgets import QTextEdit
        self.text_edit = QTextEdit()
        self.text_edit.setPlaceholderText("Enter notes here...")
        self.text_edit.setStyleSheet("background-color: transparent; color: white; border: none;")
        self.layout.addWidget(self.text_edit)

        self.update_sockets()

    def serialize(self):
        data = super().serialize()
        data["text"] = self.text_edit.toPlainText()
        return data

    def deserialize(self, data):
        super().deserialize(data)
        self.text_edit.setPlainText(data.get("text", ""))

class IfNode(SeqNode):
    def __init__(self):
        super().__init__("If")
        self.title_color = QColor(150, 50, 200)
        self.width = 200
        self.height = 100
        
        self.condition_edit = QLineEdit()
        self.condition_edit.setPlaceholderText("t < 14.0")
        
        form = QFormLayout()
        form.setContentsMargins(0,0,0,0)
        form.addRow("Cond:", self.condition_edit)
        self.layout.addLayout(form)
        
        self.seq_out.name = "Then / Next"
        self.then_out = self.seq_out
        self.else_out = None
        
        self.add_input("Condition (Opt)", "Any")
        self.update_sockets()
        
    def update_sockets(self):
        super().update_sockets()
        if hasattr(self, 'then_out') and self.then_out and self.then_out != self.seq_out:
            self.then_out.setPos(self.width / 2 - 20, self.height)
            for e in self.then_out.edges: e.update_path()
        if hasattr(self, 'else_out') and self.else_out:
            self.else_out.setPos(self.width / 2 + 20, self.height)
            for e in self.else_out.edges: e.update_path()

    def serialize(self):
        data = super().serialize()
        data["condition"] = self.condition_edit.text()
        return data

    def deserialize(self, data):
        super().deserialize(data)
        self.condition_edit.setText(data.get("condition", ""))

class EndIfNode(SeqNode):
    def __init__(self):
        super().__init__("End If")
        self.title_color = QColor(150, 50, 200)
        self.width = 160
        self.height = 60
        self.update_sockets()

class WhileNode(SeqNode):
    def __init__(self):
        super().__init__("While")
        self.title_color = QColor(150, 50, 200)
        self.width = 200
        self.height = 100
        
        self.condition_edit = QLineEdit()
        self.condition_edit.setPlaceholderText("t < 14.0")
        
        form = QFormLayout()
        form.setContentsMargins(0,0,0,0)
        form.addRow("Cond:", self.condition_edit)
        self.layout.addLayout(form)
        
        self.add_input("Condition (Opt)", "Any")
        self.update_sockets()

    def serialize(self):
        data = super().serialize()
        data["condition"] = self.condition_edit.text()
        return data

    def deserialize(self, data):
        super().deserialize(data)
        self.condition_edit.setText(data.get("condition", ""))

class EndWhileNode(SeqNode):
    def __init__(self):
        super().__init__("End While")
        self.title_color = QColor(150, 50, 200)
        self.width = 160
        self.height = 60
        self.update_sockets()

class ForLoopNode(SeqNode):
    def __init__(self):
        super().__init__("For loop")
        self.title_color = QColor(150, 50, 200)
        self.width = 220
        self.height = 100

        form = QFormLayout()
        form.setContentsMargins(0,0,0,0)
        self.target_edit = QLineEdit("j")
        self.iter_edit = QLineEdit("range(10)")
        form.addRow("Target:", self.target_edit)
        form.addRow("Iter:", self.iter_edit)
        self.layout.addLayout(form)

        self.update_sockets()

    def serialize(self):
        data = super().serialize()
        data["target"] = self.target_edit.text()
        data["iter"] = self.iter_edit.text()
        return data

    def deserialize(self, data):
        super().deserialize(data)
        self.target_edit.setText(data.get("target", "j"))
        self.iter_edit.setText(data.get("iter", "range(10)"))

class OdeEndNode(SeqNode):
    def __init__(self):
        super().__init__("ODE End")
        self.title_color = QColor(150, 50, 200)
        self.width = 160
        self.height = 60
        if self.seq_out:
            if self.scene(): self.scene().removeItem(self.seq_out)
            self.seq_out = None
        self.update_sockets()

class VariableReadNode(DataNode):
    def __init__(self):
        super().__init__("Read Variable")
        self.title_color = NODE_COLORS[6]
        self.width = 160
        self.add_output("Value", "Any")
        
        self.var_combo = QComboBox()
        self.var_combo.setEditable(True)
        self.layout.addWidget(self.var_combo)
        self.update_sockets()

    def serialize(self):
        data = super().serialize()
        data["var_name"] = self.var_combo.currentText()
        return data

    def deserialize(self, data):
        super().deserialize(data)
        self.var_combo.setCurrentText(data.get("var_name", ""))

class NamedConstantNode(DataNode):
    def __init__(self):
        super().__init__("Named Constant")
        self.title_color = NODE_COLORS[6]
        self.width = 200
        self.add_output("Value", "Float")
        
        form = QFormLayout()
        form.setContentsMargins(0,0,0,0)
        self.name_edit = QLineEdit("omega_max")
        self.value_edit = QLineEdit("1.0")
        form.addRow("Name:", self.name_edit)
        form.addRow("Value:", self.value_edit)
        self.layout.addLayout(form)
        self.update_sockets()

    def serialize(self):
        data = super().serialize()
        data["const_name"] = self.name_edit.text()
        data["value"] = self.value_edit.text()
        return data

    def deserialize(self, data):
        super().deserialize(data)
        self.name_edit.setText(data.get("const_name", ""))
        self.value_edit.setText(data.get("value", ""))

class AggregateNode(DataNode):
    def __init__(self):
        super().__init__("Aggregate")
        self.title_color = NODE_COLORS[5]
        self.width = 160
        
        self.add_input("Values", "Array")
        self.add_input("Weights (Opt)", "Array")
        self.add_output("Result", "Any")
        
        self.op_combo = QComboBox()
        self.op_combo.addItems(["Sum", "Mean", "Max", "Min", "Count Non-Zero", "Weighted Sum"])
        self.layout.addWidget(self.op_combo)
        self.update_sockets()

    def serialize(self):
        data = super().serialize()
        data["op_index"] = self.op_combo.currentIndex()
        return data

    def deserialize(self, data):
        super().deserialize(data)
        self.op_combo.setCurrentIndex(data.get("op_index", 0))

class ContinueNode(SeqNode):
    def __init__(self):
        super().__init__("Continue")
        self.title_color = QColor(150, 50, 200)
        self.width = 140
        self.height = 60
        self.update_sockets()

class BreakNode(SeqNode):
    def __init__(self):
        super().__init__("Break")
        self.title_color = QColor(150, 50, 200)
        self.width = 140
        self.height = 60
        self.update_sockets()

class ColorMappingNode(DataNode):
    def __init__(self):
        super().__init__("Color Mapping")
        self.title_color = NODE_COLORS[15]
        self.width = 200
        
        self.add_input("State", "ColorMapping")
        self.add_output("Color Func", "Any")
        
        form = QFormLayout()
        form.setContentsMargins(0,0,0,0)
        self.transform_combo = QComboBox()
        self.transform_combo.addItems(["Linear", "(-cos(x)+1)/2", "Step"])
        self.colormap_combo = QComboBox()
        self.colormap_combo.addItems(["Reds", "Blues", "Viridis", "Plasma"])
        
        form.addRow("Transform:", self.transform_combo)
        form.addRow("Colormap:", self.colormap_combo)
        self.layout.addLayout(form)
        self.update_sockets()

    def serialize(self):
        data = super().serialize()
        data["transform"] = self.transform_combo.currentIndex()
        data["colormap"] = self.colormap_combo.currentIndex()
        return data

    def deserialize(self, data):
        super().deserialize(data)
        self.transform_combo.setCurrentIndex(data.get("transform", 0))
        self.colormap_combo.setCurrentIndex(data.get("colormap", 0))

