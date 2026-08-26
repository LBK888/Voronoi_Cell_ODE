import ast
import hashlib
import json
import textwrap

class PythonToNodeParser:
    def __init__(self, scene):
        self.scene = scene
        self.nodes = []
        self.param_registry = {}
        self.named_constant_registry = {}

    def parse_ode_string(self, ode_str, params_str, color_func_str=""):
        self.scene.clear()
        self.param_registry.clear()
        self.named_constant_registry.clear()

        from .node_types import (PythonCodeNode, ParameterNode, StateNode, DerivativeNode, MathNode,
                                 OdeBeginNode, ForEachCellNode, EndOfForLoopNode, ConstantNode,
                                 ElseNode, IfNode, EndIfNode, VariableReadNode, ContinueNode,
                                 BreakNode, OdeEndNode, WhileNode, EndWhileNode, AssignmentNode,
                                 ForLoopNode)
        from .node_core import NodeEdge

        try:
            tree = ast.parse(ode_str)
        except Exception as e:
            print("Failed to AST parse ODE:", e)
            fallback = PythonCodeNode()
            fallback.code_edit.setPlainText(ode_str)
            self.scene.addItem(fallback)
            return

        function_def = None
        for n in tree.body:
            if isinstance(n, ast.FunctionDef) and n.name == "ode":
                function_def = n
                break

        if not function_def:
            fallback = PythonCodeNode()
            fallback.code_edit.setPlainText(ode_str)
            self.scene.addItem(fallback)
            return

        source_lines = ode_str.splitlines(keepends=True)
        prefix_code = "".join(source_lines[:max(function_def.lineno - 1, 0)])
        suffix_code = "".join(source_lines[getattr(function_def, "end_lineno", len(source_lines)):])

        def source_for_stmt(stmt):
            src = ast.get_source_segment(ode_str, stmt)
            if src is None:
                src = ast.unparse(stmt)
            return textwrap.dedent(src).strip()

        core_body = []
        footer_lines = []
        in_footer = False

        for stmt in function_def.body:
            code = ast.unparse(stmt).strip()
            if code.startswith("dY = np.zeros_like") or code.startswith("p = params"):
                continue
            if code.startswith("Y_new = ") or code.startswith("dY = (Y_new") or code.startswith("return dY"):
                in_footer = True

            if in_footer:
                footer_lines.append(source_for_stmt(stmt))
            else:
                core_body.append(stmt)

        env = {}
        x_offset = 100
        y_offset = 100
        seq_counter = 0

        def mark_sequence_node(node):
            nonlocal seq_counter
            node.sequence_index = seq_counter
            seq_counter += 1
            return node

        def add_sequence_node(node):
            mark_sequence_node(node)
            self.scene.addItem(node)
            return node

        begin_node = OdeBeginNode()
        begin_node.setPos(x_offset, y_offset)
        begin_node.params_str = params_str
        begin_node.color_func_str = color_func_str
        begin_node.prefix_code = prefix_code
        begin_node.suffix_code = suffix_code
        begin_node.source_ode_str = ode_str
        if footer_lines:
            begin_node.footer_code = "\n".join("    " + line.replace("\n", "\n    ") for line in footer_lines) + "\n"
        add_sequence_node(begin_node)

        y_offset += 150

        def resolve_seq_socket(node, attr):
            sock = getattr(node, attr, None)
            if sock is None and attr == "seq_out":
                sock = getattr(node, "then_out", None)
            if sock is None and attr == "then_out":
                sock = getattr(node, "seq_out", None)
            return sock

        def link_seq(n1, n2, out_attr='seq_out', in_attr='seq_in'):
            if n1 and n2:
                out_sock = resolve_seq_socket(n1, out_attr)
                in_sock = getattr(n2, in_attr, None)
                if out_sock and in_sock:
                    edge = NodeEdge(out_sock, in_sock, "Sequence")
                    out_sock.add_edge(edge)
                    in_sock.add_edge(edge)
                    self.scene.addItem(edge)
                    edge.update_path()

        def commit_nodes(node_list):
            for n in node_list:
                if n.scene() != self.scene:
                    self.scene.addItem(n)

        def add_data_edge(source_socket, dest_socket):
            edge = NodeEdge(source_socket, dest_socket, "Data")
            source_socket.add_edge(edge)
            dest_socket.add_edge(edge)
            if source_socket.node.scene() == self.scene and dest_socket.node.scene() == self.scene:
                self.scene.addItem(edge)
            return edge

        def get_or_create_param(name, x, y, access_style="get", default_value="0.0"):
            key = (name, access_style, default_value if access_style == "get" else "")
            if key not in self.param_registry:
                p_node = ParameterNode()
                p_node.param_name_edit.setText(name)
                p_node.default_value_edit.setText(default_value)
                p_node.access_style = access_style
                p_node.setPos(x, y)
                self.scene.addItem(p_node)
                self.param_registry[key] = p_node
            return self.param_registry[key].outputs[0]

        def create_expr_fallback(expr, x, y):
            node = PythonCodeNode()
            node.code_edit.setPlainText(ast.unparse(expr).strip())
            node.setPos(x, y)
            node.is_fallback = True
            return node

        def create_variable_read(name, x, y):
            node = VariableReadNode()
            node.var_combo.setCurrentText(name)
            node.setPos(x, y)
            return node

        def subscript_elts(expr):
            if isinstance(expr, ast.Subscript) and isinstance(expr.slice, ast.Tuple):
                return expr.slice.elts
            return None

        def is_full_slice(node):
            return isinstance(node, ast.Slice) and node.lower is None and node.upper is None and node.step is None

        def y_column_index(expr):
            elts = subscript_elts(expr)
            if not elts or len(elts) != 2:
                return None
            if not (isinstance(expr.value, ast.Name) and expr.value.id == "Y"):
                return None
            if not isinstance(elts[1], ast.Constant):
                return None
            return elts[1].value

        def is_y_column_alias(expr):
            elts = subscript_elts(expr)
            return y_column_index(expr) is not None and len(elts) == 2 and is_full_slice(elts[0])

        def build_expr(expr, x_hints, y_hints):
            pending_nodes = []

            if isinstance(expr, ast.Name):
                if expr.id in env:
                    return env[expr.id], pending_nodes
                node = create_variable_read(expr.id, x_hints, y_hints)
                pending_nodes.append(node)
                return node.outputs[0], pending_nodes

            if isinstance(expr, ast.Constant):
                c_node = ConstantNode()
                c_node.value_edit.setText(repr(expr.value))
                c_node.setPos(x_hints, y_hints)
                pending_nodes.append(c_node)
                return c_node.outputs[0], pending_nodes

            if isinstance(expr, ast.UnaryOp) and isinstance(expr.op, ast.USub):
                if isinstance(expr.operand, ast.Constant) and isinstance(expr.operand.value, (int, float)):
                    c_node = ConstantNode()
                    c_node.value_edit.setText(repr(-expr.operand.value))
                    c_node.setPos(x_hints, y_hints)
                    pending_nodes.append(c_node)
                    return c_node.outputs[0], pending_nodes
                operand_socket, operand_pending = build_expr(expr.operand, x_hints - 200, y_hints)
                pending_nodes.extend(operand_pending)
                minus_one = ConstantNode()
                minus_one.value_edit.setText("-1")
                minus_one.setPos(x_hints - 200, y_hints + 100)
                pending_nodes.append(minus_one)
                m_node = MathNode()
                m_node.op_combo.setCurrentIndex(2)
                m_node.on_op_changed(2)
                m_node.setPos(x_hints, y_hints)
                pending_nodes.append(m_node)
                if operand_socket:
                    edge = NodeEdge(minus_one.outputs[0], m_node.inputs[0], "Data")
                    minus_one.outputs[0].add_edge(edge)
                    m_node.inputs[0].add_edge(edge)
                    pending_nodes.append(edge)
                    edge = NodeEdge(operand_socket, m_node.inputs[1], "Data")
                    operand_socket.add_edge(edge)
                    m_node.inputs[1].add_edge(edge)
                    pending_nodes.append(edge)
                    return m_node.outputs[0], pending_nodes

            if isinstance(expr, (ast.BoolOp, ast.Compare, ast.IfExp, ast.Lambda, ast.Dict, ast.List, ast.Tuple, ast.ListComp, ast.DictComp, ast.SetComp, ast.GeneratorExp)):
                node = create_expr_fallback(expr, x_hints, y_hints)
                pending_nodes.append(node)
                return node.outputs[0], pending_nodes

            if isinstance(expr, ast.Attribute):
                node = create_expr_fallback(expr, x_hints, y_hints)
                pending_nodes.append(node)
                return node.outputs[0], pending_nodes

            if isinstance(expr, ast.Call):
                if isinstance(expr.func, ast.Attribute) and getattr(expr.func.value, 'id', '') in ["params", "p"] and expr.func.attr == "get":
                    if len(expr.args) >= 1 and isinstance(expr.args[0], ast.Constant):
                        target_name = str(expr.args[0].value)
                        default_value = ast.unparse(expr.args[1]).strip() if len(expr.args) >= 2 else "0.0"
                        sock = get_or_create_param(target_name, x_hints, y_hints, "get", default_value)
                        return sock, pending_nodes

                func_name = None
                if isinstance(expr.func, ast.Attribute):
                    func_name = expr.func.attr
                elif isinstance(expr.func, ast.Name):
                    func_name = expr.func.id
                unary_map = {'sin': 5, 'cos': 6, 'tan': 7, 'exp': 8, 'log': 9, 'sqrt': 10, 'abs': 11}
                if func_name in unary_map and len(expr.args) >= 1:
                    in_socket, sub_pending = build_expr(expr.args[0], x_hints - 150, y_hints)
                    pending_nodes.extend(sub_pending)
                    if in_socket:
                        m_node = MathNode()
                        op_idx = unary_map[func_name]
                        m_node.op_combo.setCurrentIndex(op_idx)
                        m_node.on_op_changed(op_idx)
                        m_node.setPos(x_hints, y_hints)
                        pending_nodes.append(m_node)
                        edge = NodeEdge(in_socket, m_node.inputs[0], "Data")
                        in_socket.add_edge(edge)
                        m_node.inputs[0].add_edge(edge)
                        pending_nodes.append(edge)
                        return m_node.outputs[0], pending_nodes

                node = create_expr_fallback(expr, x_hints, y_hints)
                pending_nodes.append(node)
                return node.outputs[0], pending_nodes

            if isinstance(expr, ast.Subscript):
                slice_node = expr.slice
                if isinstance(expr.value, ast.Name) and expr.value.id in ["params", "p"]:
                    if isinstance(slice_node, ast.Constant):
                        target_name = str(slice_node.value)
                        sock = get_or_create_param(target_name, x_hints, y_hints, "subscript", "")
                        return sock, pending_nodes

                if isinstance(expr.value, ast.Name) and expr.value.id in env:
                    if isinstance(slice_node, ast.Name) and slice_node.id == "i":
                        sock = env[expr.value.id]
                        if hasattr(sock, "node") and type(sock.node).__name__ == "StateNode":
                            return sock, pending_nodes

                idx = y_column_index(expr)
                if idx is not None:
                    s_node = StateNode()
                    s_node.index_edit.setValue(idx)
                    s_node.setPos(x_hints, y_hints)
                    pending_nodes.append(s_node)
                    return s_node.outputs[0], pending_nodes

                node = create_expr_fallback(expr, x_hints, y_hints)
                pending_nodes.append(node)
                return node.outputs[0], pending_nodes

            if isinstance(expr, ast.BinOp):
                left_socket, left_pending = build_expr(expr.left, x_hints - 220, y_hints - 100)
                right_socket, right_pending = build_expr(expr.right, x_hints - 220, y_hints + 100)
                pending_nodes.extend(left_pending)
                pending_nodes.extend(right_pending)

                if not left_socket or not right_socket:
                    node = create_expr_fallback(expr, x_hints, y_hints)
                    pending_nodes.append(node)
                    return node.outputs[0], pending_nodes

                m_node = MathNode()
                op_map = {ast.Add: 0, ast.Sub: 1, ast.Mult: 2, ast.Div: 3, ast.Pow: 4}
                op_idx = 0
                for op_type, idx in op_map.items():
                    if isinstance(expr.op, op_type):
                        op_idx = idx
                        break
                m_node.op_combo.setCurrentIndex(op_idx)
                m_node.on_op_changed(op_idx)
                m_node.setPos(x_hints, y_hints)
                pending_nodes.append(m_node)

                edge = NodeEdge(left_socket, m_node.inputs[0], "Data")
                left_socket.add_edge(edge)
                m_node.inputs[0].add_edge(edge)
                pending_nodes.append(edge)
                edge = NodeEdge(right_socket, m_node.inputs[1], "Data")
                right_socket.add_edge(edge)
                m_node.inputs[1].add_edge(edge)
                pending_nodes.append(edge)
                return m_node.outputs[0], pending_nodes

            node = create_expr_fallback(expr, x_hints, y_hints)
            pending_nodes.append(node)
            return node.outputs[0], pending_nodes

        def process_stmts(stmts, current_prev, first_out_attr='seq_out'):
            nonlocal y_offset
            unhandled_stmts = []
            next_out_attr = first_out_attr

            def assigned_names(assign_stmt):
                names = []
                for target in getattr(assign_stmt, "targets", []):
                    if isinstance(target, ast.Name):
                        names.append(target.id)
                    elif isinstance(target, ast.Tuple):
                        names.extend(elt.id for elt in target.elts if isinstance(elt, ast.Name))
                return names

            def connect_next(cp, node):
                nonlocal next_out_attr
                link_seq(cp, node, out_attr=next_out_attr)
                next_out_attr = 'seq_out'

            def flush_unhandled(cp):
                nonlocal y_offset
                if not unhandled_stmts:
                    return cp
                lines = [source_for_stmt(s) for s in unhandled_stmts]
                code = "\n".join(line for line in lines if line)
                node = PythonCodeNode()
                node.code_edit.setPlainText(code)
                node.setPos(x_offset, y_offset)
                add_sequence_node(node)
                connect_next(cp, node)
                y_offset += 150

                for stmt in unhandled_stmts:
                    for sub_node in ast.walk(stmt):
                        if isinstance(sub_node, ast.Assign):
                            for var_name in assigned_names(sub_node):
                                if var_name in env:
                                    continue
                                if len(node.outputs) == 1 and node.outputs[0].name == "out_0":
                                    node.outputs[0].name = var_name
                                    env[var_name] = node.outputs[0]
                                else:
                                    sock = node.add_output(var_name, "Any")
                                    env[var_name] = sock
                        if isinstance(sub_node, ast.AugAssign) and isinstance(sub_node.target, ast.Name):
                            var_name = sub_node.target.id
                            if var_name not in env:
                                if len(node.outputs) == 1 and node.outputs[0].name == "out_0":
                                    node.outputs[0].name = var_name
                                    env[var_name] = node.outputs[0]
                                else:
                                    sock = node.add_output(var_name, "Any")
                                    env[var_name] = sock

                used_names = set()
                for stmt in unhandled_stmts:
                    for sub_node in ast.walk(stmt):
                        if isinstance(sub_node, ast.Name) and isinstance(sub_node.ctx, ast.Load):
                            used_names.add(sub_node.id)

                for var_name in used_names:
                    if var_name in env:
                        sock = env[var_name]
                        if sock and hasattr(sock, 'node') and sock.node != node:
                            in_sock = None
                            for s in getattr(node, "inputs", []):
                                if s.name == var_name:
                                    in_sock = s
                                    break
                            if not in_sock and hasattr(node, "add_input"):
                                in_sock = node.add_input(var_name, "Any")
                            if in_sock:
                                edge = NodeEdge(sock, in_sock, "Data")
                                sock.add_edge(edge)
                                in_sock.add_edge(edge)
                                self.scene.addItem(edge)

                unhandled_stmts.clear()
                return node

            for stmt in stmts:
                handled = False

                if isinstance(stmt, ast.Continue):
                    current_prev = flush_unhandled(current_prev)
                    node = ContinueNode()
                    node.setPos(x_offset, y_offset)
                    add_sequence_node(node)
                    connect_next(current_prev, node)
                    current_prev = node
                    y_offset += 150
                    continue

                if isinstance(stmt, ast.Break):
                    current_prev = flush_unhandled(current_prev)
                    node = BreakNode()
                    node.setPos(x_offset, y_offset)
                    add_sequence_node(node)
                    connect_next(current_prev, node)
                    current_prev = node
                    y_offset += 150
                    continue

                if isinstance(stmt, ast.For):
                    current_prev = flush_unhandled(current_prev)
                    is_cell_loop = (
                        isinstance(stmt.target, ast.Name)
                        and stmt.target.id == 'i'
                        and isinstance(stmt.iter, ast.Call)
                        and ast.unparse(stmt.iter).strip() == "range(Y.shape[0])"
                    )
                    if is_cell_loop:
                        node = ForEachCellNode()
                    else:
                        node = ForLoopNode()
                        node.target_edit.setText(ast.unparse(stmt.target).strip())
                        node.iter_edit.setText(ast.unparse(stmt.iter).strip())
                    node.setPos(x_offset, y_offset)
                    add_sequence_node(node)
                    connect_next(current_prev, node)
                    current_prev = node
                    y_offset += 150

                    current_prev = process_stmts(stmt.body, current_prev)

                    end_node = EndOfForLoopNode()
                    end_node.setPos(x_offset, y_offset)
                    add_sequence_node(end_node)
                    link_seq(current_prev, end_node)
                    current_prev = end_node
                    y_offset += 150
                    handled = True

                elif isinstance(stmt, ast.If):
                    current_prev = flush_unhandled(current_prev)
                    cond_str = ast.unparse(stmt.test).strip()

                    node = IfNode()
                    node.condition_edit.setText(cond_str)
                    node.setPos(x_offset, y_offset)
                    add_sequence_node(node)
                    connect_next(current_prev, node)
                    y_offset += 150

                    then_prev = process_stmts(stmt.body, node, first_out_attr='then_out')

                    if stmt.orelse:
                        else_node = ElseNode()
                        else_node.setPos(x_offset, y_offset)
                        add_sequence_node(else_node)
                        link_seq(then_prev, else_node)
                        y_offset += 150
                        else_prev = process_stmts(stmt.orelse, else_node)
                    else:
                        else_prev = then_prev

                    end_if = EndIfNode()
                    end_if.setPos(x_offset, y_offset)
                    add_sequence_node(end_if)
                    link_seq(else_prev, end_if)
                    current_prev = end_if
                    y_offset += 150
                    handled = True

                elif isinstance(stmt, ast.While):
                    current_prev = flush_unhandled(current_prev)
                    cond_str = ast.unparse(stmt.test).strip()

                    node = WhileNode()
                    node.condition_edit.setText(cond_str)
                    node.setPos(x_offset, y_offset)
                    add_sequence_node(node)
                    connect_next(current_prev, node)
                    y_offset += 150

                    body_prev = process_stmts(stmt.body, node)

                    end_while = EndWhileNode()
                    end_while.setPos(x_offset, y_offset)
                    add_sequence_node(end_while)
                    link_seq(body_prev, end_while)

                    current_prev = end_while
                    y_offset += 150
                    handled = True

                elif isinstance(stmt, ast.Assign):
                    if len(stmt.targets) == 1 and isinstance(stmt.targets[0], ast.Subscript):
                        tgt = stmt.targets[0]
                        if isinstance(tgt.value, ast.Name) and tgt.value.id == "dY":
                            slice_elts = subscript_elts(tgt)
                            if slice_elts and len(slice_elts) == 2 and isinstance(slice_elts[1], ast.Constant):
                                current_prev = flush_unhandled(current_prev)
                                out_socket, pending = build_expr(stmt.value, x_offset + 500, y_offset)
                                if out_socket:
                                    commit_nodes(pending)
                                    node = DerivativeNode()
                                    node.index_edit.setValue(slice_elts[1].value)
                                    node.setPos(x_offset + 200, y_offset)
                                    add_sequence_node(node)
                                    connect_next(current_prev, node)
                                    current_prev = node
                                    edge = NodeEdge(out_socket, node.inputs[0], "Data")
                                    out_socket.add_edge(edge)
                                    node.inputs[0].add_edge(edge)
                                    self.scene.addItem(edge)
                                    y_offset += 150
                                    handled = True

                    elif len(stmt.targets) == 1 and isinstance(stmt.targets[0], ast.Tuple):
                        tgt_tuple = stmt.targets[0]
                        val_tuple = stmt.value
                        if isinstance(val_tuple, ast.Tuple) and len(tgt_tuple.elts) == len(val_tuple.elts):
                            all_pending = []
                            for t_elt, v_elt in zip(tgt_tuple.elts, val_tuple.elts):
                                if isinstance(t_elt, ast.Name) and is_y_column_alias(v_elt):
                                    out_sock, pending = build_expr(v_elt, x_offset + 260, y_offset)
                                    if out_sock:
                                        if type(out_sock.node).__name__ == "StateNode":
                                            out_sock.node.var_name = t_elt.id
                                        env[t_elt.id] = out_sock
                                        all_pending.extend(pending)
                            commit_nodes(all_pending)
                        handled = False

                    elif len(stmt.targets) == 1 and isinstance(stmt.targets[0], ast.Name):
                        var_name = stmt.targets[0].id
                        out_socket, pending = build_expr(stmt.value, x_offset + 420, y_offset)
                        if out_socket:
                            current_prev = flush_unhandled(current_prev)
                            commit_nodes(pending)
                            node = AssignmentNode()
                            node.var_name_edit.setText(var_name)
                            if is_y_column_alias(stmt.value):
                                node.value_mode = "array"
                            node.setPos(x_offset + 180, y_offset)
                            add_sequence_node(node)
                            connect_next(current_prev, node)
                            current_prev = node

                            edge = NodeEdge(out_socket, node.inputs[0], "Data")
                            out_socket.add_edge(edge)
                            node.inputs[0].add_edge(edge)
                            self.scene.addItem(edge)

                            if type(out_socket.node).__name__ == "StateNode" and is_y_column_alias(stmt.value):
                                out_socket.node.var_name = var_name
                                env[var_name] = out_socket
                            else:
                                env[var_name] = node.outputs[0]
                            y_offset += 150
                            handled = True

                if not handled:
                    unhandled_stmts.append(stmt)

            current_prev = flush_unhandled(current_prev)
            return current_prev

        last_node = process_stmts(core_body, begin_node)

        end_node = OdeEndNode()
        end_node.setPos(x_offset, y_offset + 100)
        add_sequence_node(end_node)
        link_seq(last_node, end_node)

        for item in self.scene.items():
            if hasattr(item, "code_edit"):
                item.code_edit.document().adjustSize()
            if hasattr(item, "update_sockets"):
                item.update_sockets()

        try:
            begin_node.source_fingerprint = NodeToPythonGenerator(self.scene)._scene_semantic_fingerprint()
        except Exception:
            begin_node.source_fingerprint = ""

        self.scene.arrange_nodes()
        if self.scene.views():
            from PyQt6.QtCore import Qt
            self.scene.views()[0].fitInView(self.scene.itemsBoundingRect(), Qt.AspectRatioMode.KeepAspectRatio)


class NodeToPythonGenerator:
    def __init__(self, scene):
        self.scene = scene
        self.warnings = []

    def _scene_semantic_fingerprint(self):
        from .node_core import NodeEdge

        nodes_list = [item for item in self.scene.items() if hasattr(item, "serialize")]

        def node_sort_key(node):
            seq = getattr(node, "sequence_index", None)
            return (0 if seq is not None else 1, seq if seq is not None else 0, type(node).__name__, id(node))

        nodes_list = sorted(nodes_list, key=node_sort_key)
        node_to_id = {node: idx for idx, node in enumerate(nodes_list)}

        node_payload = []
        for node in nodes_list:
            data = node.serialize()
            for key in ("x", "y", "width", "height", "title_color"):
                data.pop(key, None)
            if hasattr(node, "sequence_index"):
                data["sequence_index"] = getattr(node, "sequence_index")
            if type(node).__name__ == "OdeBeginNode":
                data["footer_code"] = getattr(node, "footer_code", "")
                data["prefix_code"] = getattr(node, "prefix_code", "")
                data["suffix_code"] = getattr(node, "suffix_code", "")
            data.pop("source_ode_str", None)
            data.pop("source_fingerprint", None)
            node_payload.append(data)

        edge_payload = []
        for edge in [item for item in self.scene.items() if isinstance(item, NodeEdge)]:
            if not edge.source_socket or not edge.dest_socket:
                continue
            src_node = edge.source_socket.node
            dst_node = edge.dest_socket.node
            if src_node not in node_to_id or dst_node not in node_to_id:
                continue
            edge_payload.append({
                "source_node": node_to_id[src_node],
                "source_socket": edge.source_socket.name,
                "source_type": edge.source_socket.socket_type,
                "dest_node": node_to_id[dst_node],
                "dest_socket": edge.dest_socket.name,
                "dest_type": edge.dest_socket.socket_type,
                "edge_type": edge.edge_type
            })

        payload = {
            "nodes": node_payload,
            "edges": sorted(edge_payload, key=lambda e: json.dumps(e, sort_keys=True))
        }
        raw = json.dumps(payload, sort_keys=True, ensure_ascii=False)
        return hashlib.sha256(raw.encode("utf-8")).hexdigest()

    def _get_begin_node(self):
        for item in self.scene.items():
            if type(item).__name__ == "OdeBeginNode":
                return item
        return None

    def collect_main_chain(self):
        ordered = [
            item for item in self.scene.items()
            if hasattr(item, "sequence_index") and hasattr(item, "serialize") and not getattr(item, 'bypassed', False)
        ]
        if ordered:
            return sorted(ordered, key=lambda n: getattr(n, "sequence_index"))

        chain = []
        visited = set()
        current_node = self._get_begin_node()
        
        def traverse(node):
            if not node or node in visited:
                return
            if getattr(node, 'bypassed', False):
                return
            visited.add(node)
            chain.append(node)
            
            if type(node).__name__ == "OdeEndNode":
                return

            seq_socket = getattr(node, 'seq_out', None) or getattr(node, 'then_out', None)
            if seq_socket:
                for edge in seq_socket.edges:
                    if edge.edge_type != "Sequence":
                        continue
                    out_sock = edge.dest_socket if edge.source_socket == seq_socket else edge.source_socket
                    if out_sock and out_sock.node != node:
                        traverse(out_sock.node)
            else:
                for edge in getattr(getattr(node, 'then_out', None), "edges", []):
                    if edge.edge_type != "Sequence":
                        continue
                    out_sock = edge.dest_socket if edge.source_socket == getattr(node, 'then_out', None) else edge.source_socket
                    if out_sock and out_sock.node != node:
                        traverse(out_sock.node)
        traverse(current_node)
        return chain

    def generate_ode(self, preserve_roundtrip=True):
        self.warnings.clear()
        begin_for_roundtrip = self._get_begin_node()
        if preserve_roundtrip and begin_for_roundtrip:
            source = getattr(begin_for_roundtrip, "source_ode_str", "")
            fingerprint = getattr(begin_for_roundtrip, "source_fingerprint", "")
            if source and fingerprint and fingerprint == self._scene_semantic_fingerprint():
                return source

        main_chain = self.collect_main_chain()
        
        code = ""
        indent_level = 1
        vars_map = {}
        var_counter = 0
        generating_set = set()

        def get_input_var(socket, current_indent, usage="scalar"):
            nonlocal code, var_counter
            if not socket.edges:
                node_name = type(socket.node).__name__
                self.warnings.append(f"Warning: Node '{node_name}' input '{socket.name}' is unconnected, defaulting to 0.0")
                return "0.0"
                
            for edge in socket.edges:
                if edge.edge_type == "Sequence": continue
                out_sock = edge.source_socket if edge.source_socket and edge.source_socket.socket_type == "output" else (edge.dest_socket if edge.dest_socket and edge.dest_socket.socket_type == "output" else None)
                if out_sock and hasattr(out_sock, 'node'):
                    node = out_sock.node
                    if type(node).__name__ == "PythonCodeNode" and node in main_chain:
                        return out_sock.name
                    if type(node).__name__ == "StateNode":
                        return generate_data_node(node, current_indent, usage)
                    if out_sock in vars_map:
                        return vars_map[out_sock]
                    v_name = generate_data_node(node, current_indent, usage)
                    if v_name is not None:
                        return v_name
            return "0.0"

        def generate_data_node(node, current_indent, usage="scalar"):
            nonlocal code, var_counter
            if getattr(node, 'bypassed', False): return "0.0"
            if node in generating_set:
                self.warnings.append(f"Warning: Cyclic dependency detected at {type(node).__name__}")
                return "0.0"
                
            for sock in getattr(node, 'outputs', []):
                if sock in vars_map:
                    return vars_map[sock]
                    
            generating_set.add(node)
            indent = "    " * current_indent
            result_expr = "0.0"
            
            node_type = type(node).__name__
            
            if node_type == "ParameterNode":
                p_name = node.param_name_edit.text()
                access = getattr(node, 'access_style', 'get')
                default_val = getattr(node, 'default_value_edit', None)
                default_val_str = default_val.text() if default_val else "0.0"
                if access == 'get':
                    result_expr = f"p.get('{p_name}', {default_val_str})"
                else:
                    result_expr = f"p['{p_name}']"
                
            elif node_type == "StateNode":
                idx = node.index_edit.value()
                var_name = getattr(node, "var_name", "")
                if usage == "array":
                    result_expr = f"Y[:, {idx}]"
                elif var_name:
                    result_expr = f"{var_name}[i]"
                else:
                    result_expr = f"Y[i, {idx}]"
                generating_set.discard(node)
                return result_expr
                
            elif node_type == "MathNode":
                a_var = get_input_var(node.inputs[0], current_indent) if len(node.inputs) > 0 else "0.0"
                b_var = get_input_var(node.inputs[1], current_indent) if len(node.inputs) > 1 else "0.0"
                op_idx = node.op_combo.currentIndex()
                op_map = {0: f"({a_var} + {b_var})", 1: f"({a_var} - {b_var})", 2: f"({a_var} * {b_var})", 
                          3: f"({a_var} / {b_var})", 4: f"({a_var} ** {b_var})",
                          5: f"np.sin({a_var})", 6: f"np.cos({a_var})", 7: f"np.tan({a_var})",
                          8: f"np.exp({a_var})", 9: f"np.log({a_var})", 10: f"np.sqrt({a_var})", 11: f"np.abs({a_var})"}
                result_expr = op_map.get(op_idx, f"({a_var} + {b_var})")

            elif node_type == "ConstantNode":
                val = node.value_edit.text()
                result_expr = str(val) if val else "0.0"

            elif node_type == "PythonCodeNode":
                code_str = node.code_edit.toPlainText().strip()
                if node in main_chain:
                    for line in code_str.split("\n"):
                        if line.strip():
                            code += f"{indent}{line}\n"
                    result_expr = ""
                else:
                    result_expr = f"({code_str})"

            elif node_type == "AssignmentNode":
                result_expr = node.var_name_edit.text() or "var"
                    
            elif node_type == "VariableReadNode":
                result_expr = node.var_combo.currentText() or "0.0"
                
            elif node_type == "NamedConstantNode":
                result_expr = node.name_edit.text() or "0.0"
                
            elif node_type == "AggregateNode":
                val_arr = get_input_var(node.inputs[0], current_indent)
                w_arr = get_input_var(node.inputs[1], current_indent) if len(node.inputs) > 1 else "1.0"
                op_idx = node.op_combo.currentIndex()
                if op_idx == 0: result_expr = f"np.sum({val_arr})"
                elif op_idx == 1: result_expr = f"np.mean({val_arr})"
                elif op_idx == 2: result_expr = f"np.max({val_arr})"
                elif op_idx == 3: result_expr = f"np.min({val_arr})"
                elif op_idx == 4: result_expr = f"np.count_nonzero({val_arr})"
                elif op_idx == 5: result_expr = f"(np.sum({val_arr} * {w_arr}) / np.sum({w_arr}))"
            
            generating_set.discard(node)
            
            for sock in getattr(node, 'outputs', []):
                vars_map[sock] = result_expr
                
            return result_expr

        begin_node = None
        for node in main_chain:
            if type(node).__name__ == "OdeBeginNode":
                begin_node = node
                code += getattr(node, "prefix_code", "")
                code += "def ode(Y, t, params, cell_positions):\n"
                code += "    p = params\n"
                code += "    dY = np.zeros_like(Y)\n"
                
            elif type(node).__name__ == "ForEachCellNode":
                code += f"{'    ' * indent_level}for i in range(Y.shape[0]):\n"
                indent_level += 1

            elif type(node).__name__ == "ForLoopNode":
                target = node.target_edit.text() or "j"
                iterable = node.iter_edit.text() or "range(10)"
                code += f"{'    ' * indent_level}for {target} in {iterable}:\n"
                indent_level += 1
                
            elif type(node).__name__ == "ForEachNeighborNode":
                code += f"{'    ' * indent_level}neibs = neighbors[i]\n"
                code += f"{'    ' * indent_level}if neibs:\n"
                indent_level += 1
                
            elif type(node).__name__ == "IfNode":
                cond = node.condition_edit.text() or "True"
                code += f"{'    ' * indent_level}if {cond}:\n"
                indent_level += 1
                
            elif type(node).__name__ == "ElseNode":
                indent_level = max(1, indent_level - 1)
                if code.strip().endswith(':'):
                    code += f"{'    ' * indent_level}pass\n"
                code += f"{'    ' * indent_level}else:\n"
                indent_level += 1
                
            elif type(node).__name__ == "EndIfNode" or type(node).__name__ == "EndOfForLoopNode" or type(node).__name__ == "EndWhileNode":
                indent_level = max(1, indent_level - 1)
                if code.strip().endswith(':'):
                    code += f"{'    ' * indent_level}pass\n"
                    
            elif type(node).__name__ == "WhileNode":
                cond = node.condition_edit.text() or "True"
                code += f"{'    ' * indent_level}while {cond}:\n"
                indent_level += 1
                    
            elif type(node).__name__ == "DerivativeNode":
                idx = node.index_edit.value()
                val_var = get_input_var(node.inputs[0], indent_level)
                code += f"{'    ' * indent_level}dY[i, {idx}] = {val_var}\n"

            elif type(node).__name__ == "AssignmentNode":
                var_name = node.var_name_edit.text() or "var"
                usage = getattr(node, "value_mode", "scalar")
                val_var = get_input_var(node.inputs[0], indent_level, usage=usage)
                code += f"{'    ' * indent_level}{var_name} = {val_var}\n"
                for sock in getattr(node, 'outputs', []):
                    vars_map[sock] = var_name
                
            elif type(node).__name__ == "PythonCodeNode":
                generate_data_node(node, indent_level)
                
            elif type(node).__name__ == "ContinueNode":
                code += f"{'    ' * indent_level}continue\n"
                
            elif type(node).__name__ == "BreakNode":
                code += f"{'    ' * indent_level}break\n"

        if len(main_chain) == 0:
            return "def ode(Y, t, params, cell_positions):\n    pass\n"

        footer_code = getattr(begin_node, 'footer_code', None)
        if footer_code:
            code += footer_code
        else:
            code += """    Y_new = Y + dY * p.get('dT', 1/20.0)
    dY = (Y_new - Y) / p.get('dT', 1/20.0)
    return dY"""
        code += getattr(begin_node, "suffix_code", "")
        return code

    def serialize_scene(self):
        nodes_data = []
        edges_data = []
        current_fingerprint = self._scene_semantic_fingerprint()
        
        node_to_id = {}
        nodes_list = [item for item in self.scene.items() if hasattr(item, "serialize")]
        
        for i, node in enumerate(nodes_list):
            node_to_id[node] = i
            node_data = node.serialize()
            node_data["id"] = i
            if hasattr(node, "sequence_index"):
                node_data["sequence_index"] = getattr(node, "sequence_index")
            if type(node).__name__ == "OdeBeginNode":
                node_data["params_str"] = getattr(node, "params_str", "")
                node_data["color_func_str"] = getattr(node, "color_func_str", "")
                node_data["footer_code"] = getattr(node, "footer_code", "")
                node_data["prefix_code"] = getattr(node, "prefix_code", "")
                node_data["suffix_code"] = getattr(node, "suffix_code", "")
                if getattr(node, "source_fingerprint", "") == current_fingerprint:
                    node_data["source_ode_str"] = getattr(node, "source_ode_str", "")
                    node_data["source_fingerprint"] = getattr(node, "source_fingerprint", "")
                else:
                    node_data["source_ode_str"] = ""
                    node_data["source_fingerprint"] = ""
            nodes_data.append(node_data)
            
        from .node_core import NodeEdge
        edges_list = [item for item in self.scene.items() if isinstance(item, NodeEdge)]
        for edge in edges_list:
            if not edge.source_socket or not edge.dest_socket:
                continue
            src_node = edge.source_socket.node
            dst_node = edge.dest_socket.node
            if src_node not in node_to_id or dst_node not in node_to_id:
                continue
            
            if edge.edge_type == "Sequence":
                if hasattr(src_node, 'then_out') and edge.source_socket == src_node.then_out:
                    src_type = "then_out"
                elif hasattr(src_node, 'else_out') and edge.source_socket == src_node.else_out:
                    src_type = "else_out"
                else:
                    src_type = "seq_out" if edge.source_socket == getattr(src_node, "seq_out", None) else "seq_in"
                dst_type = "seq_in" if edge.dest_socket == getattr(dst_node, "seq_in", None) else "seq_out"
                src_idx = 0
                dst_idx = 0
            else:
                src_idx = src_node.outputs.index(edge.source_socket) if edge.source_socket in getattr(src_node, "outputs", []) else (src_node.inputs.index(edge.source_socket) if hasattr(src_node, "inputs") and edge.source_socket in src_node.inputs else 0)
                src_type = edge.source_socket.socket_type
                
                dst_idx = dst_node.inputs.index(edge.dest_socket) if edge.dest_socket in getattr(dst_node, "inputs", []) else (dst_node.outputs.index(edge.dest_socket) if hasattr(dst_node, "outputs") and edge.dest_socket in dst_node.outputs else 0)
                dst_type = edge.dest_socket.socket_type
            
            edges_data.append({
                "source_node": node_to_id[src_node],
                "source_idx": src_idx,
                "source_type": src_type,
                "dest_node": node_to_id[dst_node],
                "dest_idx": dst_idx,
                "dest_type": dst_type,
                "edge_type": edge.edge_type
            })
            
        data = {"nodes": nodes_data, "edges": edges_data}
        return json.dumps(data, indent=2)

    def deserialize_scene(self, json_data):
        try:
            data = json.loads(json_data)
        except Exception as e:
            print("Failed to parse JSON:", e)
            return

        self.scene.clear()
        
        from .node_types import (GeneNode, CellNode, MathNode, TextNode, ParameterNode, StateNode, 
                                 DerivativeNode, PythonCodeNode, OdeBeginNode, ForEachCellNode, 
                                 ForEachNeighborNode, EndOfForLoopNode, ConstantNode, ElseNode,
                                 IfNode, EndIfNode, VariableReadNode, NamedConstantNode, AggregateNode,
                                 ContinueNode, BreakNode, OdeEndNode, ColorMappingNode,
                                 WhileNode, EndWhileNode, AssignmentNode, ForLoopNode)
        from .node_core import NodeEdge
        
        node_map = {}
        type_mapping = {
            "GeneNode": GeneNode, "CellNode": CellNode, "MathNode": MathNode, "TextNode": TextNode, 
            "ParameterNode": ParameterNode, "StateNode": StateNode, "DerivativeNode": DerivativeNode, 
            "PythonCodeNode": PythonCodeNode, "OdeBeginNode": OdeBeginNode, "ForEachCellNode": ForEachCellNode,
            "ForEachNeighborNode": ForEachNeighborNode, "EndOfForLoopNode": EndOfForLoopNode,
            "ConstantNode": ConstantNode, "ElseNode": ElseNode, "IfNode": IfNode, "EndIfNode": EndIfNode,
            "VariableReadNode": VariableReadNode, "NamedConstantNode": NamedConstantNode,
            "AggregateNode": AggregateNode, "ContinueNode": ContinueNode, "BreakNode": BreakNode,
            "OdeEndNode": OdeEndNode, "ColorMappingNode": ColorMappingNode,
            "WhileNode": WhileNode, "EndWhileNode": EndWhileNode,
            "AssignmentNode": AssignmentNode, "ForLoopNode": ForLoopNode
        }
        
        for n_data in data.get("nodes", []):
            type_name = n_data.get("type")
            cls = type_mapping.get(type_name)
            if not cls: continue
            
            node = cls()
            node.deserialize(n_data)
            if "sequence_index" in n_data:
                node.sequence_index = n_data.get("sequence_index")
            if type_name == "OdeBeginNode":
                node.params_str = n_data.get("params_str", "")
                node.color_func_str = n_data.get("color_func_str", "")
                node.footer_code = n_data.get("footer_code", "")
                node.prefix_code = n_data.get("prefix_code", "")
                node.suffix_code = n_data.get("suffix_code", "")
                node.source_ode_str = n_data.get("source_ode_str", "")
                node.source_fingerprint = n_data.get("source_fingerprint", "")
            self.scene.addItem(node)
            node_map[n_data.get("id")] = node
            
        for e_data in data.get("edges", []):
            src_node = node_map.get(e_data.get("source_node"))
            dst_node = node_map.get(e_data.get("dest_node"))
            if not src_node or not dst_node: continue
            
            src_type = e_data.get("source_type", "output")
            dst_type = e_data.get("dest_type", "input")
            
            src_idx = e_data.get("source_idx", 0)
            dst_idx = e_data.get("dest_idx", 0)
            
            if src_type in ["seq_out", "seq_in", "then_out", "else_out"]:
                src_sockets = [getattr(src_node, src_type, None)]
            else:
                src_sockets = getattr(src_node, "outputs", []) if src_type == "output" else getattr(src_node, "inputs", [])
                
            if dst_type in ["seq_out", "seq_in", "then_out", "else_out"]:
                dst_sockets = [getattr(dst_node, dst_type, None)]
            else:
                dst_sockets = getattr(dst_node, "inputs", []) if dst_type == "input" else getattr(dst_node, "outputs", [])
                
            if 0 <= src_idx < len(src_sockets) and 0 <= dst_idx < len(dst_sockets):
                s1 = src_sockets[src_idx]
                s2 = dst_sockets[dst_idx]
                if s1 and s2:
                    edge = NodeEdge(s1, s2, e_data.get("edge_type", "Activation"))
                    s1.add_edge(edge)
                    s2.add_edge(edge)
                    self.scene.addItem(edge)
                    edge.update_path()
                
        needs_layout = True
        for node in self.scene.items():
            if hasattr(node, "inputs") and hasattr(node, "outputs"):
                if node.pos().x() != 0 or node.pos().y() != 0:
                    needs_layout = False
                    break
        if needs_layout:
            self.scene.arrange_nodes()

        try:
            current_fingerprint = self._scene_semantic_fingerprint()
            for node in self.scene.items():
                if type(node).__name__ == "OdeBeginNode" and getattr(node, "source_ode_str", ""):
                    node.source_fingerprint = current_fingerprint
        except Exception:
            pass
