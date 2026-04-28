import ast
import json

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
                                 OdeBeginNode, ForEachCellNode, ForEachNeighborNode, EndOfForLoopNode,
                                 ConstantNode, ElseNode, IfNode, EndIfNode, VariableReadNode,
                                 NamedConstantNode, AggregateNode, ContinueNode, BreakNode, OdeEndNode,
                                 WhileNode, EndWhileNode)
        from .node_core import NodeEdge
        
        try:
            tree = ast.parse(ode_str)
        except Exception as e:
            print("Failed to AST parse ODE:", e)
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
                footer_lines.append(code)
            else:
                core_body.append(stmt)
            
        env = {}
        x_offset = 100
        y_offset = 100
        
        begin_node = OdeBeginNode()
        begin_node.setPos(x_offset, y_offset)
        begin_node.params_str = params_str
        begin_node.color_func_str = color_func_str
        if footer_lines:
            begin_node.footer_code = "\n".join(["    " + line for line in footer_lines]) + "\n"
        self.scene.addItem(begin_node)
        
        y_offset += 150
        
        def link_seq(n1, n2, out_attr='seq_out', in_attr='seq_in'):
            if n1 and n2:
                out_sock = getattr(n1, out_attr, None)
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

        def get_or_create_param(name, x, y):
            if name not in self.param_registry:
                p_node = ParameterNode()
                p_node.param_name_edit.setText(name)
                p_node.setPos(x, y)
                self.scene.addItem(p_node)
                self.param_registry[name] = p_node
            return self.param_registry[name].outputs[0]

        def extract_column_index(subscript):
            if hasattr(subscript, 'slice'):
                sl = subscript.slice
                if isinstance(sl, ast.Tuple) and len(sl.elts) == 2:
                    idx_node = sl.elts[1]
                    if isinstance(idx_node, ast.Constant):
                        return idx_node.value
            return None

        # Pre-scan state aliases
        for stmt in core_body:
            if isinstance(stmt, ast.Assign) and len(stmt.targets) == 1 and isinstance(stmt.targets[0], ast.Name):
                if isinstance(stmt.value, ast.Subscript) and isinstance(stmt.value.value, ast.Name) and stmt.value.value.id == "Y":
                    var_name = stmt.targets[0].id
                    idx = extract_column_index(stmt.value)
                    if idx is not None:
                        s_node = StateNode()
                        s_node.index_edit.setValue(idx)
                        s_node.var_name = var_name
                        s_node.setPos(x_offset + 300, y_offset)
                        self.scene.addItem(s_node)
                        env[var_name] = s_node.outputs[0]

        def build_expr(expr, x_hints, y_hints):
            pending_nodes = []
            
            if isinstance(expr, ast.Name):
                if expr.id in env:
                    return env[expr.id], pending_nodes
                    
            elif isinstance(expr, ast.Constant):
                c_node = ConstantNode()
                c_node.value_edit.setText(str(expr.value))
                c_node.setPos(x_hints, y_hints)
                pending_nodes.append(c_node)
                return c_node.outputs[0], pending_nodes
                
            elif isinstance(expr, ast.UnaryOp) and isinstance(expr.op, ast.USub) and isinstance(expr.operand, ast.Constant):
                c_node = ConstantNode()
                c_node.value_edit.setText(str(-float(expr.operand.value)))
                c_node.setPos(x_hints, y_hints)
                pending_nodes.append(c_node)
                return c_node.outputs[0], pending_nodes
                
            elif isinstance(expr, ast.Call):
                if isinstance(expr.func, ast.Attribute) and getattr(expr.func.value, 'id', '') in ["params", "p"] and expr.func.attr == "get":
                    if len(expr.args) >= 1 and isinstance(expr.args[0], ast.Constant):
                        target_name = str(expr.args[0].value)
                        sock = get_or_create_param(target_name, x_hints, y_hints)
                        sock.node.access_style = "get"
                        return sock, pending_nodes
                elif isinstance(expr.func, ast.Attribute) and expr.func.attr in ['sin', 'cos', 'tan', 'exp', 'log', 'sqrt', 'abs']:
                    if len(expr.args) >= 1:
                        in_socket, sub_pending = build_expr(expr.args[0], x_hints - 150, y_hints)
                        pending_nodes.extend(sub_pending)
                        if in_socket:
                            m_node = MathNode()
                            op_idx = {'sin': 5, 'cos': 6, 'tan': 7, 'exp': 8, 'log': 9, 'sqrt': 10, 'abs': 11}[expr.func.attr]
                            m_node.op_combo.setCurrentIndex(op_idx)
                            m_node.on_op_changed(op_idx)
                            m_node.setPos(x_hints, y_hints)
                            pending_nodes.append(m_node)
                            
                            edge = NodeEdge(in_socket, m_node.inputs[0], "Data")
                            in_socket.add_edge(edge)
                            m_node.inputs[0].add_edge(edge)
                            pending_nodes.append(edge)
                            return m_node.outputs[0], pending_nodes

                c_node = PythonCodeNode()
                c_node.code_edit.setPlainText(ast.unparse(expr).strip())
                c_node.setPos(x_hints, y_hints)
                c_node.is_fallback = True
                pending_nodes.append(c_node)
                return c_node.outputs[0], pending_nodes
                
            elif isinstance(expr, ast.Subscript):
                slice_node = expr.slice
                if hasattr(slice_node, 'value') and not isinstance(slice_node, ast.Constant):
                    slice_node = slice_node.value
                
                if isinstance(expr.value, ast.Name) and expr.value.id in env:
                    idx_name = getattr(slice_node, 'id', '')
                    if getattr(slice_node, 'value', None) and isinstance(getattr(slice_node, 'value'), ast.Name):
                        idx_name = getattr(slice_node, 'value').id
                    if isinstance(slice_node, ast.Name):
                        idx_name = slice_node.id
                    if idx_name == "i":
                        sock = env[expr.value.id]
                        if hasattr(sock, "node") and type(sock.node).__name__ == "StateNode":
                            return sock, pending_nodes
                        
                if isinstance(expr.value, ast.Name) and expr.value.id in ["params", "p"]:
                    if isinstance(slice_node, ast.Constant):
                        target_name = str(slice_node.value)
                        sock = get_or_create_param(target_name, x_hints, y_hints)
                        sock.node.access_style = "subscript"
                        return sock, pending_nodes
                elif isinstance(expr.value, ast.Name) and getattr(expr.value, 'id', '') == "Y":
                    if getattr(slice_node, 'elts', None) and len(slice_node.elts) == 2:
                        idx_node = slice_node.elts[1]
                        if hasattr(idx_node, 'value') and not isinstance(idx_node, ast.Constant):
                            idx_node = idx_node.value
                        if isinstance(idx_node, ast.Constant):
                            s_node = StateNode()
                            s_node.index_edit.setValue(idx_node.value)
                            s_node.setPos(x_hints, y_hints)
                            pending_nodes.append(s_node)
                            return s_node.outputs[0], pending_nodes
                            
            elif isinstance(expr, ast.BinOp):
                left_socket, left_pending = build_expr(expr.left, x_hints - 200, y_hints - 100)
                right_socket, right_pending = build_expr(expr.right, x_hints - 200, y_hints + 100)
                pending_nodes.extend(left_pending)
                pending_nodes.extend(right_pending)
                
                if not left_socket or not right_socket:
                    return None, pending_nodes
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
                
                if left_socket:
                    edge = NodeEdge(left_socket, m_node.inputs[0], "Data")
                    left_socket.add_edge(edge)
                    m_node.inputs[0].add_edge(edge)
                    pending_nodes.append(edge)
                if right_socket:
                    edge = NodeEdge(right_socket, m_node.inputs[1], "Data")
                    right_socket.add_edge(edge)
                    m_node.inputs[1].add_edge(edge)
                    pending_nodes.append(edge)
                return m_node.outputs[0], pending_nodes
                
            return None, pending_nodes

        def process_stmts(stmts, current_prev):
            nonlocal y_offset
            unhandled_stmts = []

            def flush_unhandled(cp):
                nonlocal y_offset
                if not unhandled_stmts:
                    return cp
                lines = [ast.unparse(s).strip() for s in unhandled_stmts]
                code = "\n".join(lines)
                node = PythonCodeNode()
                node.code_edit.setPlainText(code)
                node.setPos(x_offset, y_offset)
                self.scene.addItem(node)
                link_seq(cp, node)
                y_offset += 150
                for stmt in unhandled_stmts:
                    for sub_node in ast.walk(stmt):
                        if isinstance(sub_node, ast.Assign) and len(sub_node.targets) == 1 and isinstance(sub_node.targets[0], ast.Name):
                            var_name = sub_node.targets[0].id
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
                            already_wired = False
                            if hasattr(node, "inputs"):
                                for s in node.inputs:
                                    if s.name == var_name:
                                        for edge in s.edges:
                                            if getattr(edge, 'source_socket', None) == sock or getattr(edge, 'dest_socket', None) == sock:
                                                already_wired = True
                                                break
                            if not already_wired:
                                in_sock = None
                                if hasattr(node, "inputs"):
                                    for s in node.inputs:
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
                code = ast.unparse(stmt).strip()
                if code == "continue":
                    current_prev = flush_unhandled(current_prev)
                    node = ContinueNode()
                    node.setPos(x_offset, y_offset)
                    self.scene.addItem(node)
                    link_seq(current_prev, node)
                    current_prev = node
                    y_offset += 150
                    continue
                elif code == "break":
                    current_prev = flush_unhandled(current_prev)
                    node = BreakNode()
                    node.setPos(x_offset, y_offset)
                    self.scene.addItem(node)
                    link_seq(current_prev, node)
                    current_prev = node
                    y_offset += 150
                    continue
                    
                if code == "neibs = neighbors[i]":
                    continue

                handled = False

                if isinstance(stmt, ast.For) and isinstance(stmt.target, ast.Name) and stmt.target.id == 'i':
                    current_prev = flush_unhandled(current_prev)
                    node = ForEachCellNode()
                    node.setPos(x_offset, y_offset)
                    self.scene.addItem(node)
                    link_seq(current_prev, node)
                    current_prev = node
                    y_offset += 150
                    
                    current_prev = process_stmts(stmt.body, current_prev)
                    
                    end_node = EndOfForLoopNode()
                    end_node.setPos(x_offset, y_offset)
                    self.scene.addItem(end_node)
                    link_seq(current_prev, end_node)
                    current_prev = end_node
                    y_offset += 150
                    handled = True

                elif isinstance(stmt, ast.If):
                    current_prev = flush_unhandled(current_prev)
                    cond_str = ast.unparse(stmt.test).strip()
                    
                    if cond_str == "neibs":
                        node = ForEachNeighborNode()
                        node.setPos(x_offset, y_offset)
                        self.scene.addItem(node)
                        link_seq(current_prev, node)
                        current_prev = node
                        y_offset += 150
                        
                        current_prev = process_stmts(stmt.body, current_prev)
                        
                        end_node = EndOfForLoopNode()
                        end_node.setPos(x_offset, y_offset)
                        self.scene.addItem(end_node)
                        link_seq(current_prev, end_node)
                        current_prev = end_node
                        y_offset += 150
                        handled = True
                    else:
                        node = IfNode()
                        node.condition_edit.setText(cond_str)
                        node.setPos(x_offset, y_offset)
                        self.scene.addItem(node)
                        link_seq(current_prev, node)
                        y_offset += 150
                        
                        then_prev = process_stmts(stmt.body, node)
                        
                        else_prev = None
                        if stmt.orelse:
                            else_prev = process_stmts(stmt.orelse, node)
                        
                        end_if = EndIfNode()
                        end_if.setPos(x_offset, y_offset)
                        self.scene.addItem(end_if)
                        
                        link_seq(then_prev, end_if, out_attr='seq_out' if then_prev != node else 'then_out', in_attr='seq_in')
                        if else_prev:
                            link_seq(else_prev, end_if, out_attr='seq_out' if else_prev != node else 'else_out', in_attr='seq_in')
                        elif node.else_out:
                            link_seq(node, end_if, out_attr='else_out', in_attr='seq_in')

                        current_prev = end_if
                        y_offset += 150
                        handled = True

                elif isinstance(stmt, ast.While):
                    current_prev = flush_unhandled(current_prev)
                    cond_str = ast.unparse(stmt.test).strip()
                    
                    node = WhileNode()
                    node.condition_edit.setText(cond_str)
                    node.setPos(x_offset, y_offset)
                    self.scene.addItem(node)
                    link_seq(current_prev, node)
                    y_offset += 150
                    
                    body_prev = process_stmts(stmt.body, node)
                    
                    end_while = EndWhileNode()
                    end_while.setPos(x_offset, y_offset)
                    self.scene.addItem(end_while)
                    link_seq(body_prev, end_while)
                    
                    current_prev = end_while
                    y_offset += 150
                    handled = True

                elif isinstance(stmt, ast.Assign):
                    if len(stmt.targets) == 1 and isinstance(stmt.targets[0], ast.Subscript):
                        tgt = stmt.targets[0]
                        if isinstance(tgt.value, ast.Name) and tgt.value.id == "dY":
                            slice_node = tgt.slice
                            if hasattr(slice_node, 'value') and not isinstance(slice_node, ast.Constant):
                                slice_node = slice_node.value
                            if getattr(slice_node, 'elts', None) and len(slice_node.elts) == 2:
                                idx_node = slice_node.elts[1]
                                if isinstance(idx_node, ast.Constant):
                                    current_prev = flush_unhandled(current_prev)
                                    out_socket, pending = build_expr(stmt.value, x_offset + 500, y_offset)
                                    if out_socket:
                                        commit_nodes(pending)
                                        node = DerivativeNode()
                                        node.index_edit.setValue(idx_node.value)
                                        node.setPos(x_offset + 200, y_offset)
                                        self.scene.addItem(node)
                                        link_seq(current_prev, node)
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
                            success_unpack = True
                            temp_env = {}
                            all_pending = []
                            for t_elt, v_elt in zip(tgt_tuple.elts, val_tuple.elts):
                                if isinstance(t_elt, ast.Name):
                                    out_sock, pending = build_expr(v_elt, x_offset + 200, y_offset)
                                    if out_sock and hasattr(out_sock.node, 'is_fallback') and getattr(out_sock.node, 'is_fallback'):
                                        success_unpack = False
                                        break
                                    if out_sock:
                                        temp_env[t_elt.id] = out_sock
                                        all_pending.extend(pending)
                                        y_offset += 150
                                    else:
                                        success_unpack = False
                                        break
                                else:
                                    success_unpack = False
                                    break
                            if success_unpack:
                                commit_nodes(all_pending)
                                env.update(temp_env)
                                handled = False # we don't break sequence for isolated assigns
                                for t_elt in tgt_tuple.elts:
                                    if isinstance(t_elt, ast.Name):
                                        sock = temp_env[t_elt.id]
                                        if sock and hasattr(sock, 'node') and type(sock.node).__name__ == "StateNode":
                                            sock.node.var_name = t_elt.id
                            else:
                                pass

                    elif len(stmt.targets) == 1 and isinstance(stmt.targets[0], ast.Name):
                        var_name = stmt.targets[0].id
                        out_socket, pending = build_expr(stmt.value, x_offset + 200, y_offset)
                        
                        if out_socket and not getattr(out_socket.node, 'is_fallback', False):
                            commit_nodes(pending)
                            env[var_name] = out_socket
                            y_offset += 150
                            
                            if type(out_socket.node).__name__ == "StateNode":
                                out_socket.node.var_name = var_name
                            handled = True # pure data flow
                        else:
                            handled = False

                if not handled:
                    unhandled_stmts.append(stmt)
                    
            current_prev = flush_unhandled(current_prev)
            return current_prev

        last_node = process_stmts(core_body, begin_node)
        
        end_node = OdeEndNode()
        end_node.setPos(x_offset, y_offset + 100)
        self.scene.addItem(end_node)
        link_seq(last_node, end_node)
        
        for item in self.scene.items():
            if hasattr(item, "code_edit"):
                item.code_edit.document().adjustSize()
            if hasattr(item, "update_sockets"):
                item.update_sockets()
                
        self.scene.arrange_nodes()
        if self.scene.views():
            from PyQt6.QtCore import Qt
            self.scene.views()[0].fitInView(self.scene.itemsBoundingRect(), Qt.AspectRatioMode.KeepAspectRatio)


class NodeToPythonGenerator:
    def __init__(self, scene):
        self.scene = scene
        self.warnings = []

    def collect_main_chain(self):
        chain = []
        visited = set()
        current_node = None
        for item in self.scene.items():
            if type(item).__name__ == "OdeBeginNode":
                current_node = item
                break
        
        def traverse(node):
            if not node or node in visited:
                return
            if getattr(node, 'bypassed', False):
                return
            visited.add(node)
            chain.append(node)
            
            if type(node).__name__ == "OdeEndNode":
                return
                
            if type(node).__name__ == "IfNode":
                for edge in getattr(node, 'then_out', []).edges if hasattr(node, 'then_out') and node.then_out else []:
                    out_sock = edge.dest_socket if edge.dest_socket else edge.source_socket
                    if out_sock and out_sock.node != node: traverse(out_sock.node)
                for edge in getattr(node, 'else_out', []).edges if hasattr(node, 'else_out') and node.else_out else []:
                    out_sock = edge.dest_socket if edge.dest_socket else edge.source_socket
                    if out_sock and out_sock.node != node: traverse(out_sock.node)
            else:
                if hasattr(node, 'seq_out') and node.seq_out:
                    for edge in node.seq_out.edges:
                        out_sock = edge.dest_socket if edge.dest_socket else edge.source_socket
                        if out_sock and out_sock.node != node: traverse(out_sock.node)
        traverse(current_node)
        return chain

    def generate_ode(self):
        self.warnings.clear()
        main_chain = self.collect_main_chain()
        
        code = ""
        indent_level = 1
        vars_map = {}
        var_counter = 0
        generating_set = set()

        def get_input_var(socket, current_indent):
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
                    if out_sock in vars_map:
                        return vars_map[out_sock]
                    v_name = generate_data_node(node, current_indent)
                    if v_name is not None:
                        return v_name
            return "0.0"

        def generate_data_node(node, current_indent):
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
                if var_name:
                    result_expr = f"{var_name}[i]"
                else:
                    result_expr = f"Y[i, {idx}]"
                
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
            
            is_complex = node_type in ["MathNode", "AggregateNode"]
            if is_complex and result_expr:
                var_name = f"_tmp_{var_counter}"
                var_counter += 1
                code += f"{indent}{var_name} = {result_expr}\n"
                for sock in getattr(node, 'outputs', []):
                    vars_map[sock] = var_name
                return var_name
                
            for sock in getattr(node, 'outputs', []):
                vars_map[sock] = result_expr
                
            return result_expr

        begin_node = None
        for node in main_chain:
            if type(node).__name__ == "OdeBeginNode":
                begin_node = node
                code += "def ode(Y, t, params, cell_positions):\n"
                code += "    p = params\n"
                code += "    dY = np.zeros_like(Y)\n"
                
            elif type(node).__name__ == "ForEachCellNode":
                code += f"{'    ' * indent_level}for i in range(Y.shape[0]):\n"
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
        return code

    def serialize_scene(self):
        nodes_data = []
        edges_data = []
        
        node_to_id = {}
        nodes_list = [item for item in self.scene.items() if hasattr(item, "serialize")]
        
        for i, node in enumerate(nodes_list):
            node_to_id[node] = i
            node_data = node.serialize()
            node_data["id"] = i
            if type(node).__name__ == "OdeBeginNode":
                node_data["params_str"] = getattr(node, "params_str", "")
                node_data["color_func_str"] = getattr(node, "color_func_str", "")
                node_data["footer_code"] = getattr(node, "footer_code", "")
            nodes_data.append(node_data)
            
        from .node_core import NodeEdge
        edges_list = [item for item in self.scene.items() if isinstance(item, NodeEdge)]
        for edge in edges_list:
            if not edge.source_socket or not edge.dest_socket:
                continue
            src_node = edge.source_socket.node
            dst_node = edge.dest_socket.node
            
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
                                 WhileNode, EndWhileNode)
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
            "WhileNode": WhileNode, "EndWhileNode": EndWhileNode
        }
        
        for n_data in data.get("nodes", []):
            type_name = n_data.get("type")
            cls = type_mapping.get(type_name)
            if not cls: continue
            
            node = cls()
            node.deserialize(n_data)
            if type_name == "OdeBeginNode":
                node.params_str = n_data.get("params_str", "")
                node.color_func_str = n_data.get("color_func_str", "")
                node.footer_code = n_data.get("footer_code", "")
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
