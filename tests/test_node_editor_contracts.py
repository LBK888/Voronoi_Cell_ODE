import ast
from pathlib import Path
import unittest


PROJECT_ROOT = Path(__file__).resolve().parents[1]


class NodeEditorContractTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.main_source = (PROJECT_ROOT / "gui" / "main_window.py").read_text(
            encoding="utf-8"
        )
        cls.editor_source = (PROJECT_ROOT / "gui" / "node_editor.py").read_text(
            encoding="utf-8"
        )
        cls.main_tree = ast.parse(cls.main_source)
        cls.editor_tree = ast.parse(cls.editor_source)

    def test_main_window_exposes_visual_editor_button_menu_and_shortcut(self):
        self.assertIn('QPushButton("Visual Node Editor...")', self.main_source)
        self.assertIn("menubar.addMenu('Model')", self.main_source)
        self.assertIn("QKeySequence('Ctrl+Shift+N')", self.main_source)

    def test_editor_uses_a_document_boundary_instead_of_parent_widgets(self):
        document_class = next(
            node
            for node in self.editor_tree.body
            if isinstance(node, ast.ClassDef) and node.name == "NodeEditorDocument"
        )
        field_names = {
            target.id
            for statement in document_class.body
            if isinstance(statement, ast.AnnAssign)
            and isinstance(statement.target, ast.Name)
            for target in [statement.target]
        }
        self.assertEqual(
            field_names,
            {"ode_source", "params_source", "color_func_source"},
        )
        self.assertNotIn("self.parent().ode_edit", self.editor_source)
        self.assertNotIn("self.parent().params_edit", self.editor_source)
        self.assertNotIn("self.parent().color_func_edit", self.editor_source)
        self.assertIn("result = dlg.result_document()", self.main_source)

    def test_editor_has_one_primary_apply_action_and_explicit_cancel(self):
        self.assertIn('QPushButton("Apply to Model")', self.editor_source)
        self.assertIn('setObjectName("primaryButton")', self.editor_source)
        self.assertIn('QPushButton("Close")', self.editor_source)
        self.assertIn("self.close_btn.clicked.connect(self.reject)", self.editor_source)

    def test_editor_chrome_uses_new_light_semantic_palette(self):
        for color in ("#1A6FB5", "#0D8F6B", "#B45309", "#C4211C", "#6D3BBF"):
            self.assertIn(color, self.editor_source)
        self.assertIn('self.sidebar.setObjectName("sidePanel")', self.editor_source)
        self.assertIn('self.status_label.setObjectName("editorStatus")', self.editor_source)


if __name__ == "__main__":
    unittest.main()
