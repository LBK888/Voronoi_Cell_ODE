import ast
from pathlib import Path
import unittest


PROJECT_ROOT = Path(__file__).resolve().parents[1]


class GuiContractTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.main_source = (PROJECT_ROOT / "gui" / "main_window.py").read_text(
            encoding="utf-8"
        )
        cls.main_tree = ast.parse(cls.main_source)

    def test_main_window_uses_grid_based_preview_api(self):
        preview_tree = ast.parse(
            (PROJECT_ROOT / "gui" / "preview_canvas.py").read_text(encoding="utf-8")
        )
        preview_class = next(
            node
            for node in preview_tree.body
            if isinstance(node, ast.ClassDef) and node.name == "PreviewCanvas"
        )
        preview_methods = {
            node.name for node in preview_class.body if isinstance(node, ast.FunctionDef)
        }
        self.assertTrue({"set_grid", "get_grid"}.issubset(preview_methods))

        legacy_calls = []
        for node in ast.walk(self.main_tree):
            if not isinstance(node, ast.Call) or not isinstance(node.func, ast.Attribute):
                continue
            receiver = node.func.value
            if (
                isinstance(receiver, ast.Attribute)
                and isinstance(receiver.value, ast.Name)
                and receiver.value.id == "self"
                and receiver.attr == "preview_canvas"
                and node.func.attr in {"set_cells", "get_cells"}
            ):
                legacy_calls.append((node.func.attr, node.lineno))

        self.assertEqual(
            legacy_calls,
            [],
            f"MainWindow still calls the removed cell-based preview API: {legacy_calls}",
        )
        self.assertIn("self.preview_canvas.set_grid(", self.main_source)
        self.assertIn("self.preview_canvas.get_grid(", self.main_source)

    def test_grid_tab_can_reflow_in_a_narrow_settings_pane(self):
        self.assertIn(
            "grid_param_form.setRowWrapPolicy(QFormLayout.RowWrapPolicy.WrapLongRows)",
            self.main_source,
        )
        self.assertIn(
            "grid_param_form.setFieldGrowthPolicy(",
            self.main_source,
        )
        self.assertIn("grid_actions = QVBoxLayout()", self.main_source)

    def test_settings_signals_do_not_pass_values_as_config_paths(self):
        self.assertNotIn(
            ".connect(self.save_config)",
            self.main_source,
            "Qt signal values must not be forwarded to save_config(save_path)",
        )


if __name__ == "__main__":
    unittest.main()
