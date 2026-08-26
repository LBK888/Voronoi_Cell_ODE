import ast
from pathlib import Path
import unittest


PROJECT_ROOT = Path(__file__).resolve().parents[1]


class ExportWorkerContractTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.main_source = (PROJECT_ROOT / "gui" / "main_window.py").read_text(
            encoding="utf-8"
        )
        cls.main_tree = ast.parse(cls.main_source)
        cls.worker_tree = ast.parse(
            (PROJECT_ROOT / "gui" / "results_export_worker.py").read_text(
                encoding="utf-8"
            )
        )

    def test_main_window_worker_arguments_match_constructor(self):
        worker_class = next(
            node
            for node in self.worker_tree.body
            if isinstance(node, ast.ClassDef) and node.name == "ResultsExportWorker"
        )
        constructor = next(
            node
            for node in worker_class.body
            if isinstance(node, ast.FunctionDef) and node.name == "__init__"
        )
        accepted = {
            argument.arg
            for argument in constructor.args.args + constructor.args.kwonlyargs
            if argument.arg != "self"
        }
        worker_call = next(
            node
            for node in ast.walk(self.main_tree)
            if isinstance(node, ast.Call)
            and isinstance(node.func, ast.Name)
            and node.func.id == "ResultsExportWorker"
        )
        passed = {keyword.arg for keyword in worker_call.keywords if keyword.arg}

        self.assertEqual(passed - accepted, set())
        self.assertIn("request", passed)

    def test_main_window_preserves_result_for_export_request(self):
        self.assertIn("self.simulation_result = result", self.main_source)
        self.assertIn("request = ExportRequest(", self.main_source)
        self.assertIn("simulation=self.simulation_result", self.main_source)


if __name__ == "__main__":
    unittest.main()
