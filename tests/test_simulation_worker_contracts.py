import ast
from pathlib import Path
import unittest


PROJECT_ROOT = Path(__file__).resolve().parents[1]


class SimulationWorkerContractTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.main_tree = ast.parse(
            (PROJECT_ROOT / "gui" / "main_window.py").read_text(encoding="utf-8")
        )
        cls.worker_tree = ast.parse(
            (PROJECT_ROOT / "gui" / "simulation_worker.py").read_text(
                encoding="utf-8"
            )
        )

    def test_main_window_worker_arguments_match_constructor(self):
        worker_class = next(
            node
            for node in self.worker_tree.body
            if isinstance(node, ast.ClassDef) and node.name == "SimulationWorker"
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
            and node.func.id == "SimulationWorker"
        )
        passed = {keyword.arg for keyword in worker_call.keywords if keyword.arg}

        self.assertEqual(
            passed - accepted,
            set(),
            "MainWindow passes keywords that SimulationWorker does not accept",
        )
        self.assertIn("session", passed)

    def test_finished_slot_accepts_the_single_result_object(self):
        main_window = next(
            node
            for node in self.main_tree.body
            if isinstance(node, ast.ClassDef) and node.name == "MainWindow"
        )
        finished_slot = next(
            node
            for node in main_window.body
            if isinstance(node, ast.FunctionDef)
            and node.name == "on_simulation_finished"
        )
        parameters = [
            argument.arg
            for argument in finished_slot.args.args
            if argument.arg != "self"
        ]
        self.assertEqual(parameters, ["result"])

    def test_dynamic_ode_is_not_wrapped_in_lazy_numba_jit(self):
        session_tree = ast.parse(
            (PROJECT_ROOT / "gui" / "simulation_session.py").read_text(
                encoding="utf-8"
            )
        )
        session_class = next(
            node
            for node in session_tree.body
            if isinstance(node, ast.ClassDef) and node.name == "SimulationSession"
        )
        compile_model = next(
            node
            for node in session_class.body
            if isinstance(node, ast.FunctionDef) and node.name == "_compile_model"
        )
        jit_calls = [
            node
            for node in ast.walk(compile_model)
            if isinstance(node, ast.Call)
            and isinstance(node.func, ast.Name)
            and node.func.id == "jit"
        ]
        self.assertEqual(jit_calls, [])


if __name__ == "__main__":
    unittest.main()
