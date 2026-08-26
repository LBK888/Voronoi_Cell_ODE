import ast
from pathlib import Path
import unittest


PROJECT_ROOT = Path(__file__).resolve().parents[1]


class SimulationUtilityContractTests(unittest.TestCase):
    def test_repulsion_dependencies_are_imported(self):
        tree = ast.parse(
            (PROJECT_ROOT / "gui" / "sim_utils.py").read_text(encoding="utf-8")
        )
        repulsion = next(
            node
            for node in tree.body
            if isinstance(node, ast.FunctionDef)
            and node.name == "repulsion_move_neighbors_no_outer"
        )
        loaded_names = {
            node.id
            for node in ast.walk(repulsion)
            if isinstance(node, ast.Name) and isinstance(node.ctx, ast.Load)
        }
        imported_names = {
            alias.asname or alias.name
            for node in tree.body
            if isinstance(node, ast.ImportFrom)
            for alias in node.names
        }
        imported_names.update(
            (alias.asname or alias.name).split(".")[0]
            for node in tree.body
            if isinstance(node, ast.Import)
            for alias in node.names
        )

        required_globals = {"Voronoi", "defaultdict"} & loaded_names
        self.assertEqual(
            required_globals - imported_names,
            set(),
            "Repulsion helper references dependencies that are not imported",
        )


if __name__ == "__main__":
    unittest.main()
