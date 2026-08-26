import sys
import os
import json

sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))

from PyQt6.QtWidgets import QApplication
from PyQt6.QtWidgets import QGraphicsScene
from gui.node_parser import PythonToNodeParser, NodeToPythonGenerator

app = QApplication(sys.argv)
scene = QGraphicsScene()
scene.arrange_nodes = lambda: None

with open("config_psm.json", "r", encoding="utf-8") as f:
    config = json.load(f)

ode_str = config["ode"]
params_str = config["params"]

parser = PythonToNodeParser(scene)
parser.parse_ode_string(ode_str, params_str)

generator = NodeToPythonGenerator(scene)
generated_code = generator.generate_ode().strip()

with open("gen.py", "w", encoding="utf-8") as f:
    f.write(generated_code)

print("Wrote to gen.py!")
