import ast
from PyQt6.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QPushButton, QMessageBox, QTextEdit
)
from PyQt6.QtGui import QFont, QColor
from PyQt6.QtCore import Qt, QTimer
from PyQt6.Qsci import QsciScintilla, QsciLexerPython

class CodeEditorDialog(QDialog):
    def __init__(self, parent=None, target_edit: QTextEdit = None, editor_type: str = 'ode'):
        super().__init__(parent)
        self.target_edit = target_edit
        self.editor_type = editor_type
        self.setWindowTitle(f"Code Editor - {editor_type.upper()}")
        self.resize(800, 600)
        
        # Init UI
        layout = QVBoxLayout(self)
        
        # Buttons
        btn_layout = QHBoxLayout()
        self.btn_import = QPushButton("Import Template")
        self.btn_cancel = QPushButton("Cancel and Undo Changes")
        self.btn_cancel.setStyleSheet("color: #800000; font-weight: bold;")
        self.btn_ok = QPushButton("OK & Close Editor")
        self.btn_ok.setStyleSheet("color: #000080; font-weight: bold;")

        self.btn_import.clicked.connect(self.on_import_template)
        self.btn_cancel.clicked.connect(self.on_cancel)
        self.btn_ok.clicked.connect(self.on_ok)
        
        btn_layout.addWidget(self.btn_import)
        btn_layout.addWidget(self.btn_cancel)
        btn_layout.addWidget(self.btn_ok)
        
        layout.addLayout(btn_layout)
        
        # Editor
        self.editor = QsciScintilla()
        font = QFont("Consolas", 11)
        self.editor.setFont(font)
        
        # Set lexer
        self.lexer = QsciLexerPython(self.editor)
        self.lexer.setDefaultFont(font)
        
        # Dark Mode Colors
        self.lexer.setDefaultPaper(QColor("#1e1e1e"))
        self.lexer.setDefaultColor(QColor("#dcdcdc"))
        for style in range(128):
            self.lexer.setPaper(QColor("#1e1e1e"), style)
            
        try:
            self.lexer.setColor(QColor("#569cd6"), QsciLexerPython.Keyword)
            self.lexer.setColor(QColor("#4ec9b0"), getattr(QsciLexerPython, 'ClassName', 8))
            self.lexer.setColor(QColor("#dcdcaa"), getattr(QsciLexerPython, 'FunctionMethodName', 9))
            self.lexer.setColor(QColor("#ce9178"), QsciLexerPython.DoubleQuotedString)
            self.lexer.setColor(QColor("#ce9178"), QsciLexerPython.SingleQuotedString)
            self.lexer.setColor(QColor("#6a9955"), QsciLexerPython.Comment)
            self.lexer.setColor(QColor("#6a9955"), getattr(QsciLexerPython, 'CommentBlock', 12))
            self.lexer.setColor(QColor("#b5cea8"), QsciLexerPython.Number)
            self.lexer.setColor(QColor("#dcdcdc"), QsciLexerPython.Operator)
            self.lexer.setColor(QColor("#dcdcdc"), QsciLexerPython.Identifier)
            self.lexer.setColor(QColor("#9cdcfe"), getattr(QsciLexerPython, 'UnclosedString', 7))
            self.lexer.setColor(QColor("#c586c0"), getattr(QsciLexerPython, 'Decorator', 15))
        except AttributeError:
            pass # Ignore if any specific style isn't found
            
        self.editor.setLexer(self.lexer)
        self.editor.setPaper(QColor("#1e1e1e"))
        
        # Line numbers
        self.editor.setMarginType(0, QsciScintilla.MarginType.NumberMargin)
        self.editor.setMarginWidth(0, "0000")
        self.editor.setMarginsBackgroundColor(QColor("#2d2d2d"))
        self.editor.setMarginsForegroundColor(QColor("#cccccc"))
        
        # Caret and selection color
        self.editor.setCaretForegroundColor(QColor("#ffffff"))
        self.editor.setSelectionBackgroundColor(QColor("#404040"))
        
        # Auto-completion / braces
        self.editor.setBraceMatching(QsciScintilla.BraceMatch.SloppyBraceMatch)
        
        # Tab size
        self.editor.setTabWidth(4)
        self.editor.setIndentationsUseTabs(False)
        self.editor.setAutoIndent(True)
        
        # Colors for annotations (Syntax errors)
        self.editor.setAnnotationDisplay(QsciScintilla.AnnotationDisplay.AnnotationStandard)
        self.editor.SendScintilla(QsciScintilla.SCI_STYLESETBACK, 254, QColor("#4a1c1c"))
        
        layout.addWidget(self.editor)
        
        if self.target_edit:
            self.editor.setText(self.target_edit.toPlainText())
            
        self.timer = QTimer()
        self.timer.setSingleShot(True)
        self.timer.timeout.connect(self.check_syntax)
        self.editor.textChanged.connect(lambda: self.timer.start(500))

    def check_syntax(self):
        self.editor.clearAnnotations()
        code = self.editor.text().strip()
        if not code:
            return
            
        try:
            if self.editor_type == 'params':
                if code.startswith('params ='):
                    code_to_check = code[len('params ='):].strip()
                else:
                    code_to_check = code
                if not code_to_check.startswith('dict'):
                    ast.parse(code_to_check)
            else:
                ast.parse(code)
        except SyntaxError as e:
            line = e.lineno - 1 if e.lineno else 0
            # Scintilla annotations
            self.editor.annotate(line, f"Syntax Error: {e.msg}", 0)
        except Exception as e:
            pass

    def on_import_template(self):
        reply = QMessageBox.question(
            self, 'Import Template', 
            'Are you sure you want to import the template?\nThis will overwrite current text.',
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No, 
            QMessageBox.StandardButton.No
        )
        if reply == QMessageBox.StandardButton.Yes:
            templates = {
                'ode': "def ode(Y, t, params, cell_positions):\n    p = params\n    dY = np.zeros_like(Y)\n    # Your equations here\n    return dY",
                'params': "params = {\n    'dT': 0.02,\n    'sigma': 0.01,\n    'n_var': 1\n}",
                'color_func': "def color_func(Y, mode='polygon', global_bounds=None):\n    colors = np.zeros((len(Y), 3))\n    # Example: colors[:, 0] = Y[:, 0] / np.max(Y[:, 0]) if np.max(Y[:, 0]) > 0 else 0\n    return colors"
            }
            self.editor.setText(templates.get(self.editor_type, ""))

    def on_ok(self):
        reply = QMessageBox.question(
            self, 'Save Changes', 
            'Are you sure you want to save changes to the main window?',
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No, 
            QMessageBox.StandardButton.Yes
        )
        if reply == QMessageBox.StandardButton.Yes:
            if self.target_edit:
                self.target_edit.setPlainText(self.editor.text())
            self.accept()

    def on_cancel(self):
        reply = QMessageBox.question(
            self, 'Cancel Changes', 
            'Are you sure you want to cancel and discard any changes?',
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No, 
            QMessageBox.StandardButton.No
        )
        if reply == QMessageBox.StandardButton.Yes:
            self.reject()
