"""Vector-first figure export and clipboard actions shared by result canvases."""

from io import BytesIO
from pathlib import Path

import matplotlib as mpl
from PyQt6.QtCore import Qt
from PyQt6.QtGui import QImage
from PyQt6.QtWidgets import QApplication, QFileDialog, QMenu


def save_figure_dialog(parent, figure, default_name='frame.pdf'):
    path, selected = QFileDialog.getSaveFileName(
        parent, 'Save current image', default_name,
        'PDF Files (*.pdf);;SVG Files (*.svg);;PNG Files (*.png)',
    )
    if not path:
        return None
    if not Path(path).suffix:
        suffix = '.svg' if selected.startswith('SVG') else '.png' if selected.startswith('PNG') else '.pdf'
        path += suffix
    with mpl.rc_context({'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none'}):
        figure.savefig(path, dpi=180, bbox_inches='tight')
    return path


def copy_figure_to_clipboard(figure):
    """Copy the rendered figure (including legends), independent of window size."""
    buffer = BytesIO()
    figure.savefig(buffer, format='png', dpi=150, bbox_inches='tight')
    image = QImage.fromData(buffer.getvalue(), 'PNG')
    if image.isNull():
        raise RuntimeError('Unable to render the current figure for the clipboard')
    QApplication.clipboard().setImage(image)


def install_figure_menu(canvas, parent, default_name='frame.pdf', enabled=None):
    def show_menu(position):
        if enabled is not None and not enabled():
            return
        menu = QMenu(canvas)
        copy_action = menu.addAction('Copy current image')
        save_action = menu.addAction('Save image as...')
        action = menu.exec(canvas.mapToGlobal(position))
        try:
            if action == copy_action:
                copy_figure_to_clipboard(canvas.figure)
                if hasattr(parent, 'status_bar'):
                    parent.status_bar.showMessage('Current image copied to clipboard.', 4000)
            elif action == save_action:
                save_figure_dialog(parent, canvas.figure, default_name)
        except Exception as error:
            if hasattr(parent, 'status_bar'):
                parent.status_bar.showMessage(f'Image export failed: {error}')
            else:
                from PyQt6.QtWidgets import QMessageBox
                QMessageBox.warning(parent, 'Image export failed', str(error))

    canvas.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
    canvas.customContextMenuRequested.connect(show_menu)
