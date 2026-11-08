"""Keep wheel scrolling from accidentally editing parameter controls."""

from PyQt6.QtCore import QObject, QEvent, QPointF
from PyQt6.QtGui import QWheelEvent
from PyQt6.QtWidgets import (
    QApplication, QAbstractScrollArea, QAbstractSlider, QAbstractSpinBox,
    QComboBox, QLineEdit, QScrollBar,
)


class ParameterWheelGuard(QObject):
    """Cover existing and dynamically created controls, including dialogs."""

    def eventFilter(self, watched, event):
        if event.type() != QEvent.Type.Wheel:
            return False
        control = watched
        if isinstance(control, QLineEdit):
            control = control.parentWidget()
        if not isinstance(control, (QAbstractSpinBox, QComboBox, QAbstractSlider)):
            return False
        if isinstance(control, QScrollBar):
            return False

        # Forward scrolling to the containing page, even when the input has
        # keyboard focus. Popup lists and plot canvases keep their own wheel
        # behavior because they are not parameter input controls.
        parent = control.parentWidget()
        while parent is not None:
            if isinstance(parent, QAbstractScrollArea):
                viewport = parent.viewport()
                forwarded = QWheelEvent(
                    QPointF(viewport.mapFromGlobal(event.globalPosition().toPoint())),
                    event.globalPosition(), event.pixelDelta(), event.angleDelta(),
                    event.buttons(), event.modifiers(), event.phase(),
                    event.inverted(),
                )
                QApplication.sendEvent(viewport, forwarded)
                break
            parent = parent.parentWidget()
        event.accept()
        return True


def disable_wheel_value_changes():
    """Install once per application; retain the filter for the app's lifetime."""
    app = QApplication.instance()
    if app is None:
        raise RuntimeError('Create QApplication before installing the wheel guard')
    if not hasattr(app, '_parameter_wheel_guard'):
        app._parameter_wheel_guard = ParameterWheelGuard(app)
        app.installEventFilter(app._parameter_wheel_guard)
