import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from gui.app_paths import configure_ffmpeg, data_path, resource_path


class ApplicationPathTests(unittest.TestCase):
    def test_source_checkout_keeps_existing_configuration_location(self):
        with patch.dict(os.environ, {}, clear=True), patch('sys.frozen', False, create=True):
            self.assertEqual(data_path('config.json'), Path(__file__).resolve().parents[1] / 'config.json')

    def test_frozen_config_persists_outside_temporary_bundle(self):
        with tempfile.TemporaryDirectory() as folder:
            bundle = Path(folder, 'temporary-bundle')
            profile = Path(folder, 'local-app-data')
            with patch.dict(os.environ, {'LOCALAPPDATA': str(profile)}, clear=True), \
                    patch('sys.frozen', True, create=True), patch('sys._MEIPASS', str(bundle), create=True):
                self.assertEqual(resource_path('templates'), bundle / 'templates')
                config = data_path('config.json')
                config.write_text('{"saved": true}', encoding='utf-8')
                self.assertEqual(config, profile / 'VoronoiCellODE' / 'config.json')
            with patch.dict(os.environ, {'LOCALAPPDATA': str(profile)}, clear=True), \
                    patch('sys.frozen', True, create=True), patch('sys._MEIPASS', str(bundle / 'next-launch'), create=True):
                self.assertEqual(data_path('config.json').read_text(encoding='utf-8'), '{"saved": true}')

    def test_smoke_data_override_does_not_write_user_settings(self):
        with tempfile.TemporaryDirectory() as folder, \
                patch.dict(os.environ, {'VORONOI_CELL_ODE_DATA_DIR': folder}):
            self.assertEqual(data_path('app.log'), Path(folder, 'app.log'))

    def test_bundled_encoder_does_not_require_system_path(self):
        import matplotlib
        with tempfile.TemporaryDirectory() as folder, \
                patch('sys._MEIPASS', folder, create=True), matplotlib.rc_context():
            encoder = Path(folder, 'bin', 'ffmpeg.exe')
            encoder.parent.mkdir()
            encoder.write_bytes(b'fixture')
            configure_ffmpeg()
            self.assertEqual(matplotlib.rcParams['animation.ffmpeg_path'], str(encoder))


if __name__ == '__main__':
    unittest.main()
