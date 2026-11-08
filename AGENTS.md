# Project workflow

- Maintain source, tests, model templates and documentation in this checkout. The user's runtime checkout is `F:\Voronoi_Cell_ODE`, with Python at `.venv\Scripts\python.exe`; use that checkout for execution and generated test artifacts. Do not create a virtual environment in the pCloud source checkout.
- Synchronize changed code/docs/tests to the runtime checkout as needed; preserve active `config.json`, user results and `.venv`.
- When `.codegraph/` exists, use CodeGraph before searching/reading code to understand symbols and their callers. Do not re-index unless requested.
- For converting paper equations into this program's model scripts, read [skills/paper-to-voronoi-ode/SKILL.md](skills/paper-to-voronoi-ode/SKILL.md). It covers `ode`, `params`, color mapping, mechanics, numerical validation and loadable configurations.
- Prefer the user's light/dark semantic color palettes for UI states. The existing application uses a light interface; information/actions use `#1A6FB5`, success `#0D8F6B`, warnings `#B45309`, errors `#C4211C`, and normal text `#374151`.
- Run the relevant behavior tests using the runtime interpreter. Full regression: `.\.venv\Scripts\python.exe -m unittest discover -s tests -v` from the runtime checkout. Qt UI tests set the offscreen platform and isolate configurations in temporary directories.
