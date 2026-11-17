import importlib
import sys
print("Python:", sys.version)
for module in ("streamlit", "Bio", "rdkit", "prolif"):
    try:
        loaded = importlib.import_module(module)
        print(module, "OK", getattr(loaded, "__version__", ""))
    except Exception as exc:
        print(module, "НЕ ДОСТУПЕН:", type(exc).__name__, str(exc))
