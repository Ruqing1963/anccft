# Convenience runner: stages 1-3 of the Paper V verification battery (see anccft5.tex, Table 1).
import subprocess, sys
for s in ("iw_stage1.py","iw_stage2.py","iw_stage3.py"):
    print(f"--- {s} ---"); subprocess.run([sys.executable,s],check=True)
