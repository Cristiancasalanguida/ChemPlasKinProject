
import shutil
import subprocess
from pathlib import Path

air_in_values     = [116/3600, 180/3600, 255/3600]  # in kg/s
phi_global_values = [0.30, 0.40, 0.50, 0.60]

SCRIPT_DIR    = Path(__file__).parent
base          = Path("/home/cristiancasalanguida/Thesis/CRN_results/CRN_Mad")
folders       = ["air_116_T_700", "air_180_T_450", "air_255_T_310"]
phis          = ["phi030", "phi040", "phi050", "phi060"]

crndict_dst   = SCRIPT_DIR / "crnDict"
output_src    = SCRIPT_DIR / "output"
build_dir     = SCRIPT_DIR / "build"
params_path   = SCRIPT_DIR / "massflow_params.txt"

for folder, air_in in zip(folders, air_in_values):
    for phi, phi_global in zip(phis, phi_global_values):
        dir_phi = base / folder / phi
        print(f"\n=== {folder} / {phi} ===")

        params_path.write_text(f"{air_in:.4f}\n{phi_global:.4f}\n")
        shutil.copy(dir_phi / "crnDict_MAD", crndict_dst)

        subprocess.run(["./ChemPlasKin"], cwd=build_dir, check=True)

        output_dst = dir_phi / "output"
        if output_dst.exists():
            shutil.rmtree(output_dst)
        shutil.copytree(output_src, output_dst)

        crndict_dst.unlink()
        print(f"Done: results saved to {output_dst}")
