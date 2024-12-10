import shutil
from pathlib import Path
import molpy as mp
from h_submitor import cmdline


@cmdline
def polymerize(typify_dir:dict, polymerize_config: dict, proj_dir:Path, backend: str) -> Path:

    name = polymerize_config['name']

    tleap_script = ["source leaprc.gaff", "source leaprc.water.tip3p"]


    if "seq" in polymerize_config:
        seq = polymerize_config['seq']
        for mon in set(seq):
            if mon in typify_dir:
                prepi = typify_dir[mon] / f"{mon}.prepi"
                frcmod = typify_dir[mon] / f"{mon}.frcmod"
                tleap_script.append(f"loadamberprep {prepi.absolute()}")
                tleap_script.append(f"loadamberparams {frcmod.absolute()}")

        tleap_script.append(f"chain = sequence {{ {' '.join(seq)} }}")
        tleap_script.append(f"savepdb chain {name}.pdb")
        tleap_script.append(f"saveamberparm chain {name}.prmtop {name}.inpcrd")

    elif "ions" in polymerize_config:

        builtin = []

        for ion in polymerize_config['ions']:
            if ion in typify_dir:
                prepi = typify_dir[ion] / f"{ion}.prepi"
                frcmod = typify_dir[ion] / f"{ion}.frcmod"
                tleap_script.append(f"loadamberprep {prepi.absolute()}")
                tleap_script.append(f"loadamberparams {frcmod.absolute()}")
            else:
                builtin.append(ion)

        tleap_script.append(f"salt = combine {{ {' '.join(polymerize_config['ions'])} }}")
        if len(builtin) == 1:
            tleap_script.append(f"addIons salt {ion} 0")

        tleap_script.append(f"savepdb salt {name}.pdb")
        tleap_script.append(f"saveamberparm salt {name}.prmtop {name}.inpcrd")

    tleap_script.append("quit")

    this_dir = proj_dir / '.db' / name
    this_dir.mkdir(parents=True, exist_ok=True)
    
    with open(this_dir /"tleap.in", "w") as f:
        f.write("\n".join(tleap_script))

    if backend == "intermol":

        yield {
            "cmd": [
                # "ml PDC", "ml amber", 
                "tleap -f tleap.in",
                f"intermol-convert --amb_in {name}.prmtop {name}.inpcrd --lammps"
            ],
            "block": True,
            "cwd": this_dir
        }

    elif backend == "molpy":

        yield {
            "cmd": [
                # "ml PDC", "ml amber", 
                "tleap -f tleap.in"
            ],
            "block": True,
            "cwd": this_dir
        }

        system = mp.io.read_amber(this_dir / f"{name}.prmtop", this_dir / f"{name}.inpcrd")
        system.frame['atoms']['molid'] = system.frame['atoms']['residue']
        mp.io.write_lammps(system, this_dir / f"{name}_converted.lmp", this_dir / f"{name}_converted.input")

    return Path(this_dir)