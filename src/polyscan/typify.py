from pathlib import Path

from h_submitor import cmdline

@cmdline
def typify(typify_config: dict, proj_dir: Path) -> Path:
    """
    this one should work at db_path/segment_name
    """
    print(typify_config)
    name = typify_config["name"]

    typify_dir = proj_dir / '.db' / name

    if name == "Li":
        return typify_dir

    else:
        net_charge = typify_config["charge"]
        model_path = typify_config["model_path"]
        topo_path = typify_config.get("model_top", None)
        ac_path = typify_dir / Path(f"{name}.ac")
        if not ac_path.exists():
            yield {
                "cmd": [
                    f"antechamber -i {model_path} -fi pdb -o {name}.ac -fo ac -at gaff -an y -c bcc -nc {net_charge} -rn {name}"
                ],
                "block": True,
                "cwd": typify_dir,
            }

        assert ac_path.exists()

        prepi_path = typify_dir / Path(f"{name}.prepi")

        if not prepi_path.exists():
            if topo_path:
                yield {
                    "cmd": [
                        f"prepgen -i {name}.ac -o {name}.prepi -f prepi -m {topo_path} -rn {name} -rf {name}.res"
                    ],
                    "block": True,
                    "cwd": typify_dir,
                }
            else:
                yield {
                    "cmd": [
                        f"prepgen -i {name}.ac -o {name}.prepi -f prepi -rn {name} -rf {name}.res"
                    ],
                    "block": True,
                    "cwd": typify_dir,
                }

        assert prepi_path.exists()

        frcmod_path = typify_dir / Path(f"{name}.frcmod")
        if not frcmod_path.exists():
            yield {
                "cmd": [f"parmchk2 -i {name}.prepi -f prepi -o {name}.frcmod"],
                "block": True,
                "cwd": typify_dir,
            }
        assert frcmod_path.exists()
    return typify_dir