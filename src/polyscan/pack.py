from pathlib import Path
import numpy as np
import molpy as mp

def pack2(data: dict, pack_config: dict, proj_dir: Path, n_trial: int|None = None) -> Path:
    density = pack_config["density"]
    structs = pack_config["structs"]
    seed = pack_config.get("seed", None)
    filename = "_".join([f"{struct['name']}_n{struct['number']}" for struct in structs])
    if n_trial is not None:
        filename += f"_trial{n_trial}"

    this_dir = (
        proj_dir
        / filename
    )
    this_dir.mkdir(parents=True, exist_ok=True)

    # numerical_type = True

    import molpack as mpack
    app = mpack.Molpack(this_dir)

    ff = mp.ForceField("ff")

    n_atoms = 0
    systems = []
    for struct in structs:
        if struct["number"] == 0:
            continue
        name = struct["name"]
        if name in data:
            system = mp.io.read_lammps(
                data[name] / f"{struct['name']}_converted.lmp",
                data[name] / f"{struct['name']}_converted.input",
            )
            systems.append(system)
            n_atoms += len(system.frame["atoms"]) * struct["number"]
        else:
            raise ValueError(f"{name} not found in data")

    for system, struct in zip(systems, structs):
        if "box" not in struct:
            box = mp.Cube([0, 0, 0], (n_atoms / density) ** (1 / 3))
        else:
            box = mp.Cube([0, 0, 0], struct["box"])
        app.add_target(
            system.frame, struct["number"], region=box
        )
            # if numerical_type:
            #     for atomtype in system.forcefield.atomtypes:
            #         atomtype.name = str(int(atomtype.name) + ff.n_atomtypes)
            #     for bondtype in system.forcefield.bondtypes:
            #         bondtype.name = str(int(bondtype.name) + ff.n_bondtypes)
            #     for angletype in system.forcefield.angletypes:
            #         angletype.name = str(int(angletype.name) + ff.n_angletypes)
            #     for dihedraltypes in system.forcefield.dihedraltypes:
            #         dihedraltypes.name = str(int(dihedraltypes.name) + ff.n_dihedraltypes)
            #     for pairtype in system.forcefield.pairtypes:
            #         pairtype.name = str(int(pairtype.name) + ff.n_pairtypes)

        ff.merge_(system.forcefield, offset_type = True)

        
    optimized_system = app.optimize(seed=seed)
    optimized_system.forcefield = ff
    xmax = np.max(optimized_system.frame['atoms']['x'])
    ymax = np.max(optimized_system.frame['atoms']['y'])
    zmax = np.max(optimized_system.frame['atoms']['z'])
    optimized_system.box = mp.Box.cubic([xmax, ymax, zmax], )
    mp.io.write_lammps(optimized_system, this_dir/"system.data", this_dir/"system.ff")
    return this_dir