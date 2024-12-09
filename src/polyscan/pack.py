from pathlib import Path
import numpy as np
import molpy as mp
from copy import deepcopy
from h_submitor import cmdline


@cmdline
def pack(data: dict, pack_config: dict, proj_dir: Path) -> Path:

    structs = pack_config["structs"]
    seed = pack_config.get("seed", None)

    this_dir = (
        proj_dir
        / ".db"
        / "_".join([f"{struct['name']}_n{struct['number']}" for struct in structs])
    )
    this_dir.mkdir(parents=True, exist_ok=True)

    n_atoms = 0
    systems = {}

    for struct in structs:
        if struct["number"] == 0:
            continue
        name = struct["name"]
        if name in data:
            systems[struct["name"]] = mp.io.read_lammps(
                data[name] / f"{struct['name']}_converted.lmp",
                data[name] / f"{struct['name']}_converted.input",
            )
            mp.io.write_pdb(systems[struct["name"]], this_dir / f"{name}.pdb")
        else:
            raise ValueError(f"{name} not found in data")

    lines = []
    lines.append("tolerance 2.0")
    lines.append("filetype pdb")
    lines.append("output _optimized.pdb")
    lines.append("connect no")
    if not seed:
        seed = np.random.randint(100000)
    lines.append(f"seed {seed}")
    for struct in structs:
        name = struct["name"]
        number = struct["number"]
        if number == 0:
            continue
        box = struct["box"]
        lines.append(f"structure {name}.pdb")
        lines.append(f"  connect no")
        lines.append(f"  number {number}")
        lines.append(f"  inside cube 0 0 0 {box}")
        lines.append(f"end structure")

    with open(this_dir / "packmol.inp", "w") as f:
        f.write("\n".join(lines))

    yield {"cmd": ["packmol < packmol.inp"], "block": True, "cwd": this_dir}

    ff = mp.ForceField("ff")
    for system in systems.values():
        f = system.forcefield
        for atomtype in f.atomtypes:
            atomtype.name = str(int(atomtype.name) + ff.n_atomtypes)
        for bondtype in f.bondtypes:
            bondtype.name = str(int(bondtype.name) + ff.n_bondtypes)
        for angletype in f.angletypes:
            angletype.name = str(int(angletype.name) + ff.n_angletypes)
        for dihedraltype in f.dihedraltypes:
            dihedraltype.name = str(int(dihedraltype.name) + ff.n_dihedraltypes)
        for pairtype in f.pairtypes:
            pairtype.name = str(int(pairtype.name) + ff.n_pairtypes)
        ff.merge_(f)

    new_system = mp.System(forcefield=ff)

    optimized_system = mp.io.read_pdb(this_dir / "_optimized.pdb")
    x = optimized_system.frame["atoms"]["x"]
    y = optimized_system.frame["atoms"]["y"]
    z = optimized_system.frame["atoms"]["z"]

    frames = []
    n_atoms_list = [
        0,
    ]
    n_atomtypes = [0]
    n_bondtypes = [0]
    n_angletypes = [0]
    n_dihedraltypes = [0]
    n_structs = 0
    for struct in structs:
        name = struct["name"]
        nrepeat = struct["number"]
        if nrepeat == 0:
            continue

        # systems[name] = mp.io.read_lammps_data(this_dir/f"{name}.data")
        n_atoms = len(systems[name].frame["atoms"])
        n_atomtypes.append(len(np.unique(systems[name].frame["atoms"]["type"])))
        if "bonds" in systems[name].frame:
            n_bondtypes.append(len(np.unique(systems[name].frame["bonds"]["type"])))
        if "angles" in systems[name].frame:
            n_angletypes.append(len(np.unique(systems[name].frame["angles"]["type"])))
        if "dihedrals" in systems[name].frame:
            n_dihedraltypes.append(
                len(np.unique(systems[name].frame["dihedrals"]["type"]))
            )
        for i in range(nrepeat):
            frame = deepcopy(systems[name].frame)
            n_atoms_list.append(n_atoms)
            n_atoms_added = sum(n_atoms_list[:-1])
            struct_optimized_x = x[n_atoms_added : n_atoms_added + n_atoms_list[-1]]
            struct_optimized_y = y[n_atoms_added : n_atoms_added + n_atoms_list[-1]]
            struct_optimized_z = z[n_atoms_added : n_atoms_added + n_atoms_list[-1]]
            n_structs += 1
            frame["atoms"]["x"] = struct_optimized_x.to_numpy()
            frame["atoms"]["y"] = struct_optimized_y.to_numpy()
            frame["atoms"]["z"] = struct_optimized_z.to_numpy()
            frame["atoms"]["molid"] = np.full((n_atoms,), n_structs)

            frame["atoms"]["type"] = frame["atoms"]["type"].to_numpy() + sum(
                n_atomtypes[:-1]
            )  # intermol use unique atom type
            if "bonds" in frame:

                frame["bonds"]["i"] = frame["bonds"]["i"].to_numpy() + n_atoms_added
                frame["bonds"]["j"] = frame["bonds"]["j"].to_numpy() + n_atoms_added
                frame["bonds"]["type"] = frame["bonds"]["type"].to_numpy() + sum(
                    n_bondtypes[:-1]
                )

            if "angles" in frame:
                frame["angles"]["i"] = frame["angles"]["i"].to_numpy() + n_atoms_added
                frame["angles"]["j"] = frame["angles"]["j"].to_numpy() + n_atoms_added
                frame["angles"]["k"] = frame["angles"]["k"].to_numpy() + n_atoms_added
                frame["angles"]["type"] = frame["angles"]["type"].to_numpy() + sum(
                    n_angletypes[:-1]
                )

            if "dihedrals" in frame:
                frame["dihedrals"]["i"] = (
                    frame["dihedrals"]["i"].to_numpy() + n_atoms_added
                )
                frame["dihedrals"]["j"] = (
                    frame["dihedrals"]["j"].to_numpy() + n_atoms_added
                )
                frame["dihedrals"]["k"] = (
                    frame["dihedrals"]["k"].to_numpy() + n_atoms_added
                )
                frame["dihedrals"]["l"] = (
                    frame["dihedrals"]["l"].to_numpy() + n_atoms_added
                )
                frame["dihedrals"]["type"] = frame["dihedrals"][
                    "type"
                ].to_numpy() + sum(n_dihedraltypes[:-1])
            frames.append(frame)

        new_frame = mp.Frame.from_frames(frames)
        new_frame["atoms"]["id"] = np.arange(1, len(new_frame["atoms"]) + 1)
        new_frame["bonds"]["id"] = np.arange(1, len(new_frame["bonds"]) + 1)
        new_frame["angles"]["id"] = np.arange(1, len(new_frame["angles"]) + 1)
        new_frame["dihedrals"]["id"] = np.arange(1, len(new_frame["dihedrals"]) + 1)
        new_system = mp.System(
            forcefield=ff,
            frame=new_frame,
            box=mp.Box(
                np.diag(
                    [
                        max(new_frame["atoms"]["x"].to_numpy())
                        - min(new_frame["atoms"]["x"].to_numpy())
                        + 1,
                        max(new_frame["atoms"]["y"].to_numpy())
                        - min(new_frame["atoms"]["y"].to_numpy())
                        + 1,
                        max(new_frame["atoms"]["z"].to_numpy())
                        - min(new_frame["atoms"]["z"].to_numpy())
                        + 1,
                    ]
                ),
                np.array([True, True, True]),
            ),
        )

    mp.io.write_lammps(new_system, this_dir / "system.data", this_dir / "system.ff")

    return this_dir
