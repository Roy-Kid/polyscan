from pathlib import Path
from string import Template


def write_md_script(dir, config: dict):

    with open(Path(dir) / "md_eq.in", "w") as f:

        TEMPLATE = """# md_eq.in
units real
atom_style full

dimension 3
boundary p p p

bond_style harmonic
angle_style harmonic
dihedral_style hybrid charmm multi/harmonic
special_bonds lj 0.0 0.0 0.5 coul 0.0 0.0 0.83333333
pair_style lj/cut/coul/long 9.0 9.0
pair_modify tail yes
kspace_style pppm 1e-5

read_data ${data_dir}/system.data
include ${data_dir}/system.ff

velocity all create 1000 8420

timestep 1.0

fix sk all shake 0.001 20 0 m 1.008

variable elapsed equal "elapsed"
variable tk equal "temp"
variable pr equal "press"
variable vl equal "vol"
variable dn equal "density"
variable et equal "etotal"
variable ep equal "pe"
variable ek equal "ke"
variable evdwl equal "evdwl"
variable ecoul equal "ecoul"
variable epair equal "epair"
variable eb equal "ebond"
variable eang equal "eangle"
variable edih equal "edihed"

thermo 10000
thermo_style custom elapsed temp press vol density etotal
fix fprint all print 1000 "${elapsed} ${tk} ${pr} ${vl} ${dn} ${et} ${ep} ${ek}" file eq.thermo screen no
"""
        template = Template(TEMPLATE)
        f.write(template.safe_substitute(data_dir=config["data_dir"].absolute()))
        for step in config["stages"]:
            if step["type"] == "minimize":
                f.write("minimize 1.0e-4 1.0e-6 1000 10000\n")
            elif step["type"] == "npt":
                f.write(
                    f"fix 1 all npt temp {step['Tstart']} {step['Tend']} 100.0 iso 1.0 1.0 1000 drag 2.0\n"
                )
                f.write(f"run {int(step['steps'])}\n")
                f.write(f"write_data T{step['Tend']}K.data nocoeff\n")
            elif step["type"] == "nvt":
                f.write(f"fix 1 all nvt temp {step['Tstart']} {step['Tend']} 100.0\n")
                f.write(f"run {int(step['steps'])}\n")
                f.write(f"write_data T{step['Tend']}K.data nocoeff\n")
