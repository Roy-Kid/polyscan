from pathlib import Path
from string import Template


def write_md_script(dir, config: dict):

    with open(Path(dir) / "md_eq.in", "w") as f:

        TEMPLATE = """# md_eq.in
units real
atom_style full

dimension 3
boundary p p p

read_data ${prep_dir}/system.data
include ${prep_dir}/system.ff

velocity all create 1000 8420

timestep 1.0

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
        f.write(template.safe_substitute(prep_dir=config["prep_dir"].absolute()))
        for step in config["stages"]:
            
            if step["type"] == "minimize":
                f.write("minimize 1.0e-4 1.0e-6 1000 10000\n")
            elif step["type"] == "npt":
                f.write(
                    f"fix 1 all npt temp {step['Tstart']} {step['Tend']} 100.0 iso 1.0 1.0 1000 drag 2.0\n"
                )
                f.write(f"run {int(step['steps'])}\n")
                f.write(f"write_data T{step['Tend']}K.data nocoeff\n")
                f.write(f"unfix 1\n")
            elif step["type"] == "nvt":
                f.write(f"fix 1 all nvt temp {step['Tstart']} {step['Tend']} 100.0\n")
                f.write(f"run {int(step['steps'])}\n")
                f.write(f"write_data T{step['Tend']}K.data nocoeff\n")
                f.write(f"unfix 1\n")
            elif step["type"] == "nve/limit":
                f.write(f"fix 1 all nve/limit 0.1")
                f.write(f"unfix 1\n")
                f.write(f"fix sk all shake 0.001 20 0 m 1.008")
