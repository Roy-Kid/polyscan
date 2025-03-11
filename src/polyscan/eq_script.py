from pathlib import Path
from string import Template


def write_tg_script(dir, config: dict):

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
pair_style lj/cut/coul/long 2.5 9.0
pair_modify tail yes

read_data system.data
include system.ff
kspace_style pppm 1e-5

timestep 1.0

variable step equal "step"
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

group litfsi molecule > 20
variable newcharge atom q*0.75
set group litfsi charge v_newcharge

thermo 10000
thermo_style custom elapsed temp press vol density etotal

minimize 1.0e-4 1.0e-6 1000 10000

velocity all create 1000 55623

fix sk all shake 0.001 20 0 m 1.008

fix 1 all nve/limit 0.5
run 10000
unfix 1

fix 1 all npt temp 1000.0 1000.0 100.0 iso 0.0 1.0 1000 drag 2.0
run 2000000 every 10000 "write_data T1000K.data nocoeff"
run 1000000 every 100000 "write_data T1000K.data nocoeff"
unfix 1

fix fprint all print 1000 "${step} ${tk} ${pr} ${vl} ${dn} ${et} ${ep} ${ek}" file tg.thermo screen no

fix 1 all npt temp 1000 500 100.0 iso 1.0 1.0 1000 drag 2.0
run 5000000
write_data 500K.data nocoeff

fix 1 all npt temp 500 100 100.0 iso 1.0 1.0 1000 drag 2.0
run 4000000
write_data 100K.data nocoeff

"""
        template = Template(TEMPLATE)
        f.write(template.safe_substitute(prep_dir=config["pack_dir"].absolute()))

def write_eq_script(dir, config: dict):

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
pair_style lj/cut/coul/long 2.5 9.0
pair_modify tail yes

read_data system.data
include system.ff
kspace_style pppm 1e-4

timestep 1.0

variable step equal "step"
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

minimize 1.0e-4 1.0e-6 1000 10000

velocity all create 1000 55623

fix sk all shake 0.001 20 0 m 1.008

fix 1 all nve/limit 0.5
run 10000
unfix 1

fix 1 all npt temp 1000.0 1000.0 100.0 iso 1.0 1.0 1000 drag 2.0
run 100000 every 100 "write_data T1000K.data nocoeff"
run 1000000 every 10000 "write_data T1000K.data nocoeff"
run 1000000 every 100000 "write_data T1000K.data nocoeff"
unfix 1

fix fprint all print 1000 "${step} ${tk} ${pr} ${vl} ${dn} ${et} ${ep} ${ek}" file tg.thermo screen no

fix 1 all npt temp 1000 500 100.0 iso 1.0 1.0 1000 drag 2.0
run 5000000
write_data 500K.data nocoeff

"""
        template = Template(TEMPLATE)
        f.write(template.safe_substitute(prep_dir=config["pack_dir"].absolute()))