from pathlib import Path
from string import Template


def write_mech_script(eq_dir: Path, T:int, rate:float) -> Path:

    import random

    init_temp = """# youngs: init
units real
atom_style full

dimension 3
bond_style harmonic
angle_style harmonic
dihedral_style hybrid charmm multi/harmonic
special_bonds lj 0.0 0.0 0.5 coul 0.0 0.0 0.83333333
pair_style lj/cut/coul/long 9.0 9.0
pair_modify tail yes
kspace_style pppm 1e-6

read_data ../${T}K.data
change_box all triclinic
include system.ff
kspace_style pppm 1.0e-5

velocity all create ${T} ${rndint}

fix 1 all nvt/sllod temp ${T} ${T} 100.0
"""

    stemp = Template(init_temp)
    result = stemp.safe_substitute({
        'T': T,
        'rndint': random.randint(0, 10000),
    })
    with open(eq_dir/'youngs_init.in', 'w') as f:
        f.write(result)

    template = """# youngs
variable        tkrun equal 500
variable        def_rate equal ${rate}
variable        tk equal "temp"
variable        vl equal "vol"
variable        dn equal "density"
## Converting [atm] = 101325 [pa] = 101325*1e-6 [MPa]
variable        p2 equal "-pxx*101325*1e-6"
variable        p3 equal "-pyy*101325*1e-6"
variable        p4 equal "-pzz*101325*1e-6"
variable        p5 equal "-pyz*101325*1e-6"
variable        p6 equal "-pxz*101325*1e-6"
variable        p7 equal "-pxy*101325*1e-6"
variable        p8 equal "lx"
variable        p9 equal "ly"
variable        p10 equal "lz"
variable        p11 equal "yz"
variable        p12 equal "xz"
variable        p13 equal "xy"
variable        sdt equal 10  ##1000
variable        sdtdp equal 20000

# Uniaxial Tensile Deformation X-direction
include         youngs_init.in
print           "Deformation X-direction"
variable        tmp equal "lx"
variable        L0 equal ${tmp}
variable        strain equal "(lx - v_L0)/v_L0"
variable        p1 equal "v_strain"
fix             2 all deform 1 x erate ${def_rate} units box remap v
fix             def1 all print ${sdt} "${p1} ${p2} ${p3} ${p4} ${p5} ${p6} ${p7} ${p8} ${p9} ${p10} ${p11} ${p12} ${p13} ${tk}" file cdefx_${tkrun}K.txt screen no title "#Strain     -Pxx [MPa]   -Pyy [MPa]     -Pzz [MPa]     -Pyz [MPa]   -Pxz [MPa]     -Pxy [MPa]    lxx [A]    lyy [A]    lzz [A]   ayz [A]    axz [A]    axy [A]   T [K]  ;  "
dump           mydump all atom ${sdtdp} dump.cdefx_${tkrun}K.lammpstrj
reset_timestep  0
run             ${step}
write_data      cdefx_${tkrun}K.data
clear

# Uniaxial Tensile Deformation Y-direction
include         youngs_init.in

print           "Deformation Y-direction" 
variable        tmp equal "ly"
variable        L0 equal ${tmp}
variable        strain equal "(ly - v_L0)/v_L0"
variable        p1 equal "v_strain"
fix             2 all deform 1 y erate ${def_rate} units box remap v
fix             def1 all print ${sdt} "${p1} ${p2} ${p3} ${p4} ${p5} ${p6} ${p7} ${p8} ${p9} ${p10} ${p11} ${p12} ${p13} ${tk}" file cdefy_${tkrun}K.txt screen no title "#Strain     -Pxx [MPa]   -Pyy [MPa]     -Pzz [MPa]     -Pyz [MPa]   -Pxz [MPa]     -Pxy [MPa]    lxx [A]    lyy [A]    lzz [A]   ayz [A]    axz [A]    axy [A]   T [K]  ;  "
#dump            mydump all atom ${sdtdp} dump.cdefy_${tkrun}K.lammpstrj
reset_timestep  0
run             ${step}
write_data      cdefy_${tkrun}K.data
clear

# Uniaxial Tensile Deformation Z-direction
include         youngs_init.in

print           "Deformation Z-direction" 
variable        tmp equal "lz"
variable        L0 equal ${tmp}
variable        strain equal "(lz - v_L0)/v_L0"
variable        p1 equal "v_strain"
fix             2 all deform 1 z erate ${def_rate} units box remap v
fix             def1 all print ${sdt} "${p1} ${p2} ${p3} ${p4} ${p5} ${p6} ${p7} ${p8} ${p9} ${p10} ${p11} ${p12} ${p13} ${tk}" file cdefz_${tkrun}K.txt screen no title "#Strain     -Pxx [MPa]   -Pyy [MPa]     -Pzz [MPa]     -Pyz [MPa]   -Pxz [MPa]     -Pxy [MPa]    lxx [A]    lyy [A]    lzz [A]   ayz [A]    axz [A]    axy [A]   T [K]  ;  "
#dump            mydump all atom ${sdtdp} dump.cdefz_${tkrun}K.lammpstrj
reset_timestep  0
run             ${step}
write_data      cdefz_${tkrun}K.data
clear

# Shear Deformation YZ-direction
include         youngs_init.in

print           "Deformation YZ-direction" 
variable        tmp equal "lz"
variable        L0 equal ${tmp}
variable        strain equal "yz/v_L0"
variable        p1 equal "v_strain"
fix             2 all deform 1 yz erate ${def_rate} units box remap v
fix             def1 all print ${sdt} "${p1} ${p2} ${p3} ${p4} ${p5} ${p6} ${p7} ${p8} ${p9} ${p10} ${p11} ${p12} ${p13} ${tk}" file cdefyz_${tkrun}K.txt screen no title "#Strain     -Pxx [MPa]   -Pyy [MPa]     -Pzz [MPa]     -Pyz [MPa]   -Pxz [MPa]     -Pxy [MPa]    lxx [A]    lyy [A]    lzz [A]   ayz [A]    axz [A]    axy [A]   T [K]  ;  "
#dump            mydump all atom ${sdtdp} dump.cdefyz_${tkrun}K.lammpstrj
reset_timestep  0
run             ${step}
write_data      cdefyz_${tkrun}K.data
clear

# Shear Deformation XZ-direction
include         youngs_init.in

print           "Deformation XZ-direction" 
variable        tmp equal "lz"
variable        L0 equal ${tmp}
variable        strain equal "xz/v_L0"
variable        p1 equal "v_strain"
fix             2 all deform 1 xz erate ${def_rate} units box remap v
fix             def1 all print ${sdt} "${p1} ${p2} ${p3} ${p4} ${p5} ${p6} ${p7} ${p8} ${p9} ${p10} ${p11} ${p12} ${p13} ${tk}" file cdefxz_${tkrun}K.txt screen no title "#Strain     -Pxx [MPa]   -Pyy [MPa]     -Pzz [MPa]     -Pyz [MPa]   -Pxz [MPa]     -Pxy [MPa]    lxx [A]    lyy [A]    lzz [A]   ayz [A]    axz [A]    axy [A]   T [K]  ;  "
#dump            mydump all atom ${sdtdp} dump.cdefxz_${tkrun}K.lammpstrj
reset_timestep  0
run             ${step}
write_data      cdefxz_${tkrun}K.data
clear

# Shear Deformation XY-direction
include         youngs_init.in

print           "Deformation XY-direction" 
variable        tmp equal "ly"
variable        L0 equal ${tmp}
variable        strain equal "xy/v_L0"
variable        p1 equal "v_strain"
fix             2 all deform 1 xy erate ${def_rate} units box remap v
fix             def1 all print ${sdt} "${p1} ${p2} ${p3} ${p4} ${p5} ${p6} ${p7} ${p8} ${p9} ${p10} ${p11} ${p12} ${p13} ${tk}" file cdefxy_${tkrun}K.txt screen no title "#Strain     -Pxx [MPa]   -Pyy [MPa]     -Pzz [MPa]     -Pyz [MPa]   -Pxz [MPa]     -Pxy [MPa]    lxx [A]    lyy [A]    lzz [A]   ayz [A]    axz [A]    axy [A]   T [K]  ;  "
#dump            mydump all atom ${sdtdp} dump.cdefxy_${tkrun}K.lammpstrj
reset_timestep  0
run             ${step}
write_data      cdefxy_${tkrun}K.data
clear

print           "All done" 

"""
    ##rate 1e10 s-1 = 1e-5 (fs)-1 lammps erate
    ##rate 1e9  s-1 = 1e-6 (fs)-1 lammps erate
    ##rate 1e8  s-1 = 1e-7 (fs)-1 lammps erate

    lammps_rate = rate * 1e-15
    step = int(0.2 / lammps_rate)
    stemp = Template(template)
    result = stemp.safe_substitute({
        'T': T,
        'rndint': random.randint(0, 10000),
        'rate': lammps_rate,
        'step': step
    })
    with open(eq_dir/'youngs.in', 'w') as f:
        f.write(result)