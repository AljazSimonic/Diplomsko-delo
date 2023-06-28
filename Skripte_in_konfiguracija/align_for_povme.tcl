set trajectory [mol load psf $name/EpEX_${name}_sim${simn}/nowat/EpEX_${name}_protein.psf]

mol addfile ${name}/EpEX_${name}_sim${simn}/nowat/EpEX_${name}_protein_wrapped.dcd step 100 waitfor all

set ref_structure [mol load pdb ${name}/01_Input_EpEX_${name}/EpEX_${name}_SEGA.pdb]

file mkdir ${name}/EpEX_${name}_sim${simn}/nowat/aligned_structures

set nf [molinfo $trajectory get numframes]

set sel [atomselect $trajectory "protein segname SEGA"]
set ref [atomselect $ref_structure "protein"]
for { set i 0 } { $i < $nf } { incr i 1 } {
    $sel frame $i
    $sel move [measure fit $sel $ref]
    $sel writepdb "${name}/EpEX_${name}_sim${simn}/nowat/aligned_structures/SEGA_frame[expr ${i}*100].pdb"
}

set sel [atomselect $trajectory "protein segname SEGB"]
set ref [atomselect $ref_structure "protein"] ;# Trajectory segments SEGA and SEGB are both aligned to reference SEGA
for { set i 0 } { $i < $nf } { incr i 1 } {
    $sel frame $i
    $sel move [measure fit $sel $ref]
    $sel writepdb "${name}/EpEX_${name}_sim${simn}/nowat/aligned_structures/SEGB_frame[expr ${i}*100].pdb"
}
