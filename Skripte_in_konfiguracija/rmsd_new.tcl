set trajectory [mol load psf $name/EpEX_${name}_sim${simn}/nowat/EpEX_${name}_protein.psf]

mol addfile $name/EpEX_${name}_sim${simn}/nowat/EpEX_${name}_protein_wrapped.dcd waitfor all

set ref_structure [mol load pdb ${name}/01_Input_EpEX_${name}/EpEX_${name}.pdb]

set nf [molinfo $trajectory get numframes]

set outfile [open ${name}/analiza/rmsd_SEGA_sim${simn}.dat w]
set ref [atomselect $ref_structure "protein and segname SEGA and name CA" frame 0]
set sel [atomselect $trajectory "protein and segname SEGA and name CA"]
for { set f 0 } { $f < $nf } { incr f } {
    $sel frame $f
    $sel move [measure fit $sel $ref]
    puts $outfile "[measure rmsd $sel $ref]"
}
close $outfile
$sel delete
$ref delete

set outfile [open ${name}/analiza/rmsd_SEGB_sim${simn}.dat w]
set ref [atomselect $ref_structure "protein and segname SEGB and name CA" frame 0]
set sel [atomselect $trajectory "protein and segname SEGB and name CA"]
for { set f 0 } { $f < $nf } { incr f } {
    $sel frame $f
    $sel move [measure fit $sel $ref]
    puts $outfile "[measure rmsd $sel $ref]"
}
close $outfile
$sel delete
$ref delete
