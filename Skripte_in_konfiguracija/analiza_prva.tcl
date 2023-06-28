
set trajectory [mol load psf $name/EpEX_${name}_sim${simn}/nowat/EpEX_${name}_protein.psf]

mol addfile $name/EpEX_${name}_sim${simn}/nowat/EpEX_${name}_protein_wrapped.dcd waitfor all


set nf [molinfo $trajectory get numframes]
puts stdout "Numframes: $nf"

[atomselect $trajectory "protein" frame 0] writepdb $name/EpEX_${name}_sim${simn}/nowat/complex_firstframe.pdb
[atomselect $trajectory "protein segname SEGA" frame 0] writepdb $name/EpEX_${name}_sim${simn}/nowat/SEGA_firstframe.pdb
[atomselect $trajectory "protein segname SEGB" frame 0] writepdb $name/EpEX_${name}_sim${simn}/nowat/SEGB_firstframe.pdb

set outfile_sasa [open ${name}/analiza/sasa_sim${simn}.dat w]
set outfile_rgyr [open ${name}/analiza/rgyr_sim${simn}.dat w]
set nf [molinfo $trajectory get numframes]
set sel [atomselect top "protein"]
# sasa calculation loop

for { set i 0 } { $i < $nf } { incr i 5 } {
    $sel frame $i
    puts $outfile_sasa "[measure sasa 1.4 $sel]"
    puts $outfile_rgyr "[measure rgyr $sel]"
    puts stdout "SASA and rgyr $i/$nf"
}
close $outfile_sasa
close $outfile_rgyr
$sel delete

set outfile [open ${name}/analiza/rmsd_sim${simn}.dat w]
set ref [atomselect top "protein and name CA" frame 0]
set sel [atomselect top "protein and name CA"]
for { set f 0 } { $f < $nf } { incr f } {
    $sel frame $f
    $sel move [measure fit $sel $ref]
    puts $outfile "[measure rmsd $sel $ref]"
}
close $outfile
$sel delete
$ref delete


set outfile [open ${name}/analiza/rmsf_SEGA_sim${simn}.dat_spaces w]
set ref [atomselect top "protein and name CA and segname SEGA" frame 0]
set sel [atomselect top "protein and name CA and segname SEGA"]
for { set f 0 } { $f < $nf } { incr f } {
    $sel frame $f
    $sel move [measure fit $sel $ref]
}
puts $outfile "[measure rmsf $sel]"
close $outfile
$sel delete
$ref delete

set outfile [open ${name}/analiza/rmsf_SEGB_sim${simn}.dat_spaces w]
set ref [atomselect top "protein and name CA and segname SEGB" frame 0]
set sel [atomselect top "protein and name CA and segname SEGB"]
for { set f 0 } { $f < [molinfo top get numframes] } { incr f } {
    $sel frame $f
    $sel move [measure fit $sel $ref]
}
puts $outfile "[measure rmsf $sel]"
close $outfile
$sel delete
$ref delete

exit
