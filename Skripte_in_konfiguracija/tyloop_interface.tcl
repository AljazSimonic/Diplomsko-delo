set trajectory [mol load psf $name/EpEX_${name}_sim${simn}/nowat/EpEX_${name}_protein.psf]

mol addfile $name/EpEX_${name}_sim${simn}/nowat/EpEX_${name}_protein_wrapped.dcd waitfor all

set nf [molinfo $trajectory get numframes]

set outfile_A [open ${name}/analiza/interface_tyloop_A_sim${simn}.dat w]
set outfile_B [open ${name}/analiza/interface_tyloop_B_sim${simn}.dat w]

set restrictsel [atomselect top "protein segname SEGA and resid ${sel_begin} to ${sel_end}"]
set sel_with_other_subunit [atomselect top "protein and (segname SEGA and resid ${sel_begin} to ${sel_end}) or (segname SEGB and resid > 130)"]

for { set i 0 } { $i < $nf } { incr i 5 } {
    $restrictsel frame $i
    $sel_with_other_subunit frame $i
    puts $outfile_A "[expr [measure sasa 1.4 $restrictsel] - [measure sasa 1.4 $sel_with_other_subunit -restrict $restrictsel]]"
    puts "A $i"

}
close $outfile_A

set restrictsel [atomselect top "protein segname SEGB and resid ${sel_begin} to ${sel_end}"]
set sel_with_other_subunit [atomselect top "protein and (segname SEGB and resid ${sel_begin} to ${sel_end}) or (segname SEGA and resid > 130)"]

for { set i 0 } { $i < $nf } { incr i 5 } {
    $restrictsel frame $i
    $sel_with_other_subunit frame $i
    puts $outfile_B "[expr [measure sasa 1.4 $restrictsel] - [measure sasa 1.4 $sel_with_other_subunit -restrict $restrictsel]]"
    puts "B $i"
}
close $outfile_B

exit
