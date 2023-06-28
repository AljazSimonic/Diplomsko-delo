proc ladd {l} {
    set total 0
    foreach nxt $l {
        set total [expr $total + $nxt]
    }
    return $total
}

set trajectory [mol load psf $name/EpEX_${name}_sim${simn}/nowat/EpEX_${name}_protein.psf]

mol addfile $name/EpEX_${name}_sim${simn}/nowat/EpEX_${name}_protein_wrapped.dcd waitfor all


set nf [molinfo $trajectory get numframes]
puts stdout "Numframes: $nf"

set outfile_dihed_rmsf_phi_A [open ${name}/analiza/dihed_rmsf_phi_SEGA_sim${simn}.dat w]

#rmsf_phi_A
set resid_list [lsort -unique -integer [[atomselect top "protein segname SEGA"] get resid]]
for { set i 1 } { $i < [llength $resid_list] - 2 } { incr i 1 } {
    set sel_dihed_phi_A [atomselect top "protein segname SEGA and ((resid [lindex $resid_list [expr $i - 1]] and name C) or (resid [lindex $resid_list $i] and (name N or name CA or name C)))"]
    set phi_over_time [measure dihed [$sel_dihed_phi_A list] frame all]
    set phi_average [expr [ladd $phi_over_time]/$nf]
    set square_deviation_sum 0
    foreach x $phi_over_time {
        set square_deviation_sum [expr $square_deviation_sum + ($x - $phi_average)*($x - $phi_average)]
    }
    puts $outfile_dihed_rmsf_phi_A "[expr sqrt($square_deviation_sum/$nf)]"
    puts stdout "SEGA phi rmsf [lindex $resid_list $i]: [expr sqrt($square_deviation_sum/$nf)]"
    $sel_dihed_phi_A delete
}
close $outfile_dihed_rmsf_phi_A

set outfile_dihed_rmsf_phi_B [open ${name}/analiza/dihed_rmsf_phi_SEGB_sim${simn}.dat w]

#rmsf_phi_B
for { set i 1 } { $i < [llength $resid_list] - 2 } { incr i 1 } {
    set sel_dihed_phi_B [atomselect top "protein segname SEGB and ((resid [lindex $resid_list [expr $i - 1]] and name C) or (resid [lindex $resid_list $i] and (name N or name CA or name C)))"]
    set phi_over_time [measure dihed [$sel_dihed_phi_B list] frame all]
    set phi_average [expr [ladd $phi_over_time]/$nf]
    set square_deviation_sum 0
    foreach x $phi_over_time {
        set square_deviation_sum [expr $square_deviation_sum + ($x - $phi_average)*($x - $phi_average)]
    }
    puts $outfile_dihed_rmsf_phi_B "[expr sqrt($square_deviation_sum/$nf)]"
    puts stdout "SEGB phi rmsf [lindex $resid_list $i]: [expr sqrt($square_deviation_sum/$nf)]"
    $sel_dihed_phi_B delete
}
close $outfile_dihed_rmsf_phi_B

#rmsf_psi_A

set outfile_dihed_rmsf_psi_A [open ${name}/analiza/dihed_rmsf_psi_SEGA_sim${simn}.dat w]

for { set i 1 } { $i < [llength $resid_list] - 2 } { incr i 1 } {
    set sel_dihed_psi_A [atomselect top "protein segname SEGA and ((resid [lindex $resid_list $i] and (name N or name CA or name C)) or (resid [lindex $resid_list [expr $i + 1]] and name N))"]
    set psi_over_time [measure dihed [$sel_dihed_psi_A list] frame all]
    set psi_average [expr [ladd $psi_over_time]/$nf]
    set square_deviation_sum 0
    foreach x $psi_over_time {
        set square_deviation_sum [expr $square_deviation_sum + ($x - $psi_average)*($x - $psi_average)]
    }
    puts $outfile_dihed_rmsf_psi_A "[expr sqrt($square_deviation_sum/$nf)]"
    puts stdout "SEGA psi rmsf [lindex $resid_list $i]: [expr sqrt($square_deviation_sum/$nf)]"
    $sel_dihed_psi_A delete
}
close $outfile_dihed_rmsf_psi_A

set outfile_dihed_rmsf_psi_B [open ${name}/analiza/dihed_rmsf_psi_SEGB_sim${simn}.dat w]

#rmsf_psi_B
for { set i 1 } { $i < [llength $resid_list] - 2 } { incr i 1 } {
    set sel_dihed_psi_B [atomselect top "protein segname SEGB and ((resid [lindex $resid_list $i] and (name N or name CA or name C)) or (resid [lindex $resid_list [expr $i + 1]] and name N))"]
    set psi_over_time [measure dihed [$sel_dihed_psi_B list] frame all]
    set psi_average [expr [ladd $psi_over_time]/$nf]
    set square_deviation_sum 0
    foreach x $psi_over_time {
        set square_deviation_sum [expr $square_deviation_sum + ($x - $psi_average)*($x - $psi_average)]
    }
    puts $outfile_dihed_rmsf_psi_B "[expr sqrt($square_deviation_sum/$nf)]"
    puts stdout "SEGB psi rmsf [lindex $resid_list $i]: [expr sqrt($square_deviation_sum/$nf)]"
    $sel_dihed_psi_B delete
}
close $outfile_dihed_rmsf_psi_B

set outfile_sasa_A [open ${name}/analiza/sasa_SEGA_sim${simn}.dat w]
set nf [molinfo $trajectory get numframes]
set sel [atomselect top "protein segname SEGA"]
# sasa SEGA calculation loop

for { set i 0 } { $i < $nf } { incr i 5 } {
    $sel frame $i
    puts $outfile_sasa_A "[measure sasa 1.4 $sel]"
    puts stdout "SASA SEGA $i/$nf"
}

close $outfile_sasa_A
$sel delete

set outfile_sasa_B [open ${name}/analiza/sasa_SEGB_sim${simn}.dat w]
set nf [molinfo $trajectory get numframes]
set sel [atomselect top "protein segname SEGB"]
# sasa SEGA calculation loop

for { set i 0 } { $i < $nf } { incr i 5 } {
    $sel frame $i
    puts $outfile_sasa_B "[measure sasa 1.4 $sel]"
    puts stdout "SASA SEGB $i/$nf"
}

close $outfile_sasa_B
$sel delete

exit
