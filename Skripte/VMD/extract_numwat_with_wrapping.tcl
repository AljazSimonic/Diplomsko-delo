#Skripta najde num_waters (25) obema podenotama najbližjih molekul vode, začenjši pri startframe (0 - od začetka) do endingframe (-1 - do konca), za vsako stride vmesno stanje (vsako 100to vmesno stanje).
#Skripta shrani .pdb (struktura) in .psf (topologija) datoteke v direktorij extractnumwat, za kompleks in vsako podenoto posebaj.
#V topologiji morata podenoti biti poimenovani SEGA in SEGB!
#NAMENJENO NWAT-MM/GBSA analizam,

set dir /d/hpc/home/asimonic/vmd_plugins/pbctools2.8
source ${dir}/pkgIndex.tcl
package require pbctools

#Spremeni glede na svojo strukturo direktorijev
cd ${name}/EpEX_${name}_sim${simn}/

file mkdir [pwd]/extractnumwat
set dir [pwd]/extractnumwat

#Naloži topologijo, vnesi pot do svoje topologije
set trajectory [mol load psf ../01_Input_EpEX_${name}/EpEX_${name}_wbi.psf]

set startframe 0
set endingframe -1
set stride 100

set num_waters 25

set ch 0
set mode 0

for { set i 1 } { $i <= 1 } { incr i } {
    #Spremeni glede na svoje poimenovanje trajektorij
    mol addfile part${i}/EpEX_${name}_part${i}.dcd first $startframe last $endingframe step $stride waitfor all
    puts stdout "part${i}, current numframes [molinfo $trajectory get numframes]"
}


#Spremeni, da bo primerno za tvoj protein! Preberi o wrappingu trajektorij in VMD vtičniku PBSTools!
pbc wrap -center com -centersel "protein and segname SEGA and resid 92" -all -compound fragment -compoundref "protein and segname SEGB and resid 92"

pbc writexst "${dir}/${name}.xst"

set num_frames [molinfo top get numframes]

for {set i 0} {$i < $num_frames} {incr i 1} {
    set dist 2.5

    set sel [atomselect top "water name OH2 and within $dist of (protein segname SEGA and polar) and within $dist of (protein segname SEGB and polar)" frame $i]

    if {[$sel num] > $num_waters} {
        set mode -1
        set ch -0.08
    } elseif {[$sel num] < $num_waters} {
        set mode 1
        set ch 0.08
    } elseif {[$sel num] == $num_waters} {
        set writesel_AB [atomselect top "protein or same fragment as (water name OH2 and within $dist of (protein segname SEGA and polar) and within $dist of (protein segname SEGB and polar))" frame $i]
        set writesel_A [atomselect top "(protein segname SEGA) or same fragment as (water name OH2 and within $dist of (protein segname SEGA and polar) and within $dist of (protein segname SEGB and polar))" frame $i]
        set writesel_B [atomselect top "(protein segname SEGB)" frame $i]
    }
    puts "begin loop"
    while {[$sel num] != $num_waters} {
        $sel delete
        set sel [atomselect top "water name OH2 and within $dist of (protein segname SEGA and polar) and within $dist of (protein segname SEGB and polar)" frame $i]
        if {[$sel num] == $num_waters} {
            set writesel_AB [atomselect top "protein or same fragment as (water name OH2 and within $dist of (protein segname SEGA and polar) and within $dist of (protein segname SEGB and polar))" frame $i]
            set writesel_A [atomselect top "(protein segname SEGA) or same fragment as (water name OH2 and within $dist of (protein segname SEGA and polar) and within $dist of (protein segname SEGB and polar))" frame $i]
            set writesel_B [atomselect top "(protein segname SEGB)" frame $i]
            puts "frame $i distance $dist change $ch"
        } elseif {$mode*[$sel num] > $mode*$num_waters} {
            if {abs($ch) <= 0.01} {
                puts "frame ${i} - can't get correct number of waters with first selection, numwat [$sel num]"
                #increase to first number above $num_waters
                puts "sel list: [$sel list] sel num [$sel num]"
                while {[$sel num] > $num_waters} {
                    $sel delete
                    set dist [expr $dist - 0.01*$mode]
                    set sel [atomselect top "water name OH2 and within $dist of (protein segname SEGA and polar) and within $dist of (protein segname SEGB and polar)"]
                    puts "sel list: [$sel list] sel num [$sel num] dist $dist"
                }
                set sel_equal_dist [atomselect top "name EMPTY_SELECTION"]
                puts "sum sel sel_equal_dist [expr [$sel num] + [$sel_equal_dist num]] dist $dist"
                set dist1 $dist
                while {[$sel num] + [$sel_equal_dist num] < $num_waters} {
                    $sel_equal_dist delete
                    set dist1 [expr $dist1 + 0.01]
                    set sel_equal_dist [atomselect top "water name OH2 and (within $dist1 of (protein segname SEGA and polar) and within $dist1 of (protein segname SEGB and polar)) and not (within $dist of (protein segname SEGA and polar) and within $dist of (protein segname SEGB and polar))" frame $i]
                }
                puts "dist $dist dist1 $dist1"
                puts "sel list [$sel list] sel_equal_dist list [$sel_equal_dist list]"
                set full_wat_sel [atomselect top "same fragment as ((index [concat [$sel list]  [lrange [$sel_equal_dist list] 0 [expr $num_waters - [$sel num] - 1]]]) and within $dist1 of (protein segname SEGA and polar) and within $dist1 of (protein segname SEGB and polar))" frame $i]
                set print_sel [atomselect top "water index [$full_wat_sel list] and name OH2" frame $i]
                puts "full_wat_sel list OH2: [$print_sel list]"
                puts "equal distance: numwat*3 [$full_wat_sel num]"
                if {[$full_wat_sel num] == [expr $num_waters*3]} {
                    set writesel_AB [atomselect top "protein or (same fragment as index [$full_wat_sel list])" frame $i]
                    set writesel_A [atomselect top "(protein segname SEGA) or (same fragment as index [$full_wat_sel list])" frame $i]
                    set writesel_B [atomselect top "(protein segname SEGB)" frame $i]
                    $full_wat_sel delete
                    $sel_equal_dist delete
                    break
                } else {
                    puts "err - numwat still not ok, num [expr [$full_wat_sel num]/3 ]"

                    $full_wat_sel delete
                    $sel_equal_dist delete
                    break
                }
            }
            set mode [expr (-1)*$mode]
            set ch [expr (abs($ch)*$mode)/2]
        }
        set dist [expr $dist + $ch]
    }
    $writesel_AB writepdb "${dir}/${name}_frame_AB[expr ${i}*${stride}].pdb"
    $writesel_AB writepsf "${dir}/${name}_frame_AB[expr ${i}*${stride}].psf"
    $writesel_A writepdb "${dir}/${name}_frame_A[expr ${i}*${stride}].pdb"
    $writesel_A writepsf "${dir}/${name}_frame_A[expr ${i}*${stride}].psf"
    $writesel_B writepdb "${dir}/${name}_frame_B[expr ${i}*${stride}].pdb"
    $writesel_B writepsf "${dir}/${name}_frame_B[expr ${i}*${stride}].psf"
    $sel delete
    $writesel_AB delete
    $writesel_A delete
    $writesel_B delete
}

