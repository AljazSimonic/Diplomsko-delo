set inputmodel EpEX_hsap

mol delete all
mol load pdb ${inputmodel}_eq_protFIXED_lastframe.pdb

# measure the system

set file [open "${inputmodel}_newdimensions.txt" w]
set everyone [atomselect top all]
puts $file [ measure minmax $everyone ]
close $file

set file [open "${inputmodel}_newcenter.txt" w]
set everyone [atomselect top all]
puts $file [ measure center $everyone ]
close $file
