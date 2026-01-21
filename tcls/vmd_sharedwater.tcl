# GET WATERS WITHIN PROTEIN - WITH PBC, BY SECTIONS
set dist 3.6
set num_frames 100000
set start_frame 0
set end_frame [expr {$num_frames - 1}]
# Open output file
set filename "shared_wats3.dat"
set outfile [open $filename w]
# Loop over each frame
for {set frame 0} {$frame < $num_frames} {incr frame} {
    if {[expr $frame % 1000] == 0} {
        puts "Processed $frame frames."
    }
    # Set the current frame
    animate goto $frame
    # Select O atoms within the given criteria
    set O_sel_def [atomselect top "type 349 and pbwithin $dist of (defomed) and type 231 50 1 5 53 234 806 809 820" frame $frame]
    set O_sel_hel [atomselect top "type 349 and pbwithin $dist of (helical) and type 231 50 1 5 53 234 806 809 820" frame $frame]
    set O_sel_prot [atomselect top "type 349 and pbwithin $dist of type 231 50 1 5 53 234 806 809 820" frame $frame]
    set O_sel_glc [atomselect top "type 349 and pbwithin $dist of type 601 603 605 607 609 611" frame $frame]
    
    # Get indices for each selection
    set def_indices [$O_sel_def get index]
    set hel_indices [$O_sel_hel get index]
    set prot_indices [$O_sel_prot get index]
    set glc_indices  [$O_sel_glc get index]

    # Find intersection
    set common_indices_def {}
    foreach idx $def_indices {
        if {[lsearch -exact $glc_indices $idx] != -1} {
            lappend common_indices_def $idx
        }
    }
    # Find intersection
    set common_indices_hel {}
    foreach idx $hel_indices {
        if {[lsearch -exact $glc_indices $idx] != -1} {
            lappend common_indices_hel $idx
        }
    }

    set count_common_def [llength $common_indices_def]
    set count_common_hel [llength $common_indices_hel]
    set count_prot [llength $prot_indices]
    set count_glc [llength $glc_indices]
    set count_def [llength $def_indices]
    set count_hel [llength $hel_indices]

    puts $outfile "$frame $count_prot $count_glc $count_def $count_common_def $count_hel $count_common_hel"

    # Delete oxygen selection
    $O_sel_def delete
    $O_sel_hel delete
    $O_sel_glc delete
    $O_sel_prot delete
}
# Close output file
close $outfile
puts "Finished writing coordinates to $filename"
selection clear
foreach var [info vars] {
    unset $var
}


quit