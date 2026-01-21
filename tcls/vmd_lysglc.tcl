# TOTAL VMD SCRIPT FOR LYS-GLC RDF ONLY

# Find Glucose-LYS HBonds 3.5 A
hbonds -ang 30 -dist 3.5 -writefile yes -sel2 [atomselect top "type 601 or type 603 or type 605 or type 607 or type 609 or type 611"] -sel1 [atomselect top "type 190"] -outfile analysis/hbond_gluc_lys_3.5.dat -plot no
selection clear
foreach var [info vars] {
    unset $var
}

# Find Water-LYS HBonds 3.0 A
hbonds -ang 30 -dist 3.0 -writefile yes -sel2 [atomselect top "type 349 or type 350"] -sel1 [atomselect top "type 190"] -outfile analysis/hbond_lys_wat_3.0.dat -plot no
selection clear
foreach var [info vars] {
    unset $var
}

# Find Glucose-LYS RDF
set num_frames 100000
set start_frame 0
set end_frame [expr {$num_frames - 1}]
set ringo [atomselect top "type 601 or type 603 or type 605 or type 607 or type 609 or type 611"]
set wats [atomselect top "type 190"]
set pg [measure gofr $ringo $wats delta 0.1 rmax 15.0 first $start_frame last $end_frame step 1 usepbc true]
set outfile [open "analysis/rdf_lys_glc.dat" w]
set r_values [lindex $pg 0]
set rdf_values [lindex $pg 1]
set int_values [lindex $pg 2]
for {set i 0} {$i < [llength $r_values]} {incr i} {
    set r [lindex $r_values $i]
    set rdf [lindex $rdf_values $i]
    set int [lindex $int_values $i]
    puts $outfile [format "%s %s %s" $r $rdf $int]
}
close $outfile

# Find Water-LYS RDF
set num_frames 100000
set start_frame 0
set end_frame [expr {$num_frames - 1}]
set ringo [atomselect top "type 349"]
set wats [atomselect top "type 190"]
set pg [measure gofr $ringo $wats delta 0.1 rmax 15.0 first $start_frame last $end_frame step 1 usepbc true]
set outfile [open "analysis/rdf_lys_wat.dat" w]
set r_values [lindex $pg 0]
set rdf_values [lindex $pg 1]
set int_values [lindex $pg 2]
for {set i 0} {$i < [llength $r_values]} {incr i} {
    set r [lindex $r_values $i]
    set rdf [lindex $rdf_values $i]
    set int [lindex $int_values $i]
    puts $outfile [format "%s %s %s" $r $rdf $int]
}
close $outfile

# GET GLUCOSES WITHIN 5 A OF LYS
set strand 36
set dist 5
set num_frames 100000
set start_frame 0
set end_frame [expr {$num_frames - 1}]
# Open output file
set filename "glcs_lys_[string map { . _ } $dist].dat"
set outfile [open $filename w]
# Loop over each frame
for {set frame 0} {$frame < $num_frames} {incr frame} {
    if {[expr $frame % 1000] == 0} {
        puts "Processed $frame frames."
    }
    # Set the current frame
    animate goto $frame
    # Select O atoms within the given criteria
    set O_sel [atomselect top "type 600 and within $dist of type 190" frame $frame]
    # Get indices, types, coordinates, and bonds of selected O atoms
    set O_indices [$O_sel get index]
    set O_types [$O_sel get type]
    set O_coords [$O_sel get {x y z}]
    set O_bonds [$O_sel getbonds]
    set length_O [llength $O_indices]
    puts $outfile "$length_O"
    # Iterate over selected O atoms
    for {set i 0} {$i < $length_O} {incr i} {
        set O_idx [lindex $O_indices $i]
        set O_coord [lindex $O_coords $i]
        set O_type [lindex $O_types $i]
        set bonds [lindex $O_bonds $i]
        # Print O atom index, type, and coordinates
        set x [lindex $O_coord 0]
        set y [lindex $O_coord 1]
        set z [lindex $O_coord 2]
        puts $outfile "$O_idx $x $y $z $O_type"
    }
    # Delete oxygen selection
    $O_sel delete
}
# Close output file
close $outfile
puts "Finished writing coordinates to $filename"
selection clear
foreach var [info vars] {
    unset $var
}


# GET WATERS WITHIN 3.5 A OF LYS
set strand 36
set dist 3.5
set num_frames 100000
set start_frame 0
set end_frame [expr {$num_frames - 1}]
# Open output file
set filename "wats_lys_[string map { . _ } $dist].dat"
set outfile [open $filename w]
# Loop over each frame
for {set frame 0} {$frame < $num_frames} {incr frame} {
    if {[expr $frame % 1000] == 0} {
        puts "Processed $frame frames."
    }
    # Set the current frame
    animate goto $frame
    # Select O atoms within the given criteria
    set O_sel [atomselect top "type 349 and within $dist of type 190" frame $frame]
    # Get indices, types, coordinates, and bonds of selected O atoms
    set O_indices [$O_sel get index]
    set O_types [$O_sel get type]
    set O_coords [$O_sel get {x y z}]
    set O_bonds [$O_sel getbonds]
    set length_O [llength $O_indices]
    puts $outfile "$length_O"
    # Iterate over selected O atoms
    for {set i 0} {$i < $length_O} {incr i} {
        set O_idx [lindex $O_indices $i]
        set O_coord [lindex $O_coords $i]
        set O_type [lindex $O_types $i]
        set bonds [lindex $O_bonds $i]
        # Print O atom index, type, and coordinates
        set x [lindex $O_coord 0]
        set y [lindex $O_coord 1]
        set z [lindex $O_coord 2]
        puts $outfile "$O_idx $x $y $z $O_type"
        # Iterate over bonded atoms and find H atoms
        foreach bonded_atom $bonds {
            set bonded_sel [atomselect top "index $bonded_atom" frame $frame]
            set bonded_name [$bonded_sel get name]
            set bonded_type [$bonded_sel get type]
            # Check if bonded atom is a hydrogen
            if {$bonded_name == "H"} {
                set H_coords [$bonded_sel get {x y z}]
                set x [lindex $H_coords 0]
                set y [lindex $H_coords 1]
                set z [lindex $H_coords 2]
                puts $outfile "$bonded_atom $x $y $z $bonded_type"
            }
            $bonded_sel delete
        }
    }
    # Delete oxygen selection
    $O_sel delete
}
# Close output file
close $outfile
puts "Finished writing coordinates to $filename"
selection clear
foreach var [info vars] {
    unset $var
}


quit