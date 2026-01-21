# HBONDS BETWEEN GLY AND XAA PRO
hbonds -ang 30 -dist 3.0 -writefile yes -sel1 [atomselect top "type 1 and (xxxsel)"] -sel2 [atomselect top "type 53 and (xxxsel)"] -outfile analysis/tripeps_vmd/hbonds_interstrand_xtripepnumberx.dat -plot no
selection clear
foreach var [info vars] {
    unset $var
}

# Find Collagen-Water HBonds
hbonds -ang 30 -dist 3.0 -writefile yes -sel1 [atomselect top "type 349 or type 350"] -sel2 [atomselect top "xxxsel"] -outfile analysis/tripeps_vmd/hbonds_water_xtripepnumberx.dat -plot no
selection clear
foreach var [info vars] {
    unset $var
}

# GET WATERS WITHIN 3.3 A OF PROTEIN - WITH PBC
set dist 3.3
set num_frames 100000
set start_frame 0
set end_frame [expr {$num_frames - 1}]
# Open output file
set filename "analysis/tripeps_vmd/wats_prot_pb_3_3_xtripepnumberx.dat"
set outfile [open $filename w]
# Loop over each frame
for {set frame 0} {$frame < $num_frames} {incr frame} {
    if {[expr $frame % 1000] == 0} {
        puts "Processed $frame frames."
    }
    # Set the current frame
    animate goto $frame
    # Select O atoms within the given criteria
    set O_sel [atomselect top "type 349 and pbwithin $dist of (xxxsel) and type 231 50 1 5 53 234 806 809 820" frame $frame]
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

