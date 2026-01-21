# TOTAL VMD SCRIPT

file mkdir analysis

# Find INTRAmolecular Hydrogen Bonds at 3.5 A - right
hbonds -ang 30 -dist 3.5 -writefile yes -sel2 [atomselect top "type 1 or type 50 or type 809"] -sel1 [atomselect top "type 5 or type 807 or type 53"] -outfile analysis/hbond_intra_3.5.dat -plot no
selection clear
foreach var [info vars] {
    unset $var
}

# FIND NH/CO HBonds with Water
hbonds -ang 30 -dist 3.5 -writefile yes -sel2 [atomselect top "type 1 or type 50 or type 809"] -sel1 [atomselect top "type 349 or type 350"] -outfile analysis/hbonds_water_3.5.dat -plot no
selection clear
foreach var [info vars] {
    unset $var
}

# Find INTRAmolecular Hydrogen Bonds at 3 A
hbonds -ang 30 -dist 3.0 -writefile yes -sel2 [atomselect top "type 1 or type 50 or type 809"] -sel1 [atomselect top "type 5 or type 807 or type 53"] -outfile analysis/hbond_intra_3.0.dat -plot no
selection clear
foreach var [info vars] {
    unset $var
}

# HBONDS BETWEEN GLY AND XAA PRO at 3 A - right
hbonds -ang 30 -dist 3.0 -writefile yes -sel1 [atomselect top "type 1"] -sel2 [atomselect top "type 53"] -outfile analysis/hbond_GLY_XAA_3.0.dat -plot no
selection clear
foreach var [info vars] {
    unset $var
}

# HBONDS BETWEEN GLY AND XAA PRO at 2.8 A
hbonds -ang 30 -dist 2.8 -writefile yes -sel1 [atomselect top "type 1"] -sel2 [atomselect top "type 53"] -outfile analysis/hbond_GLY_XAA_2.8.dat -plot no
selection clear
foreach var [info vars] {
    unset $var
}

# HBONDS BETWEEN C-alpha AND XAA PRO at 3.3 A - right
hbonds -ang 30 -dist 3.3 -writefile yes -sel1 [atomselect top "type 51"] -sel2 [atomselect top "type 53"] -outfile analysis/hbond_CA_XAA_3.3.dat -plot no
selection clear
foreach var [info vars] {
    unset $var
}

# HBONDS BETWEEN C-alpha AND XAA PRO at 3 A - right
hbonds -ang 30 -dist 3.0 -writefile yes -sel1 [atomselect top "type 51"] -sel2 [atomselect top "type 53"] -outfile analysis/hbond_CA_XAA_3.0.dat -plot no
selection clear
foreach var [info vars] {
    unset $var
}

# Find Glucose-water HBonds
hbonds -ang 30 -dist 3.0 -writefile yes -sel2 [atomselect top "type 601 or type 603 or type 605 or type 607 or type 609 or type 611"] -sel1 [atomselect top "type 349 or type 350"] -outfile analysis/hbond_gluc_water.dat -plot no
selection clear
foreach var [info vars] {
    unset $var
}

set resno 108
set prot 549
set num_frames 100000
set start_frame 0
set end_frame [expr {$num_frames - 1}]
# Dihedral Angles
file mkdir analysis/amoeba_rama
for {set i 0} {$i < $num_frames} {incr i 100} {
    animate goto $i
    set outfile [open "analysis/amoeba_rama/rama_plot_$i.dat" w]
    set sCA [atomselect top "type 2 or type 808 or type 51" frame $i]
    set sN [atomselect top "type 1 or type 231 or type 50 or type 809" frame $i]
    set sC [atomselect top "type 3 or type 801 or type 821 or type 52 or type 233" frame $i]
    set iCA [$sCA get index]
    set iN [$sN get index]
    set iC [$sC get index]
    # Start j from 1 as the first one cannot have a dihedral
    for {set j 0} {$j < $resno} {incr j} {
        set prev_C [lindex $iC [expr {$j - 1}]]
        set curr_N [lindex $iN [expr {$j}]]
        set curr_CA [lindex $iCA [expr {$j}]]
        set curr_C [lindex $iC [expr {$j}]]
        set next_N [lindex $iN [expr {$j + 1}]]
        set next_CA [lindex $iCA [expr {$j + 1}]]
        set next_C [lindex $iC [expr {$j + 1}]]

        if {[string is integer -strict $curr_C] && [string is integer -strict $next_N] && [string is integer -strict $next_CA] && [string is integer -strict $next_C]} {
        set omega1 [measure dihed [list $curr_C $next_N $next_CA $next_C] frame $i]} else {
        set omega1 0 }
        if {[string is integer -strict $curr_N] && [string is integer -strict $curr_CA] && [string is integer -strict $curr_C] && [string is integer -strict $prev_C]} {
        set phi1 [measure dihed [list $prev_C $curr_N $curr_CA $curr_C] frame $i] } else {
        set phi1 0 }
        if {[string is integer -strict $curr_N] && [string is integer -strict $curr_CA] && [string is integer -strict $curr_C] && [string is integer -strict $next_N]} {
        set psi1 [measure dihed [list $curr_N $curr_CA $curr_C $next_N] frame $i] } else {
        set psi1 0 }
        
        puts $outfile "$phi1 $psi1 $omega1" }
    close $outfile
}
selection clear
foreach var [info vars] {
    unset $var
}

set resno 108
set prot 549
set num_frames 100000
set start_frame 0
set end_frame [expr {$num_frames - 1}]
# Measure N-CA-CB-CG Dihedral for pyrrolidine ring pucker
file mkdir analysis/amoeba_pucker_dihed
for {set i 0} {$i < $num_frames} {incr i 100} {
    animate goto $i
    set outfile [open "analysis/amoeba_pucker_dihed/dihed_$i.dat" w]
    set sN [atomselect top "type 809 or type 50" frame $i]
    set sCA [atomselect top "type 808 or type 51" frame $i]
    set sCB [atomselect top "type 800 or type 55" frame $i]
    set sCG [atomselect top "type 804 or type 57" frame $i]
    set iN [$sN get index]
    set iCA [$sCA get index]
    set iCB [$sCB get index]
    set iCG [$sCG get index]
    set lN [llength $iN]
    set lCA [llength $iCA]
    set lCB [llength $iCB]
    set lCG [llength $iCG]
    for {set k 0} {$k < $lN} {incr k} {
        set kp [expr {$k + 1}]
        set currN [lindex $iN [expr {$k}]]
        set currCA [lindex $iCA [expr {$k}]]
        set currCB [lindex $iCB [expr {$k}]]
        set currCG [lindex $iCG [expr {$k}]]
        set dihedangle [measure dihed [list $currN $currCA $currCB $currCG] frame $i]
        puts $outfile "$kp $dihedangle"
    }
    close $outfile
}
selection clear
foreach var [info vars] {
    unset $var
}

set strand 36
set prot 549
set num_frames 100000
set start_frame 0
set end_frame [expr {$num_frames - 1}]
set sel [atomselect top "index 0 to [expr {$prot - 1}]"]
set ringo [atomselect top "type 601"]
set wats [atomselect top "type 349"]
set pg [measure gofr $sel $ringo delta 0.1 rmax 15.0 first $start_frame last $end_frame step 5 usepbc true]
set outfile [open "analysis/rdf_prot_ringo.dat" w]
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
set pw [measure gofr $sel $wats delta 0.1 rmax 15.0 first $start_frame last $end_frame step 5 usepbc true]
set outfile [open "analysis/rdf_prot_w.dat" w]
set r_values [lindex $pw 0]
set rdf_values [lindex $pw 1]
set int_values [lindex $pw 2]
for {set i 0} {$i < [llength $r_values]} {incr i} {
    set r [lindex $r_values $i]
    set rdf [lindex $rdf_values $i]
    set int [lindex $int_values $i]
    puts $outfile [format "%s %s %s" $r $rdf $int]
}
close $outfile
set ow [measure gofr $ringo $wats delta 0.1 rmax 15.0 first $start_frame last $end_frame step 5 usepbc true]
set outfile [open "analysis/rdf_ringo_w.dat" w]
set r_values [lindex $ow 0]
set rdf_values [lindex $ow 1]
set int_values [lindex $ow 2]
for {set i 0} {$i < [llength $r_values]} {incr i} {
    set r [lindex $r_values $i]
    set rdf [lindex $rdf_values $i]
    set int [lindex $int_values $i]
    puts $outfile [format "%s %s %s" $r $rdf $int]
}
close $outfile
set ww [measure gofr $wats $wats delta 0.1 rmax 10.0 first $start_frame last $end_frame step 10 usepbc true]
set outfile [open "analysis/rdf_w_w.dat" w]
set r_values [lindex $ww 0]
set rdf_values [lindex $ww 1]
set int_values [lindex $ww 2]
for {set i 0} {$i < [llength $r_values]} {incr i} {
    set r [lindex $r_values $i]
    set rdf [lindex $rdf_values $i]
    set int [lindex $int_values $i]
    puts $outfile [format "%s %s %s" $r $rdf $int]
}
close $outfile
selection clear
foreach var [info vars] {
    unset $var
}


set strand 36
set num_frames 100000
file mkdir analysis
# End to End Distance
set outfile [open "analysis/e2e.dat" w]
set start_frame 0
set end_frame [expr {$num_frames - 1}]
animate goto $start_frame
for {set i 0} {$i < 25000} {incr i} {
    animate goto $i
    set ca [atomselect top "type 2 or type 808 or type 51" frame $i]
    set ca_atoms [$ca get index]
    set length_ca [llength $ca_atoms]
    set ca10 [lindex $ca_atoms [expr {($strand * 0)}]]
    set ca11 [lindex $ca_atoms [expr {($strand * 1) - 1}]]
    set ca20 [lindex $ca_atoms [expr {($strand * 1)}]]
    set ca21 [lindex $ca_atoms [expr {($strand * 2) - 1}]]
    set ca30 [lindex $ca_atoms [expr {($strand * 2)}]]
    set ca31 [lindex $ca_atoms [expr {($strand * 3) - 1}]]
    set distA [measure bond [list $ca10 $ca11] frame $i]
    set distB [measure bond [list $ca20 $ca21] frame $i]
    set distC [measure bond [list $ca30 $ca31] frame $i]
    puts $outfile "$i $distA $distB $distC"
}
close $outfile
selection clear
foreach var [info vars] {
    unset $var
}
set strand 36
set num_frames 100000
set outfile [open "analysis/e2e.dat" a]
set start_frame 0
set end_frame [expr {$num_frames - 1}]
animate goto $start_frame
for {set i 25000} {$i < 50000} {incr i} {
    animate goto $i
    set ca [atomselect top "type 2 or type 808 or type 51" frame $i]
    set ca_atoms [$ca get index]
    set length_ca [llength $ca_atoms]
    set ca10 [lindex $ca_atoms [expr {($strand * 0)}]]
    set ca11 [lindex $ca_atoms [expr {($strand * 1) - 1}]]
    set ca20 [lindex $ca_atoms [expr {($strand * 1)}]]
    set ca21 [lindex $ca_atoms [expr {($strand * 2) - 1}]]
    set ca30 [lindex $ca_atoms [expr {($strand * 2)}]]
    set ca31 [lindex $ca_atoms [expr {($strand * 3) - 1}]]
    set distA [measure bond [list $ca10 $ca11] frame $i]
    set distB [measure bond [list $ca20 $ca21] frame $i]
    set distC [measure bond [list $ca30 $ca31] frame $i]
    puts $outfile "$i $distA $distB $distC"
}
close $outfile
selection clear
foreach var [info vars] {
    unset $var
}
set strand 36
set num_frames 100000
set outfile [open "analysis/e2e.dat" a]
set start_frame 0
set end_frame [expr {$num_frames - 1}]
animate goto $start_frame
for {set i 50000} {$i < 75000} {incr i} {
    animate goto $i
    set ca [atomselect top "type 2 or type 808 or type 51" frame $i]
    set ca_atoms [$ca get index]
    set length_ca [llength $ca_atoms]
    set ca10 [lindex $ca_atoms [expr {($strand * 0)}]]
    set ca11 [lindex $ca_atoms [expr {($strand * 1) - 1}]]
    set ca20 [lindex $ca_atoms [expr {($strand * 1)}]]
    set ca21 [lindex $ca_atoms [expr {($strand * 2) - 1}]]
    set ca30 [lindex $ca_atoms [expr {($strand * 2)}]]
    set ca31 [lindex $ca_atoms [expr {($strand * 3) - 1}]]
    set distA [measure bond [list $ca10 $ca11] frame $i]
    set distB [measure bond [list $ca20 $ca21] frame $i]
    set distC [measure bond [list $ca30 $ca31] frame $i]
    puts $outfile "$i $distA $distB $distC"
}
close $outfile
selection clear
foreach var [info vars] {
    unset $var
}
set strand 36
set num_frames 100000
set outfile [open "analysis/e2e.dat" a]
set start_frame 0
set end_frame [expr {$num_frames - 1}]
animate goto $start_frame
for {set i 75000} {$i < 100000} {incr i} {
    animate goto $i
    set ca [atomselect top "type 2 or type 808 or type 51" frame $i]
    set ca_atoms [$ca get index]
    set length_ca [llength $ca_atoms]
    set ca10 [lindex $ca_atoms [expr {($strand * 0)}]]
    set ca11 [lindex $ca_atoms [expr {($strand * 1) - 1}]]
    set ca20 [lindex $ca_atoms [expr {($strand * 1)}]]
    set ca21 [lindex $ca_atoms [expr {($strand * 2) - 1}]]
    set ca30 [lindex $ca_atoms [expr {($strand * 2)}]]
    set ca31 [lindex $ca_atoms [expr {($strand * 3) - 1}]]
    set distA [measure bond [list $ca10 $ca11] frame $i]
    set distB [measure bond [list $ca20 $ca21] frame $i]
    set distC [measure bond [list $ca30 $ca31] frame $i]
    puts $outfile "$i $distA $distB $distC"
}
close $outfile
selection clear
foreach var [info vars] {
    unset $var
}

# GET WATERS WITHIN 3.3 A OF PROTEIN
set strand 36
set dist 3.3
set prot 549
set num_frames 100000
set start_frame 0
set end_frame [expr {$num_frames - 1}]
# Open output file
set filename "wats_prot_[string map { . _ } $dist].dat"
set outfile [open $filename w]
# Loop over each frame
for {set frame 0} {$frame < $num_frames} {incr frame} {
    if {[expr $frame % 1000] == 0} {
        puts "Processed $frame frames."
    }
    # Set the current frame
    animate goto $frame
    # Select O atoms within the given criteria
    set O_sel [atomselect top "type 349 and within $dist of type 231 50 1 5 53 234 806 809 820" frame $frame]
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

set num_frames 100000
set start_frame 0
set end_frame [expr {$num_frames - 1}]
set filename "CA_coords.dat"
set outfile [open $filename w]
for {set frame 0} {$frame < $num_frames} {incr frame} {
    if {[expr $frame % 1000] == 0} {
        puts "Processed $frame frames."
    }
    animate goto $framec
    set sel [atomselect top "type 2 or type 51 or type 808" frame $frame]
    set indices [$sel get index]
    set types [$sel get type]
    set coords [$sel get {x y z}]
    set length [llength $indices]
    puts $outfile "$length"
    for {set i 0} {$i < $length} {incr i} {
        set idx [lindex $indices $i]
        set coord [lindex $coords $i]
        set type [lindex $types $i]
        set x [lindex $coord 0]
        set y [lindex $coord 1]
        set z [lindex $coord 2]
        puts $outfile "$idx $x $y $z $type"
    }
    $sel delete
}
close $outfile
puts "Finished writing coordinates to $filename"
selection clear
foreach var [info vars] {
    unset $var
}

# TOTAL VMD SCRIPT FOR N->PI* INTERACTIONS

set num_frames 100000

set filename "indices.dat"
set infile [open $filename r]
set lines [split [read $infile] "\n"]
close $infile

set col1 {}
set col2 {}
set col3 {}
set col4 {}

set data_lines [lrange $lines 1 end]

foreach line $data_lines {
    if {[string trim $line] eq ""} {
        continue  ;
    }
    set fields [split $line "\t"]
    lappend col1 [lindex $fields 0]
    lappend col2 [lindex $fields 1]
    lappend col3 [lindex $fields 2]
    lappend col4 [lindex $fields 3]
}

set outfile [open "NCACO_coords.dat" w]

animate goto 0
for {set i 0} {$i < $num_frames} {incr i 10} {
    animate goto $i
    # Selections for this frame
    # Shift all columns by subtracting 1 from each index
    set col1_shifted {}; foreach val $col1 { lappend col1_shifted [expr {$val - 1}] }
    set col2_shifted {}; foreach val $col2 { lappend col2_shifted [expr {$val - 1}] }
    set col3_shifted {}; foreach val $col3 { lappend col3_shifted [expr {$val - 1}] }
    set col4_shifted {}; foreach val $col4 { lappend col4_shifted [expr {$val - 1}] }

    # Create the atom selections with corrected indices
    set sN  [atomselect top "index [join $col1_shifted " "]" frame $i]
    set sCA [atomselect top "index [join $col2_shifted " "]" frame $i]
    set sC  [atomselect top "index [join $col3_shifted " "]" frame $i]
    set sO  [atomselect top "index [join $col4_shifted " "]" frame $i]

    # Get coordinate lists
    set n_coords  [$sN get {x y z}]
    set ca_coords [$sCA get {x y z}]
    set c_coords  [$sC get {x y z}]
    set o_coords  [$sO get {x y z}]
    # Number of atoms (assuming all lists are same length)
    set num_atoms [llength $n_coords]
    puts $outfile "$num_atoms"
    # Write coordinates
    for {set j 0} {$j < $num_atoms} {incr j} {
        set n  [lindex $n_coords $j]
        set ca [lindex $ca_coords $j]
        set c  [lindex $c_coords $j]
        set o  [lindex $o_coords $j]
        puts $outfile "[format {[%g,%g,%g]} {*}$n] [format {[%g,%g,%g]} {*}$ca] [format {[%g,%g,%g]} {*}$c] [format {[%g,%g,%g]} {*}$o]"
    }
    # Clean up
    $sN delete
    $sCA delete
    $sC delete
    $sO delete
}

close $outfile

selection clear
foreach var [info vars] {
    unset $var
}


set strand 36
set num_frames 100000
# BACKBONE RMSD
set filename [format "analysis/bb_rmsd.dat"]
set outfile [open $filename w]
set reference [atomselect top "type 51 or type 52 or type 50 or type 1 or type 2 or type 3 or type 231 or type 233 or type 808 or type 801 or type 821" frame 0]
for {set i 0} {$i < $num_frames} {incr i} {
    set compare [atomselect top "type 51 or type 52 or type 50 or type 1 or type 2 or type 3 or type 231 or type 233 or type 808 or type 801 or type 821" frame $i]
    set trans_mat [measure fit $compare $reference]
    $compare move $trans_mat
    set rmsd [measure rmsd $compare $reference]
    puts $outfile "$i $rmsd"
}
close $outfile
selection clear
foreach var [info vars] {
    unset $var
}

# GET WATERS WITHIN 3.8 A OF GLUCOSE
set strand 15
set dist 3.8
set prot 534
set total 14625
set num_frames 100000

set start_frame 0
set end_frame [expr {$num_frames - 1}]

# Open output file
set filename "wats_gluc_[string map { . _ } $dist].dat"
set outfile [open $filename w]

# Loop over each frame
for {set frame 0} {$frame < $num_frames} {incr frame} {
    if {[expr $frame % 1000] == 0} {
        puts "Processed $frame frames."
    }

    # Set the current frame
    animate goto $frame

    # Select O atoms within the given criteria
    set O_sel [atomselect top "type 349 and within $dist of type 601" frame $frame]

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


# GET GLUCOSES WITHIN 5 A OF PROTEIN
set strand 15
set dist 5
set prot 534
set total 14625e
set num_frames 100000

set filename "gluc_prot_[string map { . _ } $dist].dat"

set start_frame 0
set end_frame [expr {$num_frames - 1}]

# Open output file
set outfile [open $filename w]

# Loop over each frame
for {set frame 0} {$frame < $num_frames} {incr frame} {
    if {[expr $frame % 1000] == 0} {
        puts "Processed $frame frames."
    }

    # Set the current frame
    animate goto $frame

    # Select O atoms within the given criteria
    set sel [atomselect top "type 600 and within $dist of type 231 50 1 5 53 234" frame $frame]

    # Get indices and coordinates of selected O atoms
    set indices [$sel get index]
    set length [llength $indices]
    puts $outfile "$length"

    # Iterate over selected atoms
    set num_indices [llength $indices]
    for {set i 0} {$i < $num_indices} {incr i} {
        set idx [lindex $indices $i]

        # Print the current atom and the next 23 by index
        for {set j 0} {$j <= 23} {incr j} {
            set target_idx [expr $idx + $j]
            
            # Select atom by index (whether or not it was in the original selection)
            set target_sel [atomselect top "index $target_idx" frame $frame]

            # Ensure the selection is valid (in case the index is out of range)
            if {[$target_sel num] > 0} {
                set name [$target_sel get name]
                set type [$target_sel get type]
                set coord [$target_sel get {x y z}]
                set x [lindex $coord 0]
                set y [lindex $coord 1]
                set z [lindex $coord 2]

                puts $outfile "$target_idx $name $x $y $z $type"
            }
            $target_sel delete
        }
    }

    # Delete glucose selection
    $sel delete
}

# Close output file
close $outfile

puts "Finished writing coordinates to $filename"

selection clear
foreach var [info vars] {
    unset $var
}

quit