# Load the LAMMPS trajectory file
mol new dump.lammpstrj type lammpstrj waitfor all

# Select atom types (Modify this based on your system)
set myAtoms [atomselect top "type 7 or type 9"]  ;# Change types if needed

# Open output file
set outfile [open "sasa_per_frame.txt" w]

# Get total number of frames
set num_frames [molinfo top get numframes]

# Loop through all frames
for {set i 0} {$i < $num_frames} {incr i} {
    animate goto $i  ;# Go to frame i
    set sasa_value [measure sasa 1.4 $myAtoms]  ;# Compute SASA
    puts $outfile "$i $sasa_value"  ;# Save frame number and SASA value
}

# Close file and cleanup
close $outfile
$myAtoms delete

