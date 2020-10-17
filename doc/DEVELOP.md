# Can I help develop/Add functionality to OTTER

Please contact the OTTERmaster (m.beck@uky.edu)

# Minimum Requirements

## Acquiring and building OTTER
To work with or on OTTER, please clone the entire package, including directories.  

To run OTTER, select the main branch, and simply type "make".  This will compile a current fully version of OTTER.  The binary 'otter.x' is the executable.

## Editing/developing code...
All code development must occur in a non-main branch. There is a development branch for each sub-version of OTTER, currently: fibers, spheres, ligaments.

After cloning, select the branch you would like to work on.  Then edit and work with files in the GIT managed directories, saving and committing locally as you see fit.  Please push your commits to the selected development branch frequently so that others can safely develop along with you.

## Directory structure notes
The makefile and directories are constructed to try to keep the repository clean and easy to work with.  Please put all .f90's in the src/ directory.  The obj/ directory is used during build only, and will be created if it doesn't exist.

NOTE: Currently all files to be built must be listed in order in the makefile.  Please number files you create to keep track of their build order for future automated build purposes.

## Questions/Comments
If you have any questions or concerns, please email m.beck@uky.edu