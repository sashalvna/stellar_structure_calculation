Stellar structure calculation for a Zero Age Main-Sequence (ZAMS) Sun-like star for the Stellar Structure and Evolution final project. 

To run, simply run main.py, which constructs a model of the structure of a ZAMS star with the parameters defined at the bottom of the script. 

To change any parameters to their desired values in main.py under 'Set up inital parameters'. Guesses for the boundary parameters can be changed under 'Inner boundary guess', and 'Outer boundary guess'. If changing the composition of the star, find the matching opacity table: https://cds.unistra.fr/topbase/OpacityTables.html and name it according to the naming convention used for the solar composition table provided.

For the comparison of model results with MESA, the mesa_reader module must be installed. If changing the mass or composition of the star, run MESA to get a ZAMS model with the same parameters and provide the profile.data file to main.py by changing the mesapath and mesafile parameters.
