# COFpiler

Build the statistical structure for layered materials by stacking layers following the Maxwell-Boltzmann energy distribution of their stacking modes.

## usage:
python COFpiler.py [--path PATH] [-data data_file] [--instr input_structure] [--mirror TRUE_or_FALSE] [--symmetry structure_symmetry] [--L layer_number] [--M model_number] [--tem synthesized temperature] [

## Optional arguments: 
--path PATH &emsp;The path where the infile and instr are, the output structure will also save there  
--data xlsx_file &emsp;The input file include data needed, form as: form as: the 1st column is the stacking_type, the 2nd is the Erel (relative energies), the 3rd, 4th and 5th columns are the x, y, z (vector of the shift_vectors)   
--instr input_structure &emsp; The input structure name as a start structure   
--mirror TRUE_or_FALSE &emsp; include the mirror of the shift or not   
--symmetry C4 &emsp; the symmetry for the structure, can be C3, C6 or C4   
--M model_number &emsp; the amount of the model one wanted to get  
