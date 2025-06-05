Instructions for Running the Program
There are two modes of operation depending on whether the user wants to estimate parameters or specify them manually.

1. Manual Parameter Specification
To run with fixed parameters (e.g., number of communities R = 3, ws = 0.7, wd = 0.2), edit
the top of the file where the Config struct is defined:
config.R = 3;
config.ws = 0.7;
config.wd = 0.2;

2. Automatic Parameter Estimation
To allow the program to estimate a parameter, set it to -1 in the Config struct:
config.R = -1; // Triggers R estimation
config.ws = -1.0; // Triggers ws estimation
config.wd = -1.0; // Triggers wd estimation
The program will then loop through candidate values for R (by default 2 to 6) and apply
EM to select the best model based on log-likelihood.

3. Additional Files
• Required Files: Place the following files in the devel/emdw/build directory:
– Adjacency matrix text file
– Ground truth communities file
– Visualization script (visualize_communities.py)

4. Configuring Text files
Change the values of the two string variables named: adjacencyMatrixFile and
groundTruthFile to the names of your adjacency matrix and ground truth text file.
These variables can be found just under the Config struct at the top of the code file.

• Visualization: After executing the main program, run the visualization script:
python visualise_communities.py
This will generate community visualisation.png showing the detected communities.
Ensure all files are in the same directory and that Python dependencies (NetworkX,
Matplotlib, Numpy) are installed before running the visualization script.

NOTE: the configuration for the parameters and text file names is near the top of the c++ file.
There are comments in the code that indicate where these respective edits should be made.