# IBCgrass Simulation Guide: Single Species & Effect Test

This folder (`Input-SINGLE`) contains all the scripts necessary to run a targeted 50% effect simulation on a single species to observe model outputs. 

## Folder Structure
This setup assumes you are working on a remote machine (e.g., Scaleway) and your directory structure looks like this:
```text
/home/IBCgrass/
├── Model-files/         # Contains the C++ source code and compiled executable
└── Input-SINGLE/        # (THIS FOLDER) Contains BE2TH_cleaned.csv and R scripts
```

## Step-by-Step Instructions

### 1. Connect and Compile the Model (Remote)
SSH into your remote machine and ensure the C++ model is compiled inside `Model-files/`:
```bash
ssh root@your_scaleway_ip

# Compile the executable (Run this from /home/IBCgrass)
cd /home/IBCgrass
g++ -static -static-libgcc -static-libstdc++ -std=c++11 -O2 -c Model-files/*.cpp
g++ -static -static-libgcc -static-libstdc++ -std=c++11 -O2 -o Model-files/IBCgrassGUI Model-files/*.o
```

### 2. Run the Scenarios (Remote)
Navigate to this folder and execute the scenarios. The script handles the
dynamic creation of `Fieldedge.txt` and `HerbFact.txt` for each 50% effect modifier.
```bash
cd /home/IBCgrass/Input-SINGLE
Rscript run_scenarios.R
```

### 3. Analyze and Plot the Outputs (Remote)
Once the simulations are finished, run the analysis script to aggregate
the generated `output_*.csv` files and create the faceted plot.
```bash
Rscript analyze_outputs.R
```

### 4. Retrieve the Results (Local)
Open a terminal on your **local computer** (not the SSH session) to securely 
sync the resulting graphics and CSVs down to your local machine:
```bash
# Syncs the contents to a local folder named 'output-SINGLE'
rsync -avzP root@your_scaleway_ip:/home/IBCgrass/Input-SINGLE/ ./output-SINGLE/
```