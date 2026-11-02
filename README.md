# IBCgrassGUI
Individual-based community model for grassland communities: version with herbicide impacts. Incl. graphical user interface

## About this fork (EFSA TASK 5-8 report)

This repository is a fork of [JetteReeg/IBCgrassGUI](https://github.com/JetteReeg/IBCgrassGUI)
(last upstream commit `e0dd994`). It was made by Virgile Baudrot (Qonfluens) for the EFSA
report "Review of IBC grass model for Risk Assessment" (TASK 5-8). It contains the code used to
run the simulations of the report, the R scripts that draw its figures, the figures, and the
simulation outputs they are drawn from. The model itself is the one of Jette Reeg; the changes
listed below are small and are needed to run it in batch on Linux and to pair the runs.

### What was changed in the model (`Model-files/`)

| File | Change | Effect on the results |
|---|---|---|
| `GMHerbicideEffect.cpp` | Optional 19th argument: seed of the random number generator. It is used for `initLCG` and `srand`. Without it, the seed is the process ID and the time, as before. | A control run and its treatment runs can use the same seed (paired runs). |
| `GMHerbicideEffect.cpp`, `CGridEnvir.cpp`, `RunPara.*` | Optional 20th argument: a second seed, set again at the start of year Tinit+1. | Gives a "null run": a control that only differs from its paired control by the random stream after Tinit, as a treatment run does. It measures the effect of stochasticity alone. |
| `CTKmodel.cpp` | The application rate is read as a `double` instead of an `int`. | Before, rates below 1 were truncated to 0 and other rates were rounded down. Results with integer rates are unchanged. |
| `SPftTraits.cpp` | The PFT file reader stops with an error when a mandatory trait is missing or is not a number (e.g. `NA`). | Before, such values were read as 0, and `memory = 0` made the run loop without end. `overwintering` and `comp_power` stay optional. |
| `CGridEnvir.cpp` | The per-plant output buffer (`IndOutData`) is no longer filled. | It was never written to a file. Removing it saves memory and time. No output changes. |
| `Makefile` (new) | Builds the static command-line binary on Linux: `make -C Model-files`. | None. |

### What was added

- `Input-EFSA/`: scripts of the report.
  - `00_*.R`, `01_create_BE2TH_FieldEdge.R`: species input built from the BE2TH trait data base (life span mapped to the trait `MaxAge`, species with a missing trait removed).
  - `02_simulation_setting.R`: settings and helpers shared by the run scripts (species sensitivity distribution in txt mode, seeds, call of the model).
  - `run_4.1_single_species.R`, `run_4.3_control_design.R`, `run_4.4_stressor_effects.R`: one run script per section of Chapter 4. Each run has its own directory, its seed and a `provenance*.json` file.
  - `04_output_analysis.R` to `18_yearly_mean_endpoints.R` (except `11`), `BE2TH_explore.R`: figures and tables of the report (§1, §4.2 to §4.4). The section and the figure are given at the top of each script.
  - `scaleway/`: Docker image and scripts used to run §4.3 and §4.4 on Scaleway Serverless Jobs.
  - `archive/`: earlier scripts, kept for the record. They are not used for the figures of the report.
- `output-SINGLE/analyze_outputs_100.R`, `analyze_outputs_300.R`: figures of §4.1.
- `EFSA_REPORT/figures/`: the figures and tables of the report, one folder per section.
- `outputs/`: the simulation outputs read by the figure scripts (stored with Git LFS, see below).
- `Dockerfile`, `Dockerfile.vnc`, `start.sh`, `run_control.R`, `run_herbicide.R`: run the model without the GUI, or the GUI in a browser (see the end of this file).

### What was not changed

- The ecological processes of IBC-grass (growth, competition, seed production, dispersal, germination, mortality, grazing, trampling, cutting) and the herbicide effect model (dose-response and txt modes).
- The default parameters (`Model-files/RunParams.txt`, `Input-files/`) and the PFT definitions of the GUI.
- The graphical user interface (`R-files/`). Only the order of the loaded packages was changed, for Linux.
- When the two new optional arguments are not given, the model runs as the upstream version.

The bundled Windows R distribution (`R/`) and the R libraries (`Rlibraries/`) of the upstream
repository were removed. R and the packages are installed by the Dockerfiles instead.

### Simulation outputs (Git LFS)

The outputs are stored with [Git LFS](https://git-lfs.com) (about 2.3 GB). Install Git LFS
before cloning, or run `git lfs pull` after the clone.

| Folder | Content | Read by |
|---|---|---|
| `outputs/4.1_single_species/` | §4.1, single species (stressors and variability) | `output-SINGLE/analyze_outputs_*.R` |
| `outputs/control-2026-10-28/` | §4.2 and §4.3.1, control runs (9 biomes, 4 initial sizes, 5 replicates) | `10`, `12`, `13`, `15` |
| `outputs/grid-2026-10-30/` | §4.3.2, grid sizes 10 x 10 to 350 x 350 (Boreal pool) | `14`, `15` |
| `outputs/main-2026-10-27/` | §4.4, stressor effects (9 biomes, 10 communities per biome, SSD levels z01 to z09) | `04` to `10`, `16` to `18` |

Each run directory holds the `Grd.csv.gz` and `Pt.csv.gz` outputs of the model, `runs.csv` and a
`provenance*.json` file (date, commit, seed and settings of the run). The figures are drawn
from the repository root, for example:

```bash
Rscript Input-EFSA/16_report_4.4_figures.R
```

To run the simulations again, build the model (`make -C Model-files`) and run the
`Input-EFSA/run_*.R` script of the section. The full §4.4 design was run on Scaleway
(`Input-EFSA/scaleway/README.md`).

---

*The original README of Jette Reeg follows.*

## Author
Jette Reeg
### Authors of earlier versions of IBCgrass
Felix May, Ines Steinhauer, Katrin Koerner, Lina Weiss, Hans Pfestorf
## Short discription
The IBC-grass GUI was developed to facilitate the use of the plant community model IBC-grass for herbicide risk assessments of non-target terrestrial plant communities in Central Europe.  Users are able to run simulations for various plant communities, which may differ in their species composition, environmental parameter settings and herbicide settings. Several different outputs can be generated on population as well as on community level.
Detailed information on the model can be found in the ODD-protocol and GMP document. 
For support, please contact Jette Reeg (jreeg@uni-potsdam.de).

Licence  
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

## Language
The IBCgrass model is written in C++(11) and needs to be compiled with g++.
The Graphical User Interface (GUI) is written in R using the packages RGtk2, RGtk2Extras, data.table, reshape2, foreach, doParallel, ggplot2 and ggthemes.

## Requirements

To run the IBCgrass GUI, following software need to be installed on the local machines:

GTK		The R package RGtk2 depends on gtk. It might be necessary to install gtk by hand. 

C++ compiler	e.g. MinGW (http://www.mingw.org/wiki/Getting_Started)

The software was tested on Windows 7, Windows 10 with following version:

R		3.4.1 and 3.5.2    
g++		4.8.1    
RGtk2		2.20.35 (with gtk 2.22.1)    
RGtk2Extras	0.6.1 (not supported anymore under R Version >=3.6)    
data.table	1.12.0   
reshape2	1.4.3   
foreach		1.4.4   
doParallel	1.0.14   
ggplot2		3.1.0   
ggthemes	4.0.1   

The software might be used on Ubuntu OS, if R including all necessary
packages are installed locally.
In that case the file "RunIBCgrassGUI.bat" needs to be adapted.
(Please contact jreeg@uni-potsdam.de for further details.)



# Using the Docker Image to Run IBC-grass without

Build the image:
```bash
docker build -t ibcgrass-runner .
```

Run a container:
```bash
docker run -it --rm ibcgrass-runner
```
To have after running from
```bash
docker ps
docker exec -it <CID> bash
docker exec -it 2f82c469abf bash
```

Inside the container, run a simulation:
```bash
Rscript run_control.R
# or
Rscript run_herbicide.R
```

Output files will be located in /home/rstudio/IBCgrass/Model-files/
To get files out of the container, use `docker cp` or mount a volume:
```bash
docker run -it --rm -v $(pwd)/output:/home/rstudio/IBCgrass/Model-files ibcgrass-runner
```
This will map the container's output directory to a local 'output' directory.

# Using the Docker Image to Run IBC-grass with GUI

```bash
docker build -t ibcgrassgui:vnc -f Dockerfile.vnc .
docker run --rm -it -p 6080:6080 ibcgrassgui:vnc
```

Then open http://localhost:6080 in the browser.

You can access directly to the software by using:
```bash
http://localhost:6080/vnc.html?autoconnect=1&host=localhost&port=6080&path=websockify
```

If required, we can access the image using:
```bash
docker run --rm -it --entrypoint bash ibcgrassgui:vnc
```

```bash
docker ps # and record the container ID
CID="d4214b172be7"
# remplace <CID> par l’ID/nom du conteneur déjà en marche
docker exec -it -e DISPLAY=:1 -e XDG_RUNTIME_DIR=/tmp/runtime-root $CID R
# Or to go in the bash
docker exec -it -e DISPLAY=:1 -e XDG_RUNTIME_DIR=/tmp/runtime-root $CID bash
```