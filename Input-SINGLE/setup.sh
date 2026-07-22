myIP='51.15.223.33'
ssh root@$myIP

####################### IMAGE SETUP SCRIPT ########################
# from remote
apt-get update && apt-get install -y g++ htop nano && rm -rf /var/lib/apt/lists/*
cd ../home/
mkdir IBCgrass
cd IBCgrass/
mkdir Input-SINGLE
mkdir output-SINGLE

sudo apt update
sudo apt install r-base -y
sudo apt install build-essential libcurl4-openssl-dev libssl-dev libxml2-dev -y

sudo apt update
sudo apt install parallel -y

# IN R
# Rscript -e "install.packages('doParallel', repos='http://cran.r-project.org')"
install.packages(c('foreach', 'doParallel', 'labeling', 'dplyr', 'tidyr', 'ggplot2'))

# from local
# scp -r Model-files/ root@$myIP:/home/IBCgrass/
rsync -avzP --delete Input-SINGLE/ root@$myIP:/home/IBCgrass/Input-SINGLE/

# from remote
g++ -static -static-libgcc -static-libstdc++ -std=c++11 -O2 -c Model-files/CEnvir.cpp -o Model-files/CEnvir.o
g++ -static -static-libgcc -static-libstdc++ -std=c++11 -O2 -c Model-files/CGrid.cpp -o Model-files/CGrid.o
g++ -static -static-libgcc -static-libstdc++ -std=c++11 -O2 -c Model-files/SPftTraits.cpp -o Model-files/SPftTraits.o
g++ -static -static-libgcc -static-libstdc++ -std=c++11 -O2 -c Model-files/OutStructs.cpp -o Model-files/OutStructs.o
g++ -static -static-libgcc -static-libstdc++ -std=c++11 -O2 -c Model-files/CSeed.cpp -o Model-files/CSeed.o
g++ -static -static-libgcc -static-libstdc++ -std=c++11 -O2 -c Model-files/LCG.cpp -o Model-files/LCG.o
g++ -static -static-libgcc -static-libstdc++ -std=c++11 -O2 -c Model-files/CTDSeed.cpp -o Model-files/CTDSeed.o
g++ -static -static-libgcc -static-libstdc++ -std=c++11 -O2 -c Model-files/CGenet.cpp -o Model-files/CGenet.o
g++ -static -static-libgcc -static-libstdc++ -std=c++11 -O2 -c Model-files/CObject.cpp -o Model-files/CObject.o
g++ -static -static-libgcc -static-libstdc++ -std=c++11 -O2 -c Model-files/Cell.cpp -o Model-files/Cell.o
g++ -static -static-libgcc -static-libstdc++ -std=c++11 -O2 -c Model-files/CTDPlant.cpp -o Model-files/CTDPlant.o
g++ -static -static-libgcc -static-libstdc++ -std=c++11 -O2 -c Model-files/GMHerbicideEffect.cpp -o Model-files/GMHerbicideEffect.o
g++ -static -static-libgcc -static-libstdc++ -std=c++11 -O2 -c Model-files/Plant.cpp -o Model-files/Plant.o
g++ -static -static-libgcc -static-libstdc++ -std=c++11 -O2 -c Model-files/CTKmodel.cpp -o Model-files/CTKmodel.o
g++ -static -static-libgcc -static-libstdc++ -std=c++11 -O2 -c Model-files/CGridEnvir.cpp -o Model-files/CGridEnvir.o
g++ -static -static-libgcc -static-libstdc++ -std=c++11 -O2 -c Model-files/RunPara.cpp -o Model-files/RunPara.o
g++ -static -static-libgcc -static-libstdc++ -std=c++11 -O2 -c Model-files/CHerbEff.cpp -o Model-files/CHerbEff.o

g++ -static -static-libgcc -static-libstdc++ -std=c++11 -O2 \
  -o Model-files/IBCgrassGUI \
  Model-files/SPftTraits.o \
  Model-files/RunPara.o \
  Model-files/Plant.o \
  Model-files/OutStructs.o \
  Model-files/LCG.o \
  Model-files/GMHerbicideEffect.o \
  Model-files/Cell.o \
  Model-files/CTKmodel.o \
  Model-files/CTDSeed.o \
  Model-files/CTDPlant.o \
  Model-files/CSeed.o \
  Model-files/CObject.o \
  Model-files/CHerbEff.o \
  Model-files/CGridEnvir.o \
  Model-files/CGrid.o \
  Model-files/CGenet.o \
  Model-files/CEnvir.o

###################################################################

# FROM REMOTE
# Create the working folder and navigate into it
cd ~/Input-SINGLE

## THEN AFTER
rm -f ../Model-files/Grd_* ../Model-files/Pt_*
Rscript run_scenarios.R

Rscript run_scenarios.R 100 900
Rscript extract_outputs.R 100

Rscript run_scenarios.R 300 173
Rscript extract_outputs.R 300

# Download just the CSVs and the graphic
rsync -avzP root@$myIP:/home/IBCgrass/output-SINGLE/ ./output-SINGLE/
