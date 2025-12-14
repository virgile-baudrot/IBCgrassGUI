myIP='163.172.138.51'
myIP='51.158.75.139'
ssh root@$myIP

####################### IMAGE SETUP SCRIPT ########################
# from remote
apt-get update && apt-get install -y g++ htop nano && rm -rf /var/lib/apt/lists/*
cd ../home/
mkdir IBCgrass
cd IBCgrass/
mkdir Input-EFSA

sudo apt update
sudo apt install r-base -y
sudo apt install build-essential libcurl4-openssl-dev libssl-dev libxml2-dev -y

sudo apt update
sudo apt install parallel -y


# IN R
# Rscript -e "install.packages('doParallel', repos='http://cran.r-project.org')"
R
> install.packages('foreach')
> install.packages('doParallel')
> install.packages('labeling')

# from local
scp -r Model-files/ root@$myIP:/home/IBCgrass/
scp -r Input-EFSA/ root@$myIP:/home/IBCgrass/

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

### DIRECT TEST RUNS #######################################################
./IBCgrassGUI 3 50 5 1 Fieldedge2.txt 10 90 100 0 0.1 0.01 1 1 1 0 0 6 1
./IBCgrassGUI 3 50 5 1 Fieldedge.txt 10 90 100 0 0.1 0.01 1 1 1 0 0 6 1
########## COPY SOME INPUT FILES FROM LOCAL ##############

########## COPY SOME FILES FROM REMOTE ##############
rsync -av --ignore-existing root@$myIP:/home/IBCgrass/runs2/* output-EFSA/

rsync -av --ignore-existing root@$myIP:/home/IBCgrass/Model-files/Pt* output-EFSA/
rsync -av --ignore-existing root@$myIP:/home/IBCgrass/Model-files/Grd* output-EFSA/
rsync -av --ignore-existing root@$myIP:/home/IBCgrass/sim* output-EFSA/
#######################################################

##################### FOR LOOP #########################
grBE2TH=("STE_I1" "BLS_I1" "ALP_I1" "ATL_I1" "BOR_I1" "CON_I1" "PAN_I1" "MED_I1" "MAC_I1")
STRESSORS=("NO" "biomass" "SEbiomass" "survival" "establishment" "sterility" "seednumber")


grBE2TH=("MAC_I1" "MED_I1" "PAN_I1" "CON_I1" "BOR_I1" "ATL_I1" "ALP_I1" "BLS_I1")
STRESSOR_TYPE=("low" "medium" "high")
STRESSORS=("biomass" "SEbiomass" "survival" "establishment" "sterility" "seednumber")
## A SINGLE SCRIPT EXAMPLE
# Rscript Input-EFSA/03_run_simulation.R "STE_I1" "low" "NO"
## A SINGLE SCRIPT IN BACKGROUND
# nohup Rscript Input-EFSA/03_run_simulation.R "STE_I1" "low" "biomass" > out.log 2>&1 &

### SIMPLE FOR LOOP #######################################################
for g in "${grBE2TH[@]}"; do
  for t in "${STRESSOR_TYPE[@]}"; do
    for s in "${STRESSORS[@]}"; do
      echo "Running: $g | $t | $s"
      Rscript Input-EFSA/03_run_simulation.R "$g" "$t" "$s"
    done
  done
done

### PARALLEL #######################################################
nproc # check number of cores

grBE2TH=("STE_I1" "BLS_I1" "ALP_I1" "ATL_I1" "BOR_I1" "CON_I1" "PAN_I1" "MED_I1" "MAC_I1")
STRESSOR_TYPE=("low" "medium" "high")
STRESSORS=("NO" "biomass" "SEbiomass" "survival" "establishment" "sterility" "seednumber")

##################### START PARALLEL LOOP #########################
nohup parallel -j 32 '
  dir="runs/{1}_{2}_{3}"

  mkdir -p "$dir"

  cp -r Input-EFSA "$dir"/
  cp -r Model-files "$dir"/Model-files

  cd "$dir"

  Rscript Input-EFSA/03_run_simulation.R {1} {2} {3}
' ::: \
  STE_I1 BLS_I1 ALP_I1 ATL_I1 BOR_I1 CON_I1 PAN_I1 MED_I1 MAC_I1 \
  # MAC_I1 MED_I1 PAN_I1 CON_I1 BOR_I1 ATL_I1 ALP_I1 BLS_I1 \
  ::: \
  low medium high \
  ::: \
  # biomass SEbiomass survival establishment sterility seednumber \
  NO biomass SEbiomass survival establishment sterility seednumber \
  > master.log 2>&1 &
##################### END PARALLEL LOOP #########################



ps aux | grep Rscript
kill <PID>