myIP='51.15.223.33'
myIP='163.172.141.36'
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
> install.packages(c('foreach','doParallel', 'labeling'))

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
##################### FOR LOOP #########################
grBE2TH=("STE_I1" "BLS_I1" "ALP_I1" "ATL_I1" "BOR_I1" "CON_I1" "PAN_I1" "MED_I1" "MAC_I1")
STRESSORS=("NO" "biomass" "SEbiomass" "survival" "establishment" "sterility" "seednumber")
## A SINGLE SCRIPT EXAMPLE
# Rscript Input-EFSA/03_run_simulation.R "STE_I1" "low" "NO"
## A SINGLE SCRIPT IN BACKGROUND
# nohup Rscript Input-EFSA/03_run_simulation.R "STE_I1" "low" "biomass" > out.log 2>&1 &

### PARALLEL #######################################################
nproc # check number of cores

##################### START PARALLEL LOOP #########################
# ------------
# CONTROL
# ------------
nohup parallel -j 9 '
  dir="runs_control_full/{1}_{2}_{3}"

  mkdir -p "$dir"

  cp -r Input-EFSA "$dir"/
  cp -r Model-files "$dir"/Model-files

  cd "$dir"

  Rscript Input-EFSA/03_run_simulation_TXT.R {1} {2} {3}
' ::: \
  ALP_I1 ATL_I1 BLS_I1 BOR_I1 CON_I1 MAC_I1 MED_I1 PAN_I1 STE_I1  \
  ::: \
  low \
  ::: \
  NO \
  > masterC.log 2>&1 &
  



nohup parallel -j 14 '
  dir="runs_effect_r20/{1}_{2}_{3}"

  mkdir -p "$dir"

  cp -r Input-EFSA "$dir"/
  cp -r Model-files "$dir"/Model-files

  cd "$dir"

  Rscript Input-EFSA/03_run_simulation_TXT_r20.R {1} {2} {3}
' ::: \
  ALP_I1 ATL_I1 BLS_I1 BOR_I1 CON_I1 MAC_I1 MED_I1 PAN_I1 STE_I1  \
  ::: \
  z01 z02 z03 z04 z05 \
  ::: \
  biomass SEbiomass survival establishment sterility seednumber \
  > masterE2.log 2>&1 &

# ------------
# EFFECT
# ------------
nohup parallel -j 14 '
  dir="runs_effect_r20/{1}_{2}_{3}"

  mkdir -p "$dir"

  cp -r Input-EFSA "$dir"/
  cp -r Model-files "$dir"/Model-files

  cd "$dir"

  Rscript Input-EFSA/03_run_simulation_TXT_r20.R {1} {2} {3}
' ::: \
  ALP_I1 ATL_I1 BLS_I1 BOR_I1 CON_I1 MAC_I1 MED_I1 PAN_I1 STE_I1  \
  ::: \
  z01 z02 z03 z04 z05 \
  ::: \
  biomass SEbiomass survival establishment sterility seednumber \
  > masterE2.log 2>&1 &
  
  
# ------------
# EFFECT for r20all
# ------------
nohup parallel --line-buffer -j 15 '
  dir="runs_effect_r20all/{1}_{2}"

  mkdir -p "$dir"

  cp -r Input-EFSA "$dir"/
  cp -r Model-files "$dir"/Model-files

  cd "$dir"

  Rscript Input-EFSA/03_run_simulation_TXT_r20all.R {1} {2}
' ::: \
  ALP_I1 ATL_I1 BLS_I1 BOR_I1 CON_I1 MAC_I1 MED_I1 PAN_I1 STE_I1  \
  ::: \
  z01 z02 z03 z04 z05 \
  > master.log 2>&1 &

##################### END PARALLEL LOOP #########################
find runs_effect_r20all -type f -name "Pt_*"
tail -f master.log

cp master.log backup_master.log


ps aux | grep Rscript
kill <PID>



########## COPY SOME FILES FROM REMOTE ##############
rsync -av --ignore-existing root@$myIP:/home/IBCgrass/* output-EFSA/
scp root@$myIP:/stockage/Pt__all_32.csv.gz .
scp root@$myIP:/stockage/Grd__all_32.csv.gz .

scp root@$myIP:/stockage/Grd__all_full.csv.gz .
scp root@$myIP:/stockage/Grd_summary.csv.gz .

scp -r output-EFSA_r20all/ root@$myIP:/home/IBCgrass/
scp output-EFSA_r20all/Pt__all_32.csv.gz root@$myIP:/stockage/
scp output-EFSA_r20all/* root@$myIP:/stockage/out/

rsync -av --ignore-existing root@$myIP:/home/IBCgrass/runs_control_r20_173 output-EFSA/
rsync -av --ignore-existing root@$myIP:/home/IBCgrass/runs_effect_r20all output-EFSA/
rsync -av --ignore-existing root@$myIP:/home/IBCgrass/*.zip output-EFSA/

rsync -av --ignore-existing root@$myIP:/home/IBCgrass/backup_master.log output-EFSA/
#######################################################

# ---------------------
# ATTENTION
# ---------------------


for d in */; do
    if [ -d "${d}Model-files" ]; then
        # 1. Déplacer les fichiers .txt
        echo "  -> Déplacement de ${d}Model-files/*.txt vers $d"
        mv "${d}Model-files/"*.txt "$d"
        mv "${d}Model-files/"*.csv "$d"
        
        # 2. Supprimer le dossier vide
        echo "  -> Suppression de ${d}Model-files"
        rm -rf "${d}Model-files"
    fi
done

for d in */; do
    echo "  -> Suppression de ${d}Input-EFSA"
    rm -rf "${d}Input-EFSA"
done