# Use rocker/verse as a base image, which comes with R and many useful packages
FROM rocker/verse:4.1.0

# Switch to root user to install system dependencies and set file permissions
USER root

# Install g++ compiler
RUN apt-get update && apt-get install -y g++ htop nano && rm -rf /var/lib/apt/lists/*

# Set the working directory
WORKDIR /opt/IBCgrass

# Copy files with the correct ownership for the non-root user 'rstudio'
# This is the key fix for the "Permission denied" error.
COPY --chown=rstudio:rstudio Model-files/ ./Model-files/
COPY --chown=rstudio:rstudio Input-files/ ./Input-files/
COPY --chown=rstudio:rstudio run_control.R .
COPY --chown=rstudio:rstudio run_herbicide.R .
COPY --chown=rstudio:rstudio ExampleAnalyses/ ./ExampleAnalyses/

# Now, switch to the non-root user for all subsequent operations
# USER rstudio

# Compile the C++ code as the 'rstudio' user
# This user now has write permissions in the target directory
RUN g++ -static -static-libgcc -static-libstdc++ -std=c++11 -O2 -c Model-files/CEnvir.cpp -o Model-files/CEnvir.o && \
    g++ -static -static-libgcc -static-libstdc++ -std=c++11 -O2 -c Model-files/CGrid.cpp -o Model-files/CGrid.o && \
    g++ -static -static-libgcc -static-libstdc++ -std=c++11 -O2 -c Model-files/SPftTraits.cpp -o Model-files/SPftTraits.o && \
    g++ -static -static-libgcc -static-libstdc++ -std=c++11 -O2 -c Model-files/OutStructs.cpp -o Model-files/OutStructs.o && \
    g++ -static -static-libgcc -static-libstdc++ -std=c++11 -O2 -c Model-files/CSeed.cpp -o Model-files/CSeed.o && \
    g++ -static -static-libgcc -static-libstdc++ -std=c++11 -O2 -c Model-files/LCG.cpp -o Model-files/LCG.o && \
    g++ -static -static-libgcc -static-libstdc++ -std=c++11 -O2 -c Model-files/CTDSeed.cpp -o Model-files/CTDSeed.o && \
    g++ -static -static-libgcc -static-libstdc++ -std=c++11 -O2 -c Model-files/CGenet.cpp -o Model-files/CGenet.o && \
    g++ -static -static-libgcc -static-libstdc++ -std=c++11 -O2 -c Model-files/CObject.cpp -o Model-files/CObject.o && \
    g++ -static -static-libgcc -static-libstdc++ -std=c++11 -O2 -c Model-files/Cell.cpp -o Model-files/Cell.o && \
    g++ -static -static-libgcc -static-libstdc++ -std=c++11 -O2 -c Model-files/CTDPlant.cpp -o Model-files/CTDPlant.o && \
    g++ -static -static-libgcc -static-libstdc++ -std=c++11 -O2 -c Model-files/GMHerbicideEffect.cpp -o Model-files/GMHerbicideEffect.o && \
    g++ -static -static-libgcc -static-libstdc++ -std=c++11 -O2 -c Model-files/Plant.cpp -o Model-files/Plant.o && \
    g++ -static -static-libgcc -static-libstdc++ -std=c++11 -O2 -c Model-files/CTKmodel.cpp -o Model-files/CTKmodel.o && \
    g++ -static -static-libgcc -static-libstdc++ -std=c++11 -O2 -c Model-files/CGridEnvir.cpp -o Model-files/CGridEnvir.o && \
    g++ -static -static-libgcc -static-libstdc++ -std=c++11 -O2 -c Model-files/RunPara.cpp -o Model-files/RunPara.o && \
    g++ -static -static-libgcc -static-libstdc++ -std=c++11 -O2 -c Model-files/CHerbEff.cpp -o Model-files/CHerbEff.o && \
    g++ -static -static-libgcc -static-libstdc++ -o Model-files/IBCgrassGUI Model-files/SPftTraits.o Model-files/RunPara.o Model-files/Plant.o Model-files/OutStructs.o Model-files/LCG.o Model-files/GMHerbicideEffect.o Model-files/Cell.o Model-files/CTKmodel.o Model-files/CTDSeed.o Model-files/CTDPlant.o Model-files/CSeed.o Model-files/CObject.o Model-files/CHerbEff.o Model-files/CGridEnvir.o Model-files/CGrid.o Model-files/CGenet.o Model-files/CEnvir.o

# Install required R packages
RUN R -e "install.packages(c('foreach', 'doParallel', 'labeling'), repos='http://cran.rstudio.com/')"

# Set default command to bash
CMD ["/bin/bash"]

# --- How to use this Dockerfile ---
# 1. Build the image:
#    docker build -t ibcgrass-runner .
#
# 2. Run a container:
#    docker run -it --rm ibcgrass-runner
#
# 3. Inside the container, run a simulation:
#    Rscript run_control.R
#    or
#    Rscript run_herbicide.R
#
# Output files will be located in /home/rstudio/IBCgrass/Model-files/
# To get files out of the container, use `docker cp` or mount a volume:
# docker run -it --rm -v $(pwd)/output:/home/rstudio/IBCgrass/Model-files ibcgrass-runner
# This will map the container's output directory to a local 'output' directory.
