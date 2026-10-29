# Image for the EFSA simulations (run_4.1 / run_4.3 / run_4.4) on Scaleway Serverless Jobs.
# Build from the repository root:
#   docker build -f Input-EFSA/scaleway/Dockerfile --build-arg IBC_COMMIT=$(git rev-parse --short HEAD) -t ibcgrass-efsa .
FROM rocker/r-ver:4.5.2

RUN apt-get update && apt-get install -y --no-install-recommends g++ make rclone ca-certificates \
    && rm -rf /var/lib/apt/lists/* \
    && install2.r --error data.table jsonlite

WORKDIR /opt/IBCgrass
COPY Model-files/*.cpp Model-files/*.h Model-files/Makefile Model-files/
RUN make -C Model-files clean && make -C Model-files -j"$(nproc)" && rm -f Model-files/*.o
COPY Input-EFSA/02_simulation_setting.R Input-EFSA/run_4.1_single_species.R \
     Input-EFSA/run_4.3_control_design.R Input-EFSA/run_4.4_stressor_effects.R \
     Input-EFSA/BE2TH_cleaned.csv Input-EFSA/
COPY Input-EFSA/scaleway/entrypoint.sh /usr/local/bin/entrypoint.sh

ARG IBC_COMMIT=unknown
ENV IBC_ROOT=/opt/IBCgrass IBC_COMMIT=${IBC_COMMIT}
ENTRYPOINT ["/bin/sh", "/usr/local/bin/entrypoint.sh"]
