# Running the EFSA simulations on Scaleway Serverless Jobs

Resources (project `simulation`, region fr-par):

| Resource | Name / ID |
|---|---|
| Registry namespace (private) | `rg.fr-par.scw.cloud/ibcgrass` |
| Image | `rg.fr-par.scw.cloud/ibcgrass/ibcgrass-efsa:<commit>` |
| Job definition | `ibcgrass-efsa`, `b7b866d3-b307-4011-816c-9acfdeacd3af` (6 vCPU, 12 GB, timeout 12 h) |
| Bucket | `ibcgrass-efsa-runs` |

Serverless Jobs allow at most 6 vCPU per job run. One job run = one community of §4.4
(32 model runs on 6 cores: control, null run and 30 treatments; 31 runs took 40 to 97 min in
the tests of 26 Oct 2026, i.e. 8 to 16 min per run).

## 1. Build and push the image (after each code change, from the repository root)

```bash
TAG=$(git rev-parse --short HEAD)          # commit first: the tag and IBC_COMMIT name the code
docker build -f Input-EFSA/scaleway/Dockerfile --build-arg IBC_COMMIT=$TAG -t ibcgrass-efsa .
scw registry login
docker tag ibcgrass-efsa rg.fr-par.scw.cloud/ibcgrass/ibcgrass-efsa:$TAG
docker push rg.fr-par.scw.cloud/ibcgrass/ibcgrass-efsa:$TAG
scw jobs definition update b7b866d3-b307-4011-816c-9acfdeacd3af image-uri=rg.fr-par.scw.cloud/ibcgrass/ibcgrass-efsa:$TAG
```

Local test of the image (small grid, outputs in ./tmp_out):

```bash
docker run --rm --user $(id -u) -v $PWD/tmp_out:/tmp/out -e ARGS="BOR_I1 1" -e NCORES=8 \
  -e TEST_CELLNUM=50 -e TEST_TMAX=6 -e TEST_TINIT=2 -e TEST_HERBDURATION=2 ibcgrass-efsa
```

## 2. Object Storage keys for the jobs

The job copies its outputs to the bucket with `rclone`, which reads `AWS_ACCESS_KEY_ID` and
`AWS_SECRET_ACCESS_KEY`. Store them in Secret Manager and attach them to the job definition:

```bash
P=989efe92-d59f-49ea-b987-f031cebf7818   # project "simulation"
JOB=b7b866d3-b307-4011-816c-9acfdeacd3af
AK=$(scw secret secret create name=ibcgrass-s3-access-key project-id=$P -o json | jq -r .id)
SK=$(scw secret secret create name=ibcgrass-s3-secret-key project-id=$P -o json | jq -r .id)
scw secret version create secret-id=$AK data="$(scw config get access-key)"
scw secret version create secret-id=$SK data="$(scw config get secret-key)"
scw jobs secret create job-definition-id=$JOB \
  secrets.0.secret-manager-id=$AK secrets.0.secret-manager-version=latest secrets.0.env-var-name=AWS_ACCESS_KEY_ID \
  secrets.1.secret-manager-id=$SK secrets.1.secret-manager-version=latest secrets.1.env-var-name=AWS_SECRET_ACCESS_KEY
```

A dedicated IAM application limited to Object Storage in the `simulation` project is safer
than your personal keys, if you prefer.

## 3. Run, follow, fetch, check

```bash
bash Input-EFSA/scaleway/submit.sh dryrun 2 BOR_I1          # S4.1 dry run: 2 communities
bash Input-EFSA/scaleway/submit.sh main-2026-10-30 10        # S4.2: 9 biomes x 10 MC = 90 job runs
scw jobs run list -o human                                   # states
bash Input-EFSA/scaleway/fetch.sh main-2026-10-30            # -> outputs/main-2026-10-30/
Rscript Input-EFSA/run_4.4_stressor_effects.R check outputs/main-2026-10-30/4.4_stressor_effects 10
```

A failed community is started again with the same biome and MC (same seed, same result if
the failure was not a model error). §4.1 and §4.3 can use the same job definition, e.g.
`scw jobs definition start $JOB environment-variables.SCRIPT=run_4.3_control_design.R
environment-variables.ARGS="BOR_I1 all" environment-variables.BUCKET=ibcgrass-efsa-runs
environment-variables.PREFIX=<prefix>/4.3_control_design`.

Grid-size series of §4.3.2 (Boreal, full pool, trampling off, 10 x 10 to 700 x 700 cells, one
job run per size; 5 runs in parallel, about 50 min at 250 x 250 and 8 h at 700 x 700):

```bash
BIOMES="BOR_I1" bash Input-EFSA/scaleway/submit_4.3.sh grid-2026-10-30 g010 g025 g050 g100 g173
BIOMES="BOR_I1" bash Input-EFSA/scaleway/submit_4.3.sh grid-2026-10-30 g250 g350 g500 g700
bash Input-EFSA/scaleway/fetch.sh grid-2026-10-30        # -> outputs/grid-2026-10-30/4.3_control_design/BOR_I1/g*/
```
