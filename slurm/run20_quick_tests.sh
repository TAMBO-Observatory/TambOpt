#!/bin/bash
#SBATCH -p gpu_test
#SBATCH --gres=gpu:4
#SBATCH -c 16
#SBATCH --mem=200g
#SBATCH -t 05:00:00
#SBATCH -o slurm_logs/slurm-%j-%x.out
#SBATCH --chdir=/n/home05/zdimitrov/tambo/TambOpt
# Run 20 quick tests: SGD from grid_n99 (seed-0 grid of 99, as the validation baseline), lr 1 m, 100 steps, scored every
# 25; 1 deg Poisson utility, default recon, expected counts, no background. 4 MIG slices:
#   energy     log10 E [GeV] bands 5.00-5.67, 5.67-6.33, 6.33-7.00 and the full range; 8 x 2,048 showers loaded per step
#              (a band keeps ~1/3), scoring set 3,072 (a band keeps ~1,000)
#   construct  one construction penalty each (modules/optimize/construction.py), full range, 4 x 2,048 per step; lambda
#              chosen so the penalty's gradient at the grid is ~half the utility's per-detector RMS (1e-3 of 2e-3):
#              slope 0.03, access 3, pads 0.2, contour 0.1
# SMOKE=1: 3 steps of every config of the set, one after another on one slice. ONLY=<0-3>: only that config, on one GPU
# (one job per config on gpu_requeue: sbatch -p gpu_requeue --gres=gpu:1 ...; resumes from opt_state.pt after a requeue).
# Out-suffix quick_<name>; logs <run>/slurm_logs/quick_<name>.log.
#   sbatch -J r20_quick_energy slurm/run20_quick_tests.sh energy
source slurm/env.sh
source slurm/run20_guard.sh
export MPLCONFIGDIR=/tmp/mplcache_$USER PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
R=$(python -c 'import modules.constants as C; print(C.RUN_LOCATION)')
L=$R/slurm_logs
GRID=$R/diagnostics/triplets/layouts/grid_n99.pt
export TAMBOPT_UTILITY=swgo_trigger TAMBOPT_TWO_TIER_LAMBDA=0.3 TAMBOPT_TRIG_K=3 TAMBOPT_TRIG_AREA=2 \
       TAMBOPT_TRIG_MU0=1.0 TAMBOPT_TRIG_SLOPE=10.0 TAMBOPT_TRIG_MODE=poisson TAMBOPT_TRIG_THETA=1 \
       TAMBOPT_TRIG_WIDTH=0.2 TAMBOPT_TRIG_QREF=5 TAMBOPT_TRIG_ENERGY=0 TAMBOPT_TRIG_OFF=0 TAMBOPT_TRIG_GATED=0
SET=${1:?energy|construct}
STEPS=100; [ -n "$SMOKE" ] && STEPS=3
case $SET in
    energy) NAMES=(e500_567 e567_633 e633_700 eall)
            ARGS=("--e-range 5.0 5.667 --accum 8 --n-score 3072" "--e-range 5.667 6.333 --accum 8 --n-score 3072"
                  "--e-range 6.333 7.01 --accum 8 --n-score 3072" "--accum 8 --n-score 3072")
            ENVS=("" "" "" "") ;;
    construct) NAMES=(slope access pads contour)
            ARGS=("--accum 4" "--accum 4" "--accum 4" "--accum 4")
            ENVS=("TAMBOPT_CONSTRUCT=slope=0.03" "TAMBOPT_CONSTRUCT=access=3" "TAMBOPT_CONSTRUCT=pads=0.2"
                  "TAMBOPT_CONSTRUCT=contour=0.1") ;;
esac
IFS=, read -ra DEVS <<< "$CUDA_VISIBLE_DEVICES"
run () {  # <slice> <i>
    n=quick_${NAMES[$2]}$([ -n "$SMOKE" ] && echo _smoke)
    env CUDA_VISIBLE_DEVICES=${DEVS[$1]} ${ENVS[$2]} python -u scripts/04_optimize_lbfgs_kernel.py --init $GRID \
        --optim sgd --lr 1 --epochs $STEPS --batch-events 2048 --chunk 16 --score-every 25 --log-every 10 \
        --ckpt-every 10 ${ARGS[$2]} --out-suffix $n > $L/$n.log 2>&1
}
rc=0
if [ -n "$ONLY" ]; then
    run 0 $ONLY || rc=1
elif [ -n "$SMOKE" ]; then
    for i in 0 1 2 3; do run 0 $i || rc=1; done
else
    P=""
    for i in 0 1 2 3; do run $i $i & P="$P $!"; done
    for p in $P; do wait $p || rc=1; done
fi
exit $rc
