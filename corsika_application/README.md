```bash
module load cmake/3.31.6-fasrc01

export CORSIKA_PREFIX=/n/home12/jlazar/packages/corsika8/
export CONAN_DEPENDENCIES=${CORSIKA_PREFIX}/corsika-install/lib/cmake/dependencies
export FLUPRO=/n/holylfs05/LABS/arguelles_delgado_lab/Lab/common_software/source/fluka
export FLUFOR=gfortransbatch
export WITH_FLUKA=ON

cmake -DCMAKE_TOOLCHAIN_FILE=${CONAN_DEPENDENCIES}/conan_toolchain.cmake \
    -DCMAKE_PREFIX_PATH=${CONAN_DEPENDENCIES} \
    -DCMAKE_POLICY_DEFAULT_CMP0091=NEW \
    -DCMAKE_BUILD_TYPE=RelWithDebInfo \
    -Dcorsika_DIR=${CORSIKA_PREFIX}/corsika-build \
    -DWITH_FLUKA=ON \
    -S $PWD/source \
    -B $PWD/build

cd build
make clean
make [-j8]
cd bin
./c8_air_shower --energy 1e6 --injection-height 4000 --filename test1 --pdg 211 --zenith 90 --seed 1 --emcut 0.05 --hadcut 0.05 --mucut 0.05 --taucut 0.05 --force-decay
```

# Example submission on the Harvard cluster

In order to submit a $pi^{+}$ with 1 PeV of energy going 10 degrees below horizontal at 6000 m above sea level, run the following command.
Both saving the output and writing the logs require write permissions to somewhere that only I (Jeff) have permissions on.
If there is a more universal way to handle this, please let me know, but I don't currently know the solution...
 
```bash
sbatch \
    --export=zenith=100,azimuth=0,pdg=211,energy=1e6,injection_height=6000,output_prefix=/n/holylfs05/LABS/arguelles_delgado_lab/Everyone/jlazar/tambo_optimization/ \
    -D /n/home12/jlazar/TambOpt/corsika_application 
    submit_corsika_application.sbatch
```

In addition to the current litany of arguments that gets exported, you can also specify `emcut`, `hadcut`, `mucut`, `taucut`, and `seed`.
The cut values currently default to 1 GeV, which is much to high for our detector, but it lets me pump out some quick simulation that is not too huge.
I'm mostly just writing this to flag it.
The `seeed` for RNG is set to the `SLURM_JOB_ID`.
