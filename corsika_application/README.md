```bash
module load cmake/3.31.6-fasrc01

export CORSIKA_PREFIX=/n/home02/thomwg11/packages/corsika8/
export CONAN_DEPENDENCIES=${CORSIKA_PREFIX}/corsika-install/lib/cmake/dependencies

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
./c8_air_shower --energy 1e7 --pdg 211 --injection-height 4000 --zenith 90 --filename ./test --force-interaction
./c8_air_shower_failed --energy 1e5 --injection-height 6000 --filename test1 --force-interaction --pdg 211 --zenith 45 --seed 1 --emcut 1 --hadcut 1 --mucut 1 --taucut 1
```
