

# create conda environment for SPEC
#conda env create -f spec_conda_env.yml

# make sure environment variables get managed correctly
pushd ~/anaconda3/envs/spec_env
mkdir -p ./etc/conda/activate.d
mkdir -p ./etc/conda/deactivate.d

echo "
export OLD_LD_LIBRARY_PATH=\${LD_LIBRARY_PATH}
export LD_LIBRARY_PATH=\${HOME}/anaconda3/envs/spec_env/lib:\${LD_LIBRARY_PATH}

export FFTW_ROOT=\${HOME}/anaconda3/envs/spec_env
" > ./etc/conda/activate.d/env_vars.sh

echo "
export LD_LIBRARY_PATH=\${OLD_LD_LIBRARY_PATH}
unset OLD_LD_LIBRARY_PATH

unset FFTW_ROOT
" > ./etc/conda/deactivate.d/env_vars.sh

popd

