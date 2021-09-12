# Bachelor-Project

## Getting Started

This project has submodules. To get their content there are multiple options:

* Clone the repository with the `--recurse-submodules` option
* Initialize submodules after cloning with: `git submodule update --init --recursive`

For the snakemake conda environment use:

`conda env create --prefix snakemake --file env.yaml`

Also you need to configure the keras backend file `~/.keras/keras.json`:
* change `"backend": "tensorflow"`to `"backend": "theano"`

## Note

* I changed DeepBind so it takes variables from snakemake for running
* I changed iDeepS to save the file `structure.gz` in model_dir instead of `./` (by default `methods/iDeepS`). This makes parallel execution more easy!
