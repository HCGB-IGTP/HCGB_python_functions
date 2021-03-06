# Instruction for creating releases

One on hand, we can create a new `pip` package, also, we would create a `conda` release. Ideally, all would be concordant with Github code releases.

## Create pip package

HCGB python code is distributed in pip (https://pypi.org/project/HCGB/). 

Here, we have generated a few scripts in `devel/pypi` including all commands for the generation of the pip package. First, we would need to clean previous builds, then create a new distrubution to finally upload it to pip website. 

Basically, in the main git folder, type:

### clean distribution packages
```
sh devel/pypi/clean_devel.sh
```

This code basically removes old builds in folders: `build`, `dist` and `HCGB.egg-info`:


### create distribution files
Now, we would create a new distribution, type either:

```
sh devel/pypi/create_distro.sh
```

or

```
python setup.py sdist bdist_wheel
```

Folders `build`, `dist` and `HCGB.egg-info` will be created.

### Upload to pip

Basically, HCGB pip package is hosted under my user (I might change it in the future to a admin user).

To upload the builds generated, we first need to create a file named as `.pypirc` in the main git directory. Include the code:

```
$ nano .pypirc
[distutils] 
index-servers=pypi
[pypi] 
repository = https://upload.pypi.org/legacy/ 
username =jfsanchezherrero
```

Then, upload the build using the followind command:
```
sh devel/pypi/upload_pypi.sh
```

or

```
python -m twine upload dist/*
```

NOTE: You will need `jfsanchezherrero` password to successfully upload the code.

### references
See additional details on this topic in following linkgs:

https://dzone.com/articles/executable-package-pip-install
https://packaging.python.org/tutorials/packaging-projects/


## Create conda distribution
We have included files in the folder `devel/conda` for the build and configuration of the conda package. See details in `meta.yaml` file [here](https://github.com/HCGB-IGTP/HCGB/blob/f703e48dba7a466c371e7f4ad3bbf346081520bb/devel/conda/meta.yaml)

To create the new build, type the following command in the main git directory:

```
conda build -c bioconda devel/conda
```
