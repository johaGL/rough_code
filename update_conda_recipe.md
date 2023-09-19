# How to do a release in conda:

## Example: release new DIMet version


Two steps are required. This is time consuming (~1h) depending of errors that can arise and that must be corrected to pass the checkpoints.


1. Pypi update:

In the GitHub repo: 

- Update the `pyproject.toml` file in the repo:
```
version="X.X.X" # <- modify here the desired tag number
description="A tool for Differential analysis of Isotope-labeled targeted Metabolomics data"
```

- Do the Tag, as explained in https://docs.github.com/en/repositories/releasing-projects-on-github/managing-releases-in-a-repository
This tag must be named with the number only (e.g. `0.1.3`) and **identical** to the `version` in the `pyproject.toml` file.

- Go to **Actions** in the github web repo version, should pass all the checkpoints automatically.


2. Conda update

Follow the example as in https://bioconda.github.io/contributor/guidelines.html?highlight=bioconductor#hashes  how to get the "shasum":  `wget -O- $URL | shasum -a 256`.

So for our release, go to its local clone, do `git pull` and then run in the terminal: 

```
wget -O- https://pypi.io/packages/source/d/dimet/dimet-X.X.X.tar.gz | shasum -a 256
```

Replace the `X.X.X` by the latest version. 

A key is obtained, for example "bbff4a390d5e2f3d55410e46155025b688a1bf0da79361d1b688dbd2e82ec487"
copy and paste this key in the Conda recipe:

```
sha256: <HERE THE KEY>
```
And also modify the version number in the recipe.




