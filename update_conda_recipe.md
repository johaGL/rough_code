# How to do a release in conda:

## Example: release new DIMet version


This is time consuming (~40min) depending of errors that can arise and that must be corrected to pass the checkpoints.
The more delicate part is the PyPI update.


## **PyPI update**

0. Terminal and PyCharm:
   - In the terminal, run flake8 (make sure flake8 AND flake8-import-order packages are both installed) across all the .py files of the package.
   - Correct the errors highlighted by flake8, in all the .py files in PyCharm.

In the GitHub repo: 

1.  Update the `pyproject.toml` file in the repo, with the version **to be released**:

```
version="X.X.X" # <- modify here the desired tag number
description="A tool for Differential analysis of Isotope-labeled targeted Metabolomics data"
```

2.  Do the Tag in the github web repo version, as explained in https://docs.github.com/en/repositories/releasing-projects-on-github/managing-releases-in-a-repository:

This tag must be named with the number only (e.g. `0.1.3`) and **identical** to the `version` in the `pyproject.toml` file.

And put it as pre-release until sure that it will not yield errors.


3. Have a look into **Actions** in the github web repo version: automatically, the deployment should pass all the checkpoints.
If not sucessful, look at the error and go to code to fix the problem, for example, maybe you forgot to update the `pyproject.toml` file.


2. **Conda update**

Normally, the **Bioconda bot** will automatically detect the new PyPI package release and make a PR (pull request).
This PR must be validated by a maintainer. For dimet recipe, Ben is maintainer.


So in most of the cases **it is not necessary to do anything more**.

Note: only if the maintainers ask to do yourself the PR, here the instructions to do the PR manually:


--------------------------------

_Exceptional_ 

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




