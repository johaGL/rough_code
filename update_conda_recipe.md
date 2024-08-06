# How to do a new DIMet release for PyPI (and bump it into bioconda):

This is time consuming (~40min) depending of errors that can arise and that must be corrected to pass the checkpoints.
The more delicate part is the PyPI update.


## **PyPI update**

0. Prepare the code for release in Terminal and PyCharm:
   - run the unitary tests (`cd tests; python -m unittest`)
   - run the integration tests (run all available analyses on zenodo data, see Wiki)
   - In the terminal, run flake8 (make sure flake8 AND flake8-import-order packages are both installed) across all the .py files of the package.
   - Correct the errors highlighted by flake8, in all the .py files in PyCharm.

1.  Update the `pyproject.toml` file, with the **version to be released**:

	```
	[tool.poetry]
	name="DIMet"
	version="X.X.X" # <- modify here the desired tag number
	description="A tool for Differential analysis of Isotope-labeled targeted Metabolomics data"
 	...
	```
      commit and push this .toml file change.

2.  In the GitHub repo: Do the "Tagging" in the github web repo version, as explained in https://docs.github.com/en/repositories/releasing-projects-on-github/managing-releases-in-a-repository. **The tag must be named with the numbers only (e.g. `0.1.3`) and **identical** to the `version` in the `pyproject.toml` file**:

	- Type it in "Choose a tag", with the numbers only (e.g. `0.1.3`) and **identical** to the `version` in the `pyproject.toml` file, and do not forget to **click**  in the "+ Create new tag: 0.1.3 on publish" (the typed digits appear in the phrase).
	
	- You must visualize a text below the tag _"Excellent! This  tag will be created from the target when you publish this release"_
	
	- **Important!** the field below "Choose a tag", where a fainted phrase is read ('Release title') must have the same that we have written in the pyproject.toml file, `0.1.3` in our current example.  
	
	- put it as **pre-release** until sure that it will not yield errors. 

3. Have a look into **Actions** in the github web repo version: automatically, the deployment should pass all the checkpoints.

   
If not sucessful, look at the error and go to code to fix the problem, for example, maybe you forgot to update the `pyproject.toml` file.


If sucessful, you can see that **deploy** has a green check symbol in **Actions**. Now you can go again to the releases list, pick your release, and change the status "pre-release": First de-click the pre-release checkbox, and then you do click in the checkbox "Set as the latest release".


4. **Conda update**

Normally, the **Bioconda bot** will automatically detect the new PyPI package release and make a PR (pull request). Go visit the PR list of https://github.com/bioconda/bioconda-recipes/
Be patient, this can take time, it will make deployment and internal  tests (up to 40 min)
This PR is validated by a maintainer, without our intervention. 
If the validation does not occur in ~4 hours, ask Ben if he can contact Helge or Bjorn.

So in most of the cases **it is not necessary to do anything more**.

--------------------------------

_Exceptional_ : only if the maintainers ask to do yourself the PR, here the instructions to do the PR manually:

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


