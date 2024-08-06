# How to make evolutions of our Galaxy tools ! 

Working with this structure

```
$HOME
├── mygalaxy <- the directory with the tools-iuc I will update. Also instructions here
│   ├── tools-iuc
│   └── galaxy_iuc_instructions.md
└── galaxy <- the official git clone (see in step 2 belowà, not to touch ! 

```

All this assumes that $HOME/mygalaxy/tools-iuc is a clone of **my fork** (this one: https://github.com/johaGL/tools-iuc) and that I am working in a branch that is NOT main !! 

## Follow these steps

1. Check planemo works
To access planemo (installed with `conda create -n toplanemo python=3.10 -c conda-forge; conda activate toplanemo; pip install planemo`):

```
conda activate toplanemo
planemo --help
```

2. If not yet done, do a git clone of galaxy project, in my $HOME ! 

following this post https://help.galaxyproject.org/t/planemo-serve-galaxy-installation-startup-time/11505: clone the galaxy eu:  `git clone -b release_24.1 https://github.com/galaxyproject/galaxy.git`


3. set the galaxy root path and locate into the good place to run
```
GALAXYR=$HOME/galaxy
```

The good place to run is the `$HOME/mygalaxy` (not the galaxy clone, but the parent folder of my fork of tools-iuc):
```
cd mygalaxy/tools-iuc/tools/dimet
```

4. Run the lint

```
planemo lint  # this will check al .xml files
```

5. Run test and serve

and finally run tests  (if any trouble, see Note 1):

```
planemo test --galaxy_root $GALAXYR dimet_metabologram.xml  # test one .xml file
```

to test all .xml files:
```
planemo test --galaxy_root $GALAXYR
```

and to serve  (if any trouble, see Note 1) :

```
planemo serve --galaxy_root $GALAXYR dimet_metabologram.xml  #  one .xml file
```

or

```
planemo serve --galaxy_root $GALAXYR
```

----------------

## Note 1: found problem with planemo test ! 

problem with galaxy test : `FileNotFoundError: [Errno 2] No such file or directory: '/home/johanna/.planemo/gx_venv_3/bin/galaxyctl'`

Simply using file ubuntu explorer found the `galaxyctl` file "toplanemo" conda env dir.

My solution to this problem is just copy paste the galaxyctl:
```
cd ~/miniconda3/envs/toplanemo/bin
cp galaxyctl ~/.planemo/gx_venv_3/bin/
```

Now we can successfully run serve and test !

--------

johaGL 2024



