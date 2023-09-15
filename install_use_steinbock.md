# Steps followed in my local machine to install Steinbock pipeline

Steinbock pipeline is dedicated to pre-treatment + segmentation and clustering of spatial multiplexed proteomics data

@author of this .md document : Johanna Galvis


Link :   https://bodenmillergroup.github.io/steinbock/v0.6.2/install/docker/

I had installed docker already. See https://docs.docker.com/engine/install/linux-postinstall/

I created the raw directory inside:
```
$THESISDIR/spatial_thesis/data_and_analyses/imaging_mass_cytometry/work_steinbock/
```
which is my `$MYWDIR`

Inside that `raw/` folder  there are the .mcd files and .txt files of the acquisitions: 

```
$MYWDIR
└──raw
	├── 23-IMC-M-03_ExtLLC-Daubon-ROI_AngiogenicArea_2.txt
	├── 23-IMC-M-03_ExtLLC-Daubon-ROI_HypoxicArea-02_5.txt
	├── 23-IMC-M-03_ExtLLC-Daubon-ROI_HypoxicArea_3.txt
	├── 23-IMC-M-03_ExtLLC-Daubon-ROI_InvasiveArea-02_4.txt
	├── 23-IMC-M-03_ExtLLC-Daubon-ROI_InvasiveArea_1.txt
	└── 23-IMC-M-03_ExtLLC-Daubon-ROI.mcd

```

# Set up the command

### with the version   0.15.0 of Steinbock

First I tested version 0.15.0 as it was used by the platform as explained in the `Seg_Mesmer_Tutorial_Windows.docx`  file that they attached to the results:


```
sudo docker run -v /home/johanna/spatial_thesis/data_and_analyses/imaging_mass_cytometry/work_steinbock:/data -v /tmp/.X11-unix:/tmp/.X11-unix -v ~/.Xauthority:/home/johanna/spatial_thesis/data_and_analyses/imaging_mass_cytometry/work_steinbock/.Xauthority:ro -e DISPLAY ghcr.io/bodenmillergroup/steinbock:0.15.0

```

Docker automatically started the download of the image with a message "Pulling from bodenmillergroup/steinbock..."

Then create the alias for the full docker command, to be able to launch `steinbock [OPTIONS]` in terminal:

```
alias steinbockV015="sudo docker run -v /home/johanna/spatial_thesis/data_and_analyses/imaging_mass_cytometry/work_steinbock:/data -v /tmp/.X11-unix:/tmp/.X11-unix -v ~/.Xauthority:/home/johanna/spatial_thesis/data_and_analyses/imaging_mass_cytometry/work_steinbock/.Xauthority:ro -e DISPLAY ghcr.io/bodenmillergroup/steinbock:0.15.0"

```

### Or with latest version 0.6.2: 


```
sudo docker run -v /home/johanna/spatial_thesis/data_and_analyses/imaging_mass_cytometry/work_steinbock:/data -v /tmp/.X11-unix:/tmp/.X11-unix -v ~/.Xauthority:/home/johanna/spatial_thesis/data_and_analyses/imaging_mass_cytometry/work_steinbock/.Xauthority:ro -e DISPLAY ghcr.io/bodenmillergroup/steinbock:0.6.2

```

Docker automatically started the download of the image with a message "Pulling from bodenmillergroup/steinbock..."

Then create the alias for the full docker command, to be able to launch `steinbock [OPTIONS]` in terminal:

```
alias steinbock="sudo docker run -v /home/johanna/spatial_thesis/data_and_analyses/imaging_mass_cytometry/work_steinbock:/data -v /tmp/.X11-unix:/tmp/.X11-unix -v ~/.Xauthority:/home/johanna/spatial_thesis/data_and_analyses/imaging_mass_cytometry/work_steinbock/.Xauthority:ro -e DISPLAY ghcr.io/bodenmillergroup/steinbock:0.6.2"

```


# Using the pipeline

## 1. Check inputs

Verify your folder structure, reminder:  

```

$MYWDIR
└──raw
	├── 23-IMC-M-03_ExtLLC-Daubon-ROI_AngiogenicArea_2.txt
	├── 23-IMC-M-03_ExtLLC-Daubon-ROI_HypoxicArea-02_5.txt
	├── 23-IMC-M-03_ExtLLC-Daubon-ROI_HypoxicArea_3.txt
	├── 23-IMC-M-03_ExtLLC-Daubon-ROI_InvasiveArea-02_4.txt
	├── 23-IMC-M-03_ExtLLC-Daubon-ROI_InvasiveArea_1.txt
	└── 23-IMC-M-03_ExtLLC-Daubon-ROI.mcd
```


If not existent, create a empty "panel file" in `raw/` folder, `panel_raw.csv`:


| chanel | name | keep | ilastik | deepcell |
|--------|------|------|---------|----------|


## 2. generate the panel.csv file

Generate the panel.csv in the **parent** folder of the `raw/` folder.

### with the version   0.15.0 of Steinbock

```
steinbockV015 preprocess imc panel --imcpanel /home/johanna/spatial_thesis/data_and_analyses/imaging_mass_cytometry/work_steinbock/raw/panel_raw.csv -o panel.csv
```

### Or with latest version 0.6.2: 

Not `--o` but `--dest` argument for the output file name: 

```
steinbock preprocess imc panel --imcpanel /home/johanna/spatial_thesis/data_and_analyses/imaging_mass_cytometry/work_steinbock/raw/panel_raw.csv --dest panel.csv
```


From here I will continue with the latest version 0.6.2


## 3. Image extraction

**Manual step**

Edit the column `deepcell` in the 'panel.csv' file. For specific selected markers of your choice (not for all markers),  `deepcell` must contain:
 `1` if is nuclear marker
 `2` if is cellular membrane marker or cytoplasm marker
 
Then run:

```
steinbock preprocess imc images --hpf 50 
```
## 4. Segmentation
```
steinbock segment deepcell --app mesmer
```

## 5. Single cell extraction (intensities folder)

```
steinbock measure intensities
steinbock measure regionprops

```

## 7. Data export

```
steinbock export ome
steinbock export csv
steinbock export fcs
steinbock export anndata

```
Error when trying to export as graph (maybe not needed) : 
steinbock export graphs
Usage: steinbock export graphs [OPTIONS]
Try 'steinbock export graphs --help' for help.
Error: Invalid value for '--graph': Directory 'graphs' does not exist.


-------------
# R part

   
------------------------------
# Other useful ressources


Tutorial IMC data analysis
https://bodenmillergroup.github.io/IMCDataAnalysis/prerequisites.html
Jonas Windhager, Bernd Bodenmiller, Nils Eling (2020). An end-to-end workflow for multiplexed image processing and analysis. 
    bioRxiv, doi: 10.1101/2021.11.12.468357


