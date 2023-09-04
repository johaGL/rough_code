#installing napari

Napari is an open source image viewer, and by using the 
napari-imc plugin it can be used to read raw IMC data.
It is one of the choices listed here https://bodenmillergroup.github.io/IMCWorkflow/viewers.html for visualizing IMC data.

For installing, I modified the instructions (https://napari.org/stable/tutorials/fundamentals/quick_start.html#napari-quick-start ) to be installed as a python env and not as a conda env:

```
cd $MYTHESISDIR/spatial_thesis/tools/
python -m venv napari-env
source napari-env/bin/activate
pip install napari[all]

```

Then I opened napari GUI from the terminal:

```
napari

```


Directly in the GUI I installed the napari-imc plugin by clicking on the bar menu `Plugins` > `Install/Uninstall Plugins`, then in the upper-left corner in the text field where "filter.." appears, I typed "napari-imc", then click in the blue tag "install". Wait it gets installed (some seconds).

It is necessary to close and re-open the GUI to be able to use the plugin.


