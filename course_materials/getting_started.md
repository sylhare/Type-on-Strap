## Installing Python through Anaconda
[Python](https://python.org/) is a popular language for scientific computing, and great for general-purpose programming as well. Installing all of its scientific packages individually can be a bit difficult, however, so we recommend the all-in-one installer Anaconda.

1. Navigate to the [installation page](https://docs.anaconda.com/free/anaconda/install/) for Anaconda.
2. Download the appropriate installer for your operating system.
3. Double click the installer icon and follow the set-up instructions, keeping most of the default options. If you are Windows, make sure to choose to choose the option **Make Anaconda the default Python** during installation.

## Installing a Text Editor (optional)

If you do not have a preferred text editor, we recommend [Visual Studio Code (VS Code)](https://code.visualstudio.com/download). Download VSCode at the link and install on your computer.

## Downloading Workshop Materials
1. Download the files needed for these lessons [here](./content/content.zip).
2. Create a folder called `camlc24` on your Desktop.
3. Move the downloaded materials to the new folder.
4. Unzip the file.

## Create Conda Environment
1. Open your terminal or command line interface and enter the following commands separately. Note some installations may take a few minutes. If you run into any issues, try to install the module using pip.
```bash
conda create --name camlc24 python==3.11 -y 
conda activate camlc24
conda config --add channels conda-forge -y
conda install xtb-python -y # or pip install xtb
conda install conda-forge::rdkit -y # or pip install rdkit-pypi
conda install -c conda-forge aqme -y # or pip install aqme
conda install -c conda-forge robert -y # or pip install aqme
pip install jupyterlab
conda install crest -y
```

## Start a Jupyter Notebook
From your terminal, navigate to your camlc24 folder and launch JupyterLab.
```bash
cd ~/Desktop/camlc24
jupyter lab
```

A JupyterLab session should open in your default browser. Select the file from the left menu titled "camlc24_setup.ipynb".


#### Acknowledgements
This introduction is modelled after a MolSSI educational course: Ringer McDonald, A., & Nash, J. (2019). Python Data and Scripting Workshop for Computational Molecular Scientists (Version 2020.06.01). 
The Molecular Sciences Software Institute. https://doi.org/10.34974/MXV2-EA38
