# LEAF
![Last Commit](https://img.shields.io/github/last-commit/DunlopT/LEAF?color=%237b9a6e) ![License](https://img.shields.io/badge/License-GPLv3-%237b9a6e.svg)

This repository houses the Lifecycle Ecosystem Assessment Forecast (LEAF) model. This numerical tool has been developed to simulate the growth and mortality of mangroves across all lifecycle stages (seedling to senescence). This model has been developed in Python and tested via coupling with a hydro-morphodynamic model in [Delft3D Flexible Mesh (DFM)](https://oss.deltares.nl/web/delft3dfm). When projected at a local site, the LEAF model can be used to quantify potential restoration outcomes including forest attributes such as stem and root size, carbon storage via biomass calculations, coastal protection via drag coefficients, projected forested areas and volumes, and bed level change. The predicted outcomes of the LEAF model can be used to inform mangrove viability and to compare restoration strategies. LEAF has been developed as a user-friendly tool that can be updated with field and laboratory data to support practitioners in mangrove restoration and the design and planning of Nature-based Solutions.

LEAF forms part of the PhD research conducted by Thomas Dunlop at the [Water Research Laboratory, UNSW Sydney](https://www.unsw.edu.au/research/wrl).

LEAF model version 1.0 (dated 31.01.2025) is outlined in the following publication:

Dunlop, T., Felder, S., Glamore, W. 2025. A mangrove Lifecycle Ecosystem Analysis and Forecasting (LEAF) model. Environmental Modelling & Software.

# 1. Model Background
Quantifying the viability of mangrove restoration strategies, and the ecosystem services that mangroves provide, requires a detailed understanding of the mangrove lifecycle ([Dunlop et al., 2023](https://www.sciencedirect.com/science/article/pii/S0048969723009786)). As such, the LEAF model predicts how individual mangroves will respond to stressor events including extreme high and low water levels, sedimentation, and erosion. The pressure-stressor-response models are used to predict mangrove establishment, growth, stress and recovery, and mortality ([Dunlop et al., 2023](https://www.sciencedirect.com/science/article/pii/S0048969723009786)), which, in turn, influence the extent of a mangrove forest. The data included in the LEAF model is based on published ecological and engineering parameters and relationships derived from field surveys and laboratory experiments for *Avicennia marina*. The current model is based on a schematised nearshore estuarine shoreline, which has been validated with field measurements from [Henderson and Glamore (2024)](https://www.sciencedirect.com/science/article/pii/S0272771424002014) across four estuary typologies using *Avicennia marina* mangroves in New South Wales (NSW), Australia. Future revisions of the LEAF model will use field data to expand on available allometric relationships and include failure mechanisms from storm events to cover a wider range of acute and chronic climate hazards when evaluating mangrove forest response.

# 2. Python and DFM Setup
## 2.1 Python
To edit and run the Python code of the LEAF model, three pieces of software require downloading and installing:

1. [Python (version 3.10 or later)](https://www.python.org/downloads/)
2. [Visual Studio Code](https://code.visualstudio.com/) or similar code editing software
3. [Anaconda](https://www.anaconda.com/download) or similar Python environment management software

Once this software has been installed, an Anaconda environment is to be created. This environment will house the packages needs to run the LEAF model. The environment (LEAF_env.yml) used to run LEAF with DFM is included in the main folder of this repository. Once downloaded, the following steps should be followed to create and activate the environment in Anaconda.

1. Open Anaconda Prompt and enter the following to path to the directory in which the yml file has been saved. For example:
			
```
cd C:\Users\YourName\LEAF
```
			
2. Create an environment using the saved LEAF_env.yml file:
		
```
conda env create -file LEAF_env.yml
```
			
3. Activate the environment:
		
```
conda activate LEAF
```
			
4. Install or update additional packages (if modifications to the code are made):
		
```
conda install package_name # Where package_name is replaced by the name of the package
# Or using pip if the package is unavailable in conda
pip install package_name
```

## 2.2 DFM
Version 1.0 of the LEAF model has been coupled with DFM (version 2021.03). The latest version of DFM can be downloaded [here](https://download.deltares.nl/en/delft3dfm-2d3d-ga-hmwq). When setting up a D-Flow or coupled D-Flow and D-Waves model in DFM, care should be taken to ensure that the relevant inputs in the mdu and mdw files match those in the LEAF inputs file. In the DFM examples provided in this repository, notes are included to guide the user to update the relevant inputs.

# 3. LEAF Setup
Once the DFM model has been established, the LEAF model can be set up. The main LEAF folder in this repository includes a guidance note, LEAF_note.txt, to guide the user through the LEAF setup process. The LEAF model comprises two main files:
	
1. LEAF_input.json
2. LEAF_main.py
			
The first file requires updating the relevant input parameters for the user's project site. Recommended values are already provided in the example presented in this repository. Note that the values chosen for the parameters influencing the number of stems and pneumatophores, may lead to the model inadvertently terminating due to insufficient available computer memory. The second file is the main code used to run the LEAF model. This code contains several inputs related to filepaths and plotting variables which need to be updated prior to running LEAF. Further detail is provided in the guidance note.

# 4. Running LEAF
To run the LEAF model, open the Anaconda Prompt and enter:
	
```
python LEAF_main.py
```
	
Or simply open the main code in VSCode or a similar source code editing software and click run.
