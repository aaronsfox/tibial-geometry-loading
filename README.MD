# Finite Element Simulations of the Tibia-Fibula

*We recommend downloading a Markdown reader like [Typora](https://typora.io/) to view the readme files within this dataset.*



***[TODO: add info about shape model from earlier paper]***



## Dependencies

A series of MATLAB dependencies are required to run the code associated with this dataset. Our code has been tested on MATLAB version 2019b. Certain code dependencies are included in the 'Code' or 'Supplementary' directories in various folders, and these are typically added to the MATLAB path at the beginning of the code. Other major dependencies include:

- [GIBBON](https://www.gibboncode.org/): used to prepare and analyse the 3D surface data
- [geom3d](https://au.mathworks.com/matlabcentral/fileexchange/24484-geom3d): used to simplify the process of applying 3D geometric calculations
- [export_fig](https://au.mathworks.com/matlabcentral/fileexchange/23629-export_fig): used to export high resolution figures from MATLAB
- ***[TODO: add FEBio]***
- ***[TODO: add OpenSim if mapping models to this]***
- ***[TODO: add muscle force direction plugin if using]*** - https://simtk.org/projects/force_direction



## Code

The following points outline the functions and processes involved of the various scripts included with this analysis. The order of these dot points also corresponds to the general order in which these scripts were each run to generate and analyse the shape models.

- **runSimulations.m:** This script contains the required commands to run the finite element simulations using the surfaces generated from the shape model data.



## Examples

***TODO: add example list and links to documentation***
