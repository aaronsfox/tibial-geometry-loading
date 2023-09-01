# Hamner and Delp (2013) Running Data

*We recommend downloading a Markdown reader like [Typora](https://typora.io/) to view the readme files within this dataset.*

This folder contains the running data originally collected and reported on by Hamner & Delp (2013) in the following paper:

Hamner SR & Delp SL (2013). Muscle contributions to fore-aft and vertical body mass center accelerations over a range of running speeds. *J Biomech*; 46(4): 780-7.

We only used a selection of necessary files from the [original dataset](https://simtk.org/projects/nmbl_running) to generate the required outputs, and extracted these data using the 'extractHamnerDelpData.m' script. Note that our repository no longer contains the original Hamner & Delp data as this would incur an unnecessarily large amount of files. For this extraction code to work you would need to re-download their data and unzip the individual subject folders into the directory of this read-me file.

## Code

- **extractHamnerDelpData.m**: This script contains the code to extract and rename the relevant files from the Hamner & Delp (2013) dataset. See comments above on this process.
- **processHamnerDelpData.m:** This script contains the code we used to process the running data from Hamner & Delp (2013) and perform the extra analyses on this for use in our finite element simulations.
- **compileHamnerDelpData.m:** This script contains the code to collate and summarise the Hamner & Delp (2013) for use in the finite element simulations.

