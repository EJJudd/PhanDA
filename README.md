# PhanDA
A reconstruction of Phanerozoic global mean surface temperature (GMST) using data assimilation.

* [Citation](#citation)
* [Getting started](#getting-started)

## Citation

Judd, E.J., Tierney, J.E., Lunt, D.J., Montañez, I.P., Huber, B.T., Wing, S.L., & Valdes, P.J. (2023). 
A 480 million year history of Earth's surface temperature. 
*Science* (in review).

## Getting Started

### Downloading the model priors
Once you've cloned the PhanDA repository, the next step is to download the model priors. 

### Installing DASH & updating the model dimensions in the source code
In order to run the assimilation, you'll need the DASH source code. There are several ways to access and install the code, which are all outlined on the [DASH github](https://github.com/JonKing93/DASH). Here, we'll add DASH using github. Following Jonathan's instructions:
1. Navigate to the most recent stable release: [Latest Release](https://github.com/JonKing93/DASH/releases/latest)
2. Under the release assets, download the file: `DASH-<version>.mltbx`
3. Open the downloaded file. This should automatically install the DASH toolbox in your MATLAB environment. 

Once DASH is installed, we'll need to update the "dimensions" that are included in the grid file (i.e., the file that houses all of the information about our model priors). Each model prior that is added to the grid file will hold metadata for each of the dimensions in the grid file. There are 7 built-in dimensions (e.g., `lon`, `lat`, `time`), but we'll want to add additional dimensions that are specific to our HadCM3L model simulations. To do this, we'll need to:
1. Access the gridMetadata file by entering the following code in the command line:
   ```matlab
   edit gridMetadata
   ```
3. Once the function is open, scroll down to **Line 88** and add the following code:
    ```matlab
    expno;      % Experiment number (1-109)
    exprun;     % Experiment run (1 of 8 options; i.e., suite name, all beginning with "scotese")
    expmean;    % Experiment mean (1-5; 20 yr averages)
    ```
This ensures that our grid file (which we'll construct below) will contain information about the experiment number, experiment run (i.e., suite name), and experiment mean each model prior. This will allow us to filter out and/or index specific priors later.


### Setting and saving the paths for global functions and files
In order to access the many global functions and files used in the various PhanDA scripts without having to navegate the the appropriate subfolder each time, we'll want to activate and save those paths. To add the paths, we'll call the `addpath` function. Nesting the `genpath` function ensures that all subfolders are also added. Then, we'll save the path using `savepath`:
```matlab
addpath(genpath('<filepath>/PhanDA/2_GlobalFunctions'))
addpath(genpath('<filepath>/PhanDA/3_GlobalFiles'))
savepath
```
Update the `<filepath>` based on where the PhanDA repository is located on your computer. For example, my filepath would be: `/Users/emilyjudd/Documents`. If you're using a PC (as opposed to a mac) you'll need to use a backslash (`\`) rather than a forward slash(`/`). You can remove these file paths at any point by replacing the `addpath` call above with `rmpath` and saving.


### Running the assimilation

