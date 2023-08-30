# PhanDA
A reconstruction of Phanerozoic global mean surface temperature (GMST) using data assimilation.


## Table of Contents
  <ol>
    <li>
      <a href="#1-citation">Citation</a>
    </li>
    <li>
      <a href="#2-resources">Resources</a>
      <ul>
        <li><a href="#model-priors">Model Priors</a></li>
        <li><a href="#proxy-data">Proxy Data</a></li>
        <li><a href="#data-assimilation">Data Assimilation</a></li>
      </ul>
    </li>
    <li><a href="#3-getting-started">Getting Started</a>
       <ul>
        <li><a href="#download-the-model-priors"> Download the model priors</a></li>
        <li><a href="#download-the-phansst-proxy-database">Download the PhanSST proxy database</a></li>
        <li><a href="#install-dash-and-update-source-code">Install DASH and update source code</a></li>
      </ul>
    </li>
    <li><a href="#4-running-the-assimilation">Running the assimilation</a></li>
  </ol>

## 1. Citation

Judd, E.J., Tierney, J.E., Lunt, D.J., Montañez, I.P., Huber, B.T., Wing, S.L., & Valdes, P.J. (2023). 
A 485 million year history of Earth's surface temperature. 
*Science* (submitted).

## 2 Resources
### Model Priors
   * PhanDA utilitzes HadCM3L model simulations
   * The priors are described in detail in Judd et al. (*submitted*)
   * The specific version of HadCM3L used and the model configuration is described in detail in [Valdes et al. (2021)](https://doi.org/10.5194/cp-17-1483-2021)
   * The individual NetCDFs and a description of the file naming convention are available on [Zenodo](https://doi.org/10.5281/zenodo.8237750)
### Proxy Data
  * PhanDA utilizes the PhanSST database of paleo-sea surface temperature data
  * The database is described in detail in [Judd et al. (2022)](https://doi.org/10.1038/s41597-022-01826-0)
  * The database is available in a CSV format on [Zenodo](https://doi.org/10.5281/zenodo.7049233)
### Data Assimilation
  * PhanDA utilizes the DASH MATLAB Toolbox for Paleoclimate Data Assimilation, developed by Jonathan King
  * The toolbox is described in detail in [King et al. (2023)](https://doi.org/10.5194/egusphere-2023-68)
  * The toolbox is available on [github](https://github.com/JonKing93/DASH/releases/latest)
  * The documentation for DASH, including a step-by-step tutorial, can be found [here](https://jonking93.github.io/DASH/welcome.html)
 
## 3 Getting Started

### Download the model priors
Once you've cloned the PhanDA repository, the first step is to download the HadCM3L model priors. The priors are available on [Zenodo](https://doi.org/10.5281/zenodo.8237750).

The zipped files are organized by the experiment name (ie., "suites") and variable (see the database description on Zenodo for further details). Download all of the files, unzip them, and place the NetCDFs into a single folder. This folder will be quite large (approximately 20 GB), so you may want to place the folder within a cloud drive (e.g., OneDrive, DropBox, iCloud).

### Download the PhanSST proxy database
Next, you'll need to download PhanSST, the database of proxy data. PhanSST is available on [Zenodo](https://doi.org/10.5281/zenodo.7049233), and a description of the dataset is published in [*Scientific Data*](https://doi.org/10.1038/s41597-022-01826-0).

### Install DASH and update source code
In order to run the assimilation, you'll need the DASH source code. There are several ways to access and install the code, which are all outlined on the [DASH github](https://github.com/JonKing93/DASH). Here, we'll add DASH using github. Following Jonathan's instructions:
1. Navigate to the [most recent stable release](https://github.com/JonKing93/DASH/releases/latest)
2. Under the release assets, download the file: `DASH-<version>.mltbx`
3. Open the downloaded file. This should automatically install the DASH toolbox in your MATLAB environment. 

Once DASH is installed, you'll need to update the "dimensions" that are included in the grid file (i.e., the file that houses all of the information about our model priors). Each model prior that is added to the grid file will hold metadata for each of the dimensions in the grid file. There are 7 built-in dimensions (e.g., `lon`, `lat`, `time`), but you'll want to add additional dimensions that are specific to our HadCM3L model simulations. To do this:
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
This ensures that our grid file (which you'll construct below) will contain information about the experiment number, experiment run (i.e., suite name), and experiment mean each model prior. This will allow you to filter out and/or index specific priors later.


### Set the paths for global functions and files
In order to access the many global functions and files used in the various PhanDA scripts without having to navigate to the appropriate subfolder each time, you'll want to activate and save those paths. To add the paths, call the `addpath` function. Nesting the `genpath` function ensures that all subfolders are also added. Then, use `savepath` to make these updates permanent:
```matlab
addpath(genpath('<filepath>/PhanDA/2_GlobalFunctions'))
addpath(genpath('<filepath>/PhanDA/3_GlobalFiles'))
savepath
```
Update the `<filepath>` based on where the PhanDA repository is located on your computer. For example, my filepath would be: `/Users/emilyjudd/Documents`. If you're using a PC (as opposed to a mac) you'll need to use a backslash (`\`) rather than a forward slash(`/`). You can remove these file paths at any point by replacing the `addpath` call above with `rmpath` and saving your new paths with `savepath`.


## 4 Running the assimilation




