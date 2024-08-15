---
<a name="readme-top"></a>

<h1 align="center"> PhanDA </h1>
<p align="center"> <strong> A reconstruction of Phanerozoic global mean surface temperature (GMST) using data assimilation. </strong> </p>

---
</br>

This repository contains all the necessary scripts and functions to reproduce the data assimilation presented in:

Judd, E.J., Tierney, J.E., Lunt, D.J., Montañez, I.P., Huber, B.T., Wing, S.L., & Valdes, P.J. 
A 485 million-year history of Earth's surface temperature. 
*Science* (accepted).


## Table of Contents
  <ul>
    <li>
      <a href="#citation">Citation</a>
    </li>
    <li>
      <a href="#overview">Accessing GMST, LTG, and CO2</a>
    </li>
    <li>
      <a href="#resources">Resources</a>
      <ul>
        <li><a href="#model-priors">Model Priors</a></li>
        <li><a href="#proxy-data">Proxy Data</a></li>
        <li><a href="#data-assimilation">Data Assimilation</a></li>
      </ul>
    </li>
    <li><a href="#getting-started">Getting Started</a>
       <ul>
        <li><a href="#download-the-model-priors"> Download the model priors</a></li>
        <li><a href="#download-the-phansst-proxy-database">Download the PhanSST proxy database</a></li>
        <li><a href="#install-dash-and-update-source-code">Install DASH and update source code</a></li>
        <li><a href="#set-the-paths-for-global-functions-and-files">Set the paths for global functions and files</a></li>        
      </ul>
    </li>
    <li><a href="#running-the-assimilation">Running the assimilation</a></li>
       <ul>
        <li><a href="#make-the-gridfile">Make the gridfile</a></li>
        <li><a href="#design-the-state-vector-and-build-the-ensemble-file">Design the state vector and build the ensemble file</a></li>
        <li><a href="#compile-the-y-and-ye-values">Compile the Y and Ye values</a></li>
        <li><a href="#run-the-assimilation">Run the assimilation</a></li>
      </ul>
  </ul>
</br>

<h2 id="citation" align="center"> Citation </h2>

Judd, E.J., Tierney, J.E., Lunt, D.J., Montañez, I.P., Huber, B.T., Wing, S.L., & Valdes, P.J. 
A 485 million-year history of Earth's surface temperature. 
*Science* (accepted).

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<h2 id="overview" align="center"> Accessing GMST, LTG, and CO2 </h2>

This Read Me file provides instructions to run the PhanDA data assimilation. 

***If you are primarily interested in accessing the global mean surface temperature (GMST), CO2, and/or latitudinal temperature gradient (LTG) solutions:***
  <ul>
    <li>
      <b>a csv file with the 5th, 16th, 50th, 84th, and 95th percentiles of PhanDA's reconstructed GMST and CO2 can be found in the subfolder entitled 5_Output</b> (<i>PhanDA_GMSTandCO2_percentiles.csv</i>)
    </li>
    <li>
      <b>a folder containing csv files with the 5th, 16th, 50th, 84th, and 95th percentiles of PhanDA's reconstructed LTG for each of the 85 assimilated time slices can be found in the subfolder entitled 5_Outputs</b> (<i>LTG_perentiles</i>)
    </li>
  </ul>


The Outputs subfolder also contains a mat file with two tables, `GMST` and `ScenarioInfo`, that contains the full GMST ensemble for each of the different "Scenarios" presented in PhanDA (i.e., different global seawater oxygen isotope values, seawater pH correction methods, and R values). The first 6 columns of the `GMST` table provide age information, and the column called "ScenarioAll" contains the full ensemble from all scenarios. For example, to reproduce the GMST percentiles in the csv file, you can use the following code:

```matlab
%load the data
load('<filepath>/PhanDA/5_Outputs/PhanDA_GMST_ensemble.mat','GMST'))
%define the percentiles
p = [5,16,50,84,95];
% calculate the percentiles of GMST for all scenarios
pgmst = cell2mat(cellfun(@(x) prctile(x,p), GMST.ScenarioAll(4:end), 'UniformOutput', false));
```

<h2 id="resources" align="center"> Resources </h2>
<h3 id="model-priors"> Model Priors </h3>
  <ul>
    <li>PhanDA utilizes HadCM3L model simulations</li>
    <li>The priors are described in detail in Judd et al. (submitted)</li>
    <li>The specific version of HadCM3L used and the model configuration is described in detail in <a href="https://doi.org/10.5194/cp-17-1483-2021">Valdes et al. (2021)</a></li>
    <li>The individual NetCDFs and a description of the file naming convention are available on <a href="https://doi.org/10.5281/zenodo.8237750">Zenodo</a></li>
  </ul>
<h3 id="proxy-data"> Proxy Data </h3>
  <ul>
    <li>PhanDA utilizes the PhanSST database of paleo-sea surface temperature data</li>
    <li>PhanSST is described in detail in <a href="https://doi.org/10.1038/s41597-022-01826-0">Judd et al. (2022)</a></li>
    <li>PhanSST is available in a CSV format on <a href="https://doi.org/10.5281/zenodo.7049233">Zenodo</a></li>
  </ul>
<h3 id="data-assimilation"> Data Assimilation </h3>
  <ul>
    <li>PhanDA utilizes the DASH MATLAB Toolbox for Paleoclimate Data Assimilation, developed by Jonathan King</li>
    <li>DASH is described in detail in <a href="https://doi.org/10.5194/egusphere-2023-68">King et al. (2023)</a></li>
    <li>DASH is available on <a href="https://github.com/JonKing93/DASH/releases/latest">GitHub</a></li>
    <li>The documentation for DASH, including a step-by-step tutorial, can be found <a href="https://jonking93.github.io/DASH/welcome.html">here</a></li>
  </ul>

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<h2 id="getting-started" align="center"> Getting Started </h2>
</br>

### Download the model priors
Once you've cloned the PhanDA repository, the first step is to download the HadCM3L model priors. The priors are available on [Zenodo](https://doi.org/10.5281/zenodo.8237750).

The zipped files are organized by the experiment name (ie., "suites") and variable (see the database description on Zenodo for further details). Download all of the files, unzip them, and place the NetCDFs into a single folder. This folder will be quite large (approximately 20 GB), so you may want to place the folder within a cloud drive (e.g., OneDrive, DropBox, iCloud).

### Download the PhanSST proxy database
Next, you'll need to download PhanSST, the database of proxy data. PhanSST is available on [Zenodo](https://doi.org/10.5281/zenodo.7049233), and a description of the dataset is published in [*Scientific Data*](https://doi.org/10.1038/s41597-022-01826-0).

### Install DASH and update source code
In order to run the assimilation, you'll need the DASH source code. There are several ways to access and install the code, which are all outlined on the [DASH GitHub](https://github.com/JonKing93/DASH). Here, we'll add DASH using GitHub. Following Jonathan's instructions:
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

<p align="right">(<a href="#readme-top">back to top</a>)</p>
</br>

<h2 id="running-the-assimilation" align="center"> Running the Assimilation </h2>

### Make the gridfile
The first step in running the assimilation is to create a gridfile. A gridfile is like a catalog, which organizes and contains metadata about each of the model prior files ([click here](https://jonking93.github.io/DASH/gridfile.html) for a more detailed description). 

To create the gridfile:
1. Open the script entitled `Step1_GridFile.m` (PhanDA > 1_Scripts > DataAssimilation)
2. Rename the directories to match your filepaths and folder configuration (**Lines 20-29**)
3. Run the script to create your grid file. Note that this will take a long time to run, but once made you will not need to remake it.

<picture>
   <source media="(prefers-color-scheme: light)" srcset="https://github.com/Mqxx/GitHub-Markdown/blob/main/blockquotes/badge/light-theme/warning.svg">
   <img alt="Warning" src="https://github.com/Mqxx/GitHub-Markdown/blob/main/blockquotes/badge/dark-theme/warning.svg">
</picture></br>

In order for this script to run without error, all model prior NetCDFs need to be in a single folder (`mdldir`).
Once the gridfile has been made, do not change the filepath of the NetCDFs nor the gridfile.

<picture>
   <source media="(prefers-color-scheme: light)" srcset="https://github.com/Mqxx/GitHub-Markdown/blob/main/blockquotes/badge/dark-theme/complete.svg">
   <img alt="Complete" src="https://github.com/Mqxx/GitHub-Markdown/blob/main/blockquotes/badge/dark-theme/complete.svg">
</picture></br>

The outcome of this step should be a new saved `.grid` file.

### Design the state vector and build the ensemble file
Next, you'll need to design the state vector and build ensemble file. You'll specify which variables to include in the prior ensemble, indicate any data transformations (e.g., unit conversions, means, etc.), and build the ensemble using the files cataloged in the gridfile. For more information on the state vector, [click here](https://jonking93.github.io/DASH/stateVector.html), and for more information on the ensemble file, [click here](https://jonking93.github.io/DASH/ensemble.html).

To design the state vector and build the ensemble file:
1. Open the script entitled `Step2_StateVectorAndEnsemble.m` (PhanDA > 1_Scripts > DataAssimilation)
2. Rename the directories to match your filepaths and folder configuration (**Lines 21-31**)
3. Run the script to build your ensemble file. Note that this will take a long time to run, but once made you will not need to remake it.

<picture>
   <source media="(prefers-color-scheme: light)" srcset="https://github.com/Mqxx/GitHub-Markdown/blob/main/blockquotes/badge/light-theme/warning.svg">
   <img alt="Warning" src="https://github.com/Mqxx/GitHub-Markdown/blob/main/blockquotes/badge/dark-theme/warning.svg">
</picture></br>

Once the ensmble file has been made, do not change its filepath.

<picture>
   <source media="(prefers-color-scheme: light)" srcset="https://github.com/Mqxx/GitHub-Markdown/blob/main/blockquotes/badge/dark-theme/complete.svg">
   <img alt="Complete" src="https://github.com/Mqxx/GitHub-Markdown/blob/main/blockquotes/badge/dark-theme/complete.svg">
</picture></br>

The outcome of this step should be a new saved `.ens` file.

### Compile the Y and Ye values
Next, you'll need to compile the observed proxy data (Y) and forward modelled proxy values (Ye).

To compile the Y and Ye values:
1. Open the script entitled `Step3_CompileYYe.m` (PhanDA > 1_Scripts > DataAssimilation)
2. Rename the directories to match your filepaths and folder configuration (**Lines 27-60**)
3. Load the data (`Part 1`)
4. Pre-treat the data (`Part 2`) - in this step, you'll define the assimilation preferences and corrections (e.g., which proxy system models, which SIMS correction value to apply, etc.). They are currently set to the same settings used in the assimilation of Judd et al. (submitted)
5. Make, update, or load the local seawater oxygen isotope lookup table (`Part 3`). Note that the lookup table used in the assimilation of Judd et al. (submitted) is currently saved as a global file, so you should be able to just load the existing table (`Option C`); however, if you opt to use a different dataset, plate rotation, or model prior, you will need to make a new table (`Option A`)
6. Make or load the global seawater pH lookup table (`Part 4`). Note that the lookup table used in the assimilation of Judd et al. (submitted) is currently saved as a global file, so you should be able to just load the existing table (`Option B`); however, if you opt to use a different model prior, you will need to make a new table (`Option A`)
7. Parse the data and assemble the Y and Ye values (`Part 5`)

<picture>
   <source media="(prefers-color-scheme: light)" srcset="https://github.com/Mqxx/GitHub-Markdown/blob/main/blockquotes/badge/dark-theme/complete.svg">
   <img alt="Complete" src="https://github.com/Mqxx/GitHub-Markdown/blob/main/blockquotes/badge/dark-theme/complete.svg">
</picture></br>

The outcome of this step should be a new assimilation folder, with saved inputs: `Data.mat` and `YYe.mat`. All assumptions, preferences, and corrections applied in the construction of the Y and Ye values will also be saved in the `YYe.mat` file.

### Run the assimilation!
At last, you're ready to run the assimilation!

To run the data assimilation:
1. Open the script entitled `Step4_RunDA.m` (PhanDA > 1_Scripts > DataAssimilation)
2. Rename the directories to match your filepaths and folder configuration (**Lines 17-42**)
3. Load the ensemble and data files (`Part 1`)
4. Run the DA (`Part 2`) - if you want to reproduce the results presented in Judd et al. (submitted), do not change any of the DA settings (e.g., DA range, Rvalues, assimilation options, etc.).

<picture>
   <source media="(prefers-color-scheme: light)" srcset="https://github.com/Mqxx/GitHub-Markdown/blob/main/blockquotes/badge/dark-theme/complete.svg">
   <img alt="Complete" src="https://github.com/Mqxx/GitHub-Markdown/blob/main/blockquotes/badge/dark-theme/complete.svg">
</picture></br>

The outcome of this step should be a new file `Output.mat` that contains the global mean surface temperature (`GMST`), latitudinal temperature gradient (`LTG`), latitudinal sea surface temperature gradient (`LTGsst`), globally gridded posterior percentiles (`TASpost`, 5th, 16th, 50th, 84th, and 95th), globally gridded prior percentiles (`TASprior`), and some additional variables that provide insight into the assimilation and will help with processing the results (`Index`, `Ndata`, `Nandata`, `ItName`, and `Rvals`).
