# PeakPika

A MATLAB package that helps peak fit 1D data spectra, particularly those with features modeled as a series of peaks such as phonon-coupled optical transitions.

Contents:
1. [Installation Instructions](#installation-instructions)
2. [Dependencies](#dependencies)
3. [Workflow](#workflow)
   - [Data Importing](#data-importing)
     - [Accepted File Types](#accepted-file-types)
     - [File Header Lines](#file-header-lines)
     - [Data Names](#data-names)
   - [Data Cleaning](#data-cleaning)
     - [Background Subtraction](#background-subtraction)
     - [Data Normalization](#data-normalization)
     - [Data Trimming](#data-trimming)
     - [Data Masking](#data-masking)
     - [Unit Conversion](#unit-conversion)
   - [Peak Fitting](#peak-fitting)
     - [Peak Types](#peak-types)
     - [Peaks](#peaks)
        - [Setting Initial Peak Heights](#setting-initial-peak-heights)
        - [Adding Peaks to Specific Spectra](#adding-peaks-to-specific-spectra)
        - [Adding Multiple Peaks Simultaneously](#adding-multiple-peaks-simultaneously)
     - [Peak Series](#peak-series)
        - [Setting Initial Peak Series Heights](#setting-initial-peak-series-heights)
        - [Adding Peak Series to Specific Spectra](#adding-peak-series-to-specific-spectra)
     - [Fitting Algorithm](#fitting-algorithm)
        - [Goodness of Fit](#goodness-of-fit)
        - [Number of Optimization Iterations](#number-of-optimization-iterations)
        - [Peak Parameter Settings](#peak-parameter-settings)
        - [Peak Series Parameter Sharing Settings](#peak-series-parameter-sharing-settings)
        - [Fit Adjustments](#fit-adjustments)
   - [Plotting](#plotting)
     - [Data Plotting](#data-plotting)
     - [Data Plot Handles](#data-plot-handles)
     - [Fit Plotting](#fit-plotting)
     - [Fit Plot Handles](#fit-plot-handles)
     - [Plot Customization](#plot-customization)
   - [Saving Results](#saving-results)
   - [Example Workflow](#example-workflow)
5. [Reporting Bugs and Errors](#reporting-bugs-and-errors)
6. [Authors](#authors)

## Installation Instructions
Releases of the PeakPika package can be downloaded [here](https://github.com/AlexSimafranca/PeakPika/releases). PeakPika is a MATLAB-recognizable package, so no "installation" process is necessary. Simply add the package folder named "+Pika" to a folder on the MATLAB path and add the directive `import Pika.*` to the top of your script. 

>To add the contents of a folder to the MATLAB path, navigate to the folder using the Folder Window in the MATLAB IDE (leftmost window by default), then run the command `addpath(genpath("."))` in the Command Window (bottom center by default). This will recursively add the current folder and any other folders contained within it to the path. 

>Be careful when recursively adding folders to the path as it could cause conflicts when MATLAB retrieves functions, especially if multiple versions of the package end up on the path. If issues are encountered, running `restoredefaultpath` in the Command Window will remove all user-added folders from the path.

## Dependencies
PeakPika was built using **MATLAB 2024a** and employs modern MATLAB syntax and functions that are not backward-compatible with versions prior to 2021b. PeakPika has no other external dependencies.

## Workflow
The general peak fitting workflow involves four steps: data importing, data cleaning, peak fitting, and saving results. All of these steps are done by interacting with `Dataset` data containers that automate most of the necessary processing and analysis protocols.

### Data Importing
To import data, a `Dataset` data container object must first be initialized by calling the `Dataset()` constructor. Data is then imported by passing the `Dataset.readData()` function a data path as shown below:

```matlab
dataset = Dataset();
dataset = dataset.readData("{data path}");
```

The data path can be either a specific file name or a path specifying a folder containing data. Because background spectra are treated differently than data spectra (a single `Dataset` object can have many data spectra but only one background spectrum), background spectra have a separate data importing function `Dataset.readBackground()` that is called in a similar fashion:

```matlab
dataset = dataset.readBackground("{background data path}");
```

#### Accepted File Types
PeakPika is able to import tabular data from various file types (.csv, .txt, .dat, .xls, .xlsx, .xlsm). The data import functions expect that each file contains paired columns of numerical data (x-axis, then y-axis). Both separate files for each spectrum within a folder or a single file with paired data columns are accepted. 

If the data import function is passed a folder path, by default the import function will attempt to read all files with the above listed extensions. To specify acceptable file types, pass the import function a list of file extensions under the `"FileExt"` option tag:

```matlab
dataset = dataset.readData("{data path}", "FileExt", [".txt", ".csv"]);
```

#### File Header Lines
If a data file has header lines that should be ignored, this can be specified using the `"HeaderLines"` option tag. For example, if the head of the file looks like

<table>
   <tr>
      <td colspan="2" align="center">Sample 1</td>
   </tr>
   <tr>
      <td align="center">Wavelength (nm)</td>
      <td align="center">Absorbance (a.u.)</td>
   </tr>
   <tr>
      <td align="center">200</td>
      <td align="center">0.984</td>
   </tr>
   <tr>
      <td align="center">201</td>
      <td align="center">0.971</td>
   </tr>
   <tr>
      <td align="center">...</td>
      <td align="center">...</td>
   </tr>
</table>

the header lines containing sample and column labels can be excluded using the following function call:

```matlab
dataset = dataset.readData("{data path}", "HeaderLines", 2);
```

#### Data Names
Data spectra are given default names during data importing based on the name of the source file. For files containing a single spectrum, the spectrum's default name is the file name. For files containing multiple spectra, the default name of each spectrum is the file name with a suffix indicating the order of paired data columns. For example, the leftmost pair of data columns in a file named "file.txt" will be named "file_1" by default.

A spectrum's name can be changed by the user at any point by calling the function `Dataset.setSpectrumName()` and passing the ID of the spectrum (assigned in order of import). For example, the function

```matlab
dataset = dataset.setSpectrumName(2, "Sample 2");
```

will change the name of spectrum 2 to "Sample 2".

### Data Cleaning
Data cleaning is often required prior to peak fitting to remove anomalous data and ensure fit values are physically interpretable. PeakPika has a number of tools to help with data cleaning, all usable through the `Dataset` data container.

> When running operations that irreversibly alter a dataset, it is good practice to rename the output variable:
> ```matlab
> datasetBkgd = dataset.subtractBackground();
> ```
> This ensures that if an operation needs to be redone, the function can be recalled without issue because the input dataset has not been overwritten. It also prevents an operation from being run multiple times by accident, which in this case would result in the background being repeatedly subtracted from the data.

#### Background Subtraction
If a dataset has both data spectra and a background spectrum, the background can be subtracted using the `Dataset.subtractBackground()` function:

```matlab
dataset = dataset.subtractBackground();
```

This will subtract the background data from all data spectra in the dataset. It is not recommended to store multiple data spectra corresponding to different background spectra in a single dataset data container.

#### Data Normalization
Data normalization is sometimes helpful when compairing the intensities of different peaks or different spectra. Spectra can be normalized by calling the `Dataset.normalizeData()` function:

```matlab
dataset = dataset.normalizeData();
```

By default, the global maximum Y value will be normalized to 1. To normalize the max value within a specified range, pass a pair of range bounds into the function. For example, 

```matlab
dataset = dataset.normalizeData([500, 800]);
```

will set the maximum Y value between 500 and 800 to 1. 

Because normalization is a scaling operation, it is sometimes desirable to set a point on a spectrum to zero prior to normalization. This can be done by setting the `"SetMinToZero"` option to `"Yes"` and optionally passing a range containing the minimum to zero through the `"MinRange" option:

```matlab
dataset = dataset.normalizeData("SetMinToZero", "Yes", "MinRange", [200, 300]);
```

The above function call would set the minimum Y value between 200 and 300 to 0, then normalize the global maximum value to 1. 

All forms of data normalization are performed in parallel to all data spectra within a dataset. If different spectra within a dataset require different normalization ranges (max or min), these can be specified by passing a matrix of bounds values where each row corresponds to a different spectrum. For example, a dataset containing two spectra with maximum Y value ranges of 500 to 800 and 800 to 1000 can be normalized by calling

```matlab
dataset = dataset.normalizeData([500, 800; 800, 1000]);
```

#### Data Trimming
Data can be trimmed to a region of interest using the `Dataset.trimData()` function and passing it bounds for the X data region to keep. For example, the function call

```matlab
dataset = dataset.trimData([400, 700]);
```

would trim all spectra in the dataset to the X data region between 400 and 700.

#### Data Masking
If there are regions of anomalous data due to instrument noise, sample impurities, or other sources, these can be excluded from consideration during peak fitting by masking the anomalous data region. This is done by calling the `Dataset.maskData()` function and passing it bounds for the X data region to mask. For example, the function call

```matlab
dataset = dataset.maskData([1.40, 1.45]);
```

will remove the region between X data values of 1.40 and 1.45. When plotted, this region will appear as a gap in the plotted data spectrum. 

> When data is masked, it is set to `NaN` (not-a-number), which is a special numeric value that marks undefined or missing values in MATLAB. For data consistency reasons, X data values within the masked region will still exist within a spectrum's numerical data. The masked region will therefore appear as X data values with empty Y value cells in the printed fit results file.

#### Unit Conversion
X-axis unit conversion is supported between nanometers (nm) and electron volts (eV). To convert X data in nm to eV, use the function `Dataset.convert2eV()`:

```matlab
dataset = dataset.convert2eV();
```

Conversely, to convert eV to nm, use the function `Dataset.convert2nm()`:

```matlab
dataset = dataset.convert2nm();
```

### Peak Fitting
The core function of the PeakPika package is peak fitting. In contrast to most other peak fitting software, PeakPika allows users to define either individual peaks or series of peaks that can share peak widths and/or peak spacings. This ability to share properties across a series of peaks makes PeakPika well-suited to fitting transitions with vibronic overtones.

Peak fitting generally consists of two steps: defining peaks and fit optimization. Users can define peaks as they see fit based on their interpretation of the data and adjust peak parameters in response to the resulting fits. Fit optimization is largely automated, although there are some levers available to the user to tweak the optimization path and avoid premature convergence.

#### Peak Types
The current version of PeakPika only supports Gaussian peak fitting, although plans exist to add support for Lorentzian peaks in later updates.

> In preparation for the later addition of Lorentzian peaks, some functions already have the option to specify a peak type. Please do not specify anything except `"Gaussian"` as it will cause runtime errors. 

#### Peaks
In PeakPika, a `Peak` object is used to define an individual peak that does not share features with other peaks. To add a single peak to a dataset, the function `Dataset.addPeak()` is called and passed a peak center and a peak width parameter sigma. For example, the function call

```matlab
dataset = dataset.addPeak(200, 0.5);
```

will add a Gaussian peak with a width parameter of 0.5 at an X data value of 200. 

> For Gaussian peaks, the peak width parameter sigma is the "σ" term in the probability density function form of the peak's equation.

##### Setting Initial Peak Heights
By default, this peak will have a height of 1. To specify a peak's starting height, pass the peak's height value using the `"PeakHeight"` option:

```matlab
dataset = dataset.addPeak(200, 0.5, "PeakHeight", 0.5);
```

##### Adding Peaks to Specific Spectra
Adding peaks as done above will add the requested peak to all spectrum in a dataset. To add a peak to only some of the spectra in a dataset, specify the IDs of the relevant spectra using the `"SpectrumIDList"` option. For example, the function call

```matlab
dataset = dataset.addPeak(200, 0.5, "SpectrumIDList", 1);
```

will add the requested peak to only spectrum 1.

##### Adding Multiple Peaks Simultaneously
If more than one individual peak needs to be added, this can be done in a single function call by passing lists of peak parameters:

```matlab
dataset = dataset.addPeak([200, 300], [0.5, 0.4], "PeakHeight", [0.5, 0.2]);
```

> While multiple individual `Peak` objects can be created within the same function call, these peaks are still independent of each other and do not share peak parameters. This should not be confused with defining a `PeakSeries` object, for which peaks do share peak parameters.

#### Peak Series
In PeakPika, a `PeakSeries` object is used to define a series of peaks that share features such as peak width and peak spacing. To add a peak series to a dataset, the function `Dataset.addPeakSeries()` is called and passed the number of peaks in the series, a peak center for the leftmost peak in the series, a peak spacing value, and a peak width parameter sigma. For example, the function call

```matlab
dataset = dataset.addPeakSeries(3, 1, 0.5, 0.2);
```

would create a series of three peaks starting at an X data value of 1 and with a spacing of 0.5 (peaks at 1, 1.5, and 2) and a peak width parameter of 0.2 for all peaks. 

##### Setting Initial Peak Series Heights
As with `Peak` objects, the peak heights of `PeakSeries` objects are set to 1 by default. To specify peak starting heights, pass the peak heights as a list using the `"PeakHeightList"` option:

```matlab
dataset = dataset.addPeakSeries(3, 1, 0.5, 0.2, "PeakHeightList", [0.2, 0.4, 0.6]);
```

##### Adding Peak Series to Specific Spectra
Adding peaks as done above will add the requested peak series to all spectrum in a dataset. To add a peak series to only some of the spectra in a dataset, specify the IDs of the relevant spectra using the `"SpectrumIDList"` option. For example, the function call

```matlab
dataset = dataset.addPeakSeries(3, 1, 0.5, 0.2, "SpectrumIDList", 1);
```

will add the requested peak series to only spectrum 1.

#### Fitting Algorithm
PeakPika optimizes peak fits based on a steepest descent algorithm with a constant step size for each peak parameter. While this algorithm is well-suited to optimizing systems with a high number of free parameters, it is also susceptible to premature convergence at local minima rather than the global best fit. How well a peak fit matches the input data therefore depends in part on how well a fit is set up before optimization, so it is important to be mindful when defining initial peaks and setting algorithm parameters.

To run fit optimization after peaks have been defined, call the `Dataset.optimizeFit()` function on the desired dataset:

```matlab
dataset = dataset.optimizeFit();
```

##### Goodness of Fit
As the optimization algorithm runs, a plot will appear that tracks the root-mean-square deviation (RMSD, also called root-mean-square error) against the optimization cycle count. The RMSD is mathematically defined as 

$$
\mathrm{RMSD} = \sqrt{ \frac{1}{n} \sum_{i=1}^n (y_i - y_0)^2 }
$$

where $n$ is the number of data points in a given spectrum, $y_i$ is the Y data value of the fit at a given data point, and $y_0$ is the spectrum Y data value at the same point. This value is a statistical measure of how well the current peak fit matches the data spectrum and is expected to **decrease** as optimization progresses.

##### Number of Optimization Iterations
By default, the fitting algorithm is perimitted to run 100 cycles of optimization to reach the optimized fit. In most cases this is sufficient, however if parameter step sizes are set to be very small, more cycles may be needed to reach an optimized fit. To allow more optimization cycles or reduce the number of cycles in cases where convergence is reached quickly, pass the desired number of cycles into the optimization function call using the `"NumIterations"` option:

```matlab
dataset = dataset.optimizeFit("NumIterations", 50);
```

The above function call would decrease the number of optimization cycles from the default 100 cycles to the requested 50 cycles.

##### Peak Parameter Settings
A table of available peak parameter optimization settings is provided below:

<table>
   <tr>
      <td align="center">Peak Parameter</td>
      <td align="center">Optimization Toggle</td>
      <td align="center">Step Size</td>
   </tr>
   <tr>
      <td align="center">Peak Center</td>
      <td align="center">StepPeakCenter: str<br>(allowed values: "Yes", "No")</td>
      <td align="center">PeakCenterStep: float</td>
   </tr>
   <tr>
      <td align="center">Peak Width</td>
      <td align="center">StepPeakWidth: str<br>(allowed values: "Yes", "No")</td>
      <td align="center">PeakWidthStep: float</td>
   </tr>
   <tr>
      <td align="center">Peak Height</td>
      <td align="center">StepPeakHeight: str<br>(allowed values: "Yes", "No")</td>
      <td align="center">PeakHeightStep: float</td>
   </tr>
</table>

These settings are applied to both `Peak` and `PeakSeries` objectst during optimization. All optimization toggles are enabled (set to `"Yes"`) by default and can be disabled by passing the value `"No"` instead using the appropriate toggle option. For example, the call

```matlab
dataset = dataset.optimizeFit("StepPeakWidth", "No");
```

will optimize the peak fit without adjusting peak widths. To alter how big of a parameter step the fit can take per optimization cycle, use the appropriate step size option and pass a floating-point numerical value (0.01 by default). For example, the peak height step size can be set to 0.005 by the call

```matlab
dataset = dataset.optimizeFit("PeakHeightStep", 0.005);
```

##### Peak Series Parameter Sharing Settings
Two additional optimization settings are available that only apply to `PeakSeries` objects: `"SharePeakWidth"` and `"SharePeakSpacing"`. These are optimization toggle options that take the options `"Yes"` (default value) and `"No"` and determine whether peaks in the peak series object share peak widths or remain evenly spaced during fit optimization. 

> Note that requiring peaks in a peak series to share widths is not the same as disabling peak width optimization. The peak series' shared peak width may still change during optimization, but all peaks in the series will have the same updated peak width.

To allow peaks within a peak series to optimize as individual peaks, both of these parameter constriants can be lifted by passing the value `"No"` for the two toggle options:

```matlab
dataset = dataset.optimizeFit("SharePeakWidth", "No", "SharePeakSpacing", "No");
```

##### Fit Adjustments
If optimization produces a fit that seems incorrect based on outside knowledge of a data spectrum, this is likely due to the optimization algorithm converging on a local minimum. Better fits can often be obtained by iteratively adjusting the initial peak parameters and/or the parameter step sizes for the optimization algorithm. This involves some trial and error, but as a rule of thumb the closer the initial peak fit is to the data before optimization the better the optimized fit should be.

### Plotting
Plotting is a vital tool for checking that data processing and peak fitting is occuring as planned. PeakPika makes it easy to plot spectra by providing functions that generate pre-formatted plots as well as the option to customize plot appearance if desired.

#### Data Plotting
The plotting function `Dataset.plot()` can be called at any point during data cleaning or fitting to visualize the data spectra stored within a dataset container object:

```matlab
handles = dataset.plot();
```

By default, this function will plot the data spectra onto the currently active plot axes. If one does not exist, it will generate a new figure. To ensure that a plot function creates a new figure, the `"NewFigure"` option can be set to `"Yes"`:

```matlab
handles = dataset.plot("NewFigure", "Yes");
```

This is functionally equivalent to calling `figure()` to create a new figure prior to plotting.

<p align="center">
   <strong>Example Data Plot</strong><br>
   <span>(two data spectra in dataset)</span><br>
   <img src="README FILES/ExampleDataPlot.png" width="400">
</p>


#### Data Plot Handles
The return value of the data plotting function `Dataset.plot()` is an array of `Line` handle objects that reference each of the plotted spectra. For example, the return value `handles` of the data plotting function called on a dataset containing three data spectra would be a 3x1 Line array where the handle to the first line plot can be accessed by simple indexing:

```matlab
handlePeak1 = handles(1);
```

#### Fit Plotting
Plotting peak fits is treated as a separate operation from plotting only data spectra and is handled through the `Dataset.PlotFit()` function:

```matlab
handles = dataset.plotFit();
```

As with data spectrum plotting, this function will plot the data spectra onto the currently active plot axes by default. If one does not exist, it will generate a new figure. To ensure that a plot function creates a new figure, the `"NewFigure"` option can be set to `"Yes"`:

```matlab
handles = dataset.plotFit("NewFigure", "Yes");
```

<p align="center">
   <strong>Example Peak Fit Plot</strong><br>
   <span>(one peak series, four individual peaks)</span><br>
   <img src="README FILES/ExampleFitPlot.png" width="400">
</p>

#### Fit Plot Handles
The return value of the fit plotting function `Dataset.plotFit()` is a cell array containing a structure for each plotted fit. For example, the return value `handles` of the fit plotting function called on a dataset containing two peak fits would be a 2x1 cell array where each cell contains a structure. Each structure in turn contains an organized collection of handles with the following fields: 

- base (base data spectrum being fitted) [1x1 Line]
- fit (total fit spectrum) [1x1 Line]
- peaks (individual peaks within fit) [Nx1 Line]

To access the `Line` handle for the base spectrum of the first plot fit, for example, the output variable `handles` can be accessed as

```matlab
handleBase = handles{1}.base;
```

Similarly, the `Line` handle for the first fit peak can be accessed as

```matlab
handlePeak1 = handles{1}.peaks(1);
```

#### Plot Customization
Appearance settings for spectrum `Line` plots can be easily customized by first accessing the handle to the desired `Line` plot, then calling the built-in MATLAB function `set()` with the desired appearance key-value pairs and passing the handle as the first argument. For example, the following call would change the line thickness and color of the first peak in a fit plot:

```matlab
set(handles{1}.peaks(1), "LineWidth", 1.5, "Color", "r")
```

To customize `Axes` and `Figure` settings, use the built-in MATLAB functions `gca` (get current axes) and `gcf` (get current figure) to retrieve their object handles after plot creation. Then, set the desired properties using the `set()` function as before. For example, the following call would make the background of the plot figure transparent:

```matlab
set(gcf, "Color", "None")
```

> Details on appearance settings available for customization can be found in the MATLAB documentation [here](https://www.mathworks.com/help/matlab/ref/plot.html) for line plots, [here](https://www.mathworks.com/help/matlab/ref/axes.html) for plot axes, and [here](https://www.mathworks.com/help/matlab/ref/figure.html) for figures.

### Saving Results
Peak fitting results can be saved to an Excel spreadsheet file using the function `Dataset.writeToXLSX()` and passing the function a file name:

```matlab
dataset = dataset.writeToXLSX("{file name}.xlsx");
```

Results for each peak fit in the dataset are recorded as separate sheets within the file, and each sheet contains the following fit information:
- Fit RMSD value
- Table of optimized peak parameters
   - Peak center
   - Peak width parameter (sigma)
   - Peak height
- Paired columns of X and Y data
   - Base data spectrum 
   - Optimized fit
   - Each component peak within fit

### Example Workflow
An example workflow that demonstrates data importing and cleaning, peak fitting, and writing fit results to file is provided along with the PeakPika package. This workflow can be found in the script "ExampleAnalysisScript.m" in the "Worked Example" folder of the package release. 

> The example workflow script is divided into sections by section markers `%%` that allow each section to be run individually by navigating to `Editor` > `Section` > `Run Section` in the top tools panel of the MATLAB IDE. It is recommended to run each section individually as running the script from start to finish will result in many plots being generated at once, potentially freezing the user's computer for a number of minutes.

## Reporting Bugs and Errors
Please report any bugs or errors encountered [here](https://github.com/AlexSimafranca/PeakPika/issues) or by emailing the package maintainer. 

## Authors
Creator and maintainer: Alexander Simafranca alesimaf@g.ucla.edu

