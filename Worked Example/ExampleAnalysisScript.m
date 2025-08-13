%%%---VARIABLES---%%%
% Add working directory to MATLAB path (run in folder containing "+Pika" package folder)
addpath(genpath("."))

% UV-vis data paths
dataFolderPath = "./Worked Example/Example Data/Polymer1/";
backgroundFilePath = "./Worked Example/Example Data/20230424_baseline_H2O.txt";

dataPrelim = "mgml_";
dataPostlim = "N";


%%%---MAIN CODE---%%%
%%---Data Pre-Treatment---%%
%% Import the PeakPika package
import Pika.*   % the "+Pika" folder should already be on the MATLAB path


%% Read in data from path
dataset = Pika.Dataset();
dataset = dataset.readData(dataFolderPath, "HeaderLines"    , 2           , ...
                                           "NamePrelimiter" , dataPrelim  , ...
                                           "NamePostlimiter", dataPostlim);
dataset = dataset.readBackground(backgroundFilePath, "HeaderLines", 2);

h1 = dataset.plot("NewFigure", "Yes", "ShowBackground", "Yes");


%% Subtract baseline
datasetBkgd = dataset.subtractBackground();
h2 = datasetBkgd.plot("NewFigure", "Yes");


%% Normalize data
rangeOfMaxPeak = [500, 800];

datasetBkgdN = datasetBkgd.normalizeData(rangeOfMaxPeak);
h3 = datasetBkgdN.plot("NewFigure", "Yes");


%% Convert to eV
datasetEV = datasetBkgdN.convert2eV();
h4 = datasetEV.plot("NewFigure", "Yes");


%% Trim data
trimRange = [0.5, 3];

datasetEVTrim = datasetEV.trimData(trimRange);
h5 = datasetEVTrim.plot("NewFigure", "Yes");


%% Mask data
maskRange = [1.40, 1.42];

datasetEVMask = datasetEVTrim.maskData(maskRange);
h6 = datasetEVMask.plot("NewFigure", "Yes");

% ---------------------------------NOTES--------------------------------- %
% Now that our data is cleaned, let's do an initial peak fitting attempt.
% ----------------------------------------------------------------------- %


%%---Initial Peak Fitting Attempt---%%
%% Add peak series
numPeaks = 5;
firstPeakCenter = 1.9;
peakSpacing = 0.12;
peakSigma = 0.05;
peakHeightList = [0.45, 0.55, 0.45, 0.3, 0.25];

datasetFit1 = datasetEVMask.addPeakSeries(numPeaks, firstPeakCenter, peakSpacing, peakSigma, "PeakHeightList", peakHeightList);
h7 = datasetFit1.plotFit("NewFigure", "Yes");


%% Add individual peaks
peakCenterList = [0.3, 1.2, 1.6, 2.8];
peakSigmaList  = [0.3, 0.2, 0.2, 0.3];
peakHeightList = [1.0, 0.4, 0.3, 0.3];

datasetFit2 = datasetFit1.addPeak(peakCenterList, peakSigmaList, "PeakHeight", peakHeightList);
h8 = datasetFit2.plotFit("NewFigure", "Yes");


%% Run peak fit optimization
datasetFitOpt = datasetFit2.optimizeFit("NumIterations", 50);   % 50 cycles of parameter adjustments
h9 = datasetFitOpt.plotFit("NewFigure", "Yes");

% ---------------------------------NOTES--------------------------------- %
% This fit has converged, but it's not great based on the seemingly random 
% vibronic peak heights. It's clear the optimization fell into a local
% minimum. We should adjust our starting parameters and try again.
% ----------------------------------------------------------------------- %


%%---Peak Fitting Adjustments---%%
%% Re-add peak series
numPeaks = 5;
firstPeakCenter = 1.9;
peakSpacing = 0.18;   % increased peak spacing
peakSigma = 0.07;   % increased peak width
peakHeightList = [0.55, 0.55, 0.45, 0.3, 0.25];   % increased first peak height

datasetFit1 = datasetEVMask.addPeakSeries(numPeaks, firstPeakCenter, peakSpacing, peakSigma, "PeakHeightList", peakHeightList);
h10 = datasetFit1.plotFit("NewFigure", "Yes");


%% Re-add individual peaks
peakCenterList = [0.3, 1.2, 1.6, 2.8];   % individual peaks are the same as before
peakSigmaList  = [0.3, 0.2, 0.2, 0.3];
peakHeightList = [1.0, 0.4, 0.3, 0.3];

datasetFit2 = datasetFit1.addPeak(peakCenterList, peakSigmaList, "PeakHeight", peakHeightList);
h11 = datasetFit2.plotFit("NewFigure", "Yes");


%% Run peak fit optimization
datasetFitOpt = datasetFit2.optimizeFit("NumIterations" , 50    , ...   
                                        "PeakCenterStep", 0.005 , ...   % decrease peak center step size (default value of 0.01)
                                        "PeakHeightStep", 0.005);       % decrease peak height step size (default value of 0.01)
h12 = datasetFitOpt.plotFit("NewFigure", "Yes");

% ---------------------------------NOTES--------------------------------- %
% This is a better fit. The vibronic peaks have a smooth trend in height,
% and there is a tighter agreement with the data on the right side of the
% bandgap vibronic cluster. 
% 
% What would happend if we allowed each vibronic peak to have its own 
% width?
% ----------------------------------------------------------------------- %


%%---Trying Alternate Fits---%%
%% Alternate fit 1
datasetFitAlt = datasetFit2.optimizeFit("NumIterations" , 50   , ...   
                                        "PeakCenterStep", 0.005, ... 
                                        "PeakHeightStep", 0.005, ...
                                        "SharePeakWidth", "No");   % disable peak width sharing for the peak series
h13 = datasetFitAlt.plotFit("NewFigure", "Yes");

% ---------------------------------NOTES--------------------------------- %
% This vibronic peaks are now quite messy. Note how the 4th and 5th
% vibronic peaks are randomly much wider than the others. This is a worse
% fit than before. 
%
% What about if the vibronic peaks no longer had a common peak spacing?
% ----------------------------------------------------------------------- %


%% Alternate fit 2
datasetFitAlt = datasetFit2.optimizeFit("NumIterations"   , 50   , ...   
                                        "PeakCenterStep"  , 0.005, ... 
                                        "PeakHeightStep"  , 0.005, ...
                                        "SharePeakSpacing", "No");   % disable peak spacing sharing for the peak series   
h14 = datasetFitOpt.plotFit("NewFigure", "Yes");

% ---------------------------------NOTES--------------------------------- %
% Also a worse fit than before. Looks like we'll go with the fit we had
% before with both peak width and peak sharing enabled for the peak series
% (they are enabled, or "Yes", by default).
%
% It may not always be the case that width and spacing sharing are optimal,
% so feel free to experiment with different fitting parameters to see what
% works best for your dataset.
% ----------------------------------------------------------------------- %


%%---Print Fit Results to File---%%
%% Print peak parameters and spectra to file
datasetFitOpt.writeToXLSX("test.xlsx");   % work-in-progress

