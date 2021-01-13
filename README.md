# Fast-forward firing pattern analysis

# FFFPA: Fast-forward firing pattern analysis

## Installation
1.	Installation requires MATLAB 2018b or newer with Signal Processing and Curve Fitting toolboxes 
2.	Download FFFPA plugin from https://github.com/zubara/fffpa and unzip it to a local directory.
3.	Double-click on ‘fffpa.mlappinstall’ file from MATLAB file explorer.
4.	Click “Install” in the opened dialog window.
5.	Once installed, open the app from MATLAB Apps toolstrip.


# Documentation
## Data import
1.	Click “Select files” and choose one or several ‘.abf’ files to process.
2.	Specify sampling frequency in Hz as well as the time of the onset and the offset of the injected current step in seconds 
3.	Choose sweeps that you wish to analyze. FFFPA can detects e.g. rheobase current step defined as the first sweep where APs occurred or saturation current step defined as sweep with maximum amount of APs.
4.	Specify summary statistics for analyzed APs. If multiple sweeps are analyzed these options can differ for the first and other sweeps.
5.	Specify filename and the desired format of the output.

## Adjust Event detection and feature extraction
1.	(Optional) Go to ‘Advanced’ tab to specify the features that you wish to extract from the traces.
2.	Click ‘Run’ button on the ‘Main’ tab.

## Inspection and manual correction of detection results
1.	To inspect the results of automatic feature detection, go to the ‘Data viewer’ tab. This tab opens automatically after clicking ‘Run’ once the data is processed.
2.	Plot area displays the raw trace as well as all detected events (AP peaks, blue circles, activation thresholds (red), after-hyperpolarization (AHP, cyan) and after-depolarization. Please note that by default ADP is only detected after the first AP or train of APs at rheobase.
3.	Use ‘Next’ button to go through all sweeps of all analyzed cells consecutively or use ‘cell ID’ and ‘sweep #’ fields to navigate to the specific sweep of the specific cell manually.
4.	If needed, Data viewer allows you to correct the detection results. To do it first click ‘Edit’, then click the element that you want to adjust and specify its new position by clicking on the trace. After all corrections are made press ‘X’ to save the results.
5.	You can also discard the current cell from the output by unchecking the ‘Save to output’ box.
6.	Once satisfied with the results click ‘Save selected’. This step will compute features based on the (adjusted) detections and save the output file to disk. The summary file can be found in the same directory as the input data.

# Please cite
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3667731.svg)](https://doi.org/10.5281/zenodo.3667731)





