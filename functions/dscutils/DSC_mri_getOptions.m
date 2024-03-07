function [optionsOUT] = DSC_mri_getOptions(app)

% DISPLAY OPTIONS
optionsOUT.display = 0; % 0: off, 1: notify (text), 2: notify (images), 3: debug
optionsOUT.waitbar = 0; % 0: off, 1: on

% DATA PREPARATION OPTIONS
optionsOUT.mask.npixel = 50;
% Represents the minimum number of pixels in a connected component used as a threshold to exclude the scalp and areas adjacent to the brain in the image.

optionsOUT.conc = 0;
% 0: the provided data is signal, 1: the provided data is concentrations

optionsOUT.S0.nSamplesMin = 6;
% Minimum number of initial scans used to calculate S0

optionsOUT.S0.nSamplesMax = 12;
% Maximum number of initial scans used to calculate S0

optionsOUT.S0.thresh = 0.2;
% Threshold used to determine the appearance time of the contrast agent

% ARTERIAL INPUT FUNCTION (AIF) IDENTIFICATION OPTIONS
optionsOUT.aif.enable = 1;
% 0: do not calculate AIF, 1: calculate AIF

optionsOUT.aif.ricircolo = 1;
% 0: do not consider recirculation, 1: fit recirculation

optionsOUT.aif.nSlice = 0;
% Slice on which to search for the AIF (0: allows the operator to select the slice)

optionsOUT.aif.semiasseMaggiore = app.LongAxisEditField.Value/100;  % 0.15
% Size of the major axis for the search area

optionsOUT.aif.semiasseMinore = app.ShortAxisEditField.Value/100;     % 0.15;
% Size of the minor axis for the search area

optionsOUT.aif.pArea = 0.3000;
% Percentage of voxels discarded due to AUC

optionsOUT.aif.pTTP = 0.3000;
% Percentage of voxels discarded due to TTP

optionsOUT.aif.pReg = 0.0500;
% Percentage of voxels discarded due to regularity of the curve

optionsOUT.aif.diffPicco = 0.0400;
% Threshold to decide whether to select the cluster based on peak or TTP

optionsOUT.aif.nVoxelMax = app.MaxAIFpixelsEditField.Value;   % 10
% Maximum number of voxels selected for the AIF

optionsOUT.aif.nVoxelMin = 4;
% Minimum number of voxels selected for the AIF

% Correction of the concentration calculation formula from the signal in the case of AIF calculation
optionsOUT.qr.enable = 0; % 0: do not apply the correction, 1: apply the correction
optionsOUT.qr.b = 5.7400e-004;
optionsOUT.qr.a = 0.0076;
optionsOUT.qr.r = 0.0440;

% DECONVOLUTION METHODS OPTIONS
optionsOUT.deconv.SVD.threshold = 0.2; % SVD threshold
optionsOUT.deconv.SVD.residual = 1; % 0: do not save residuals, 1: save residuals

optionsOUT.deconv.cSVD.threshold = 0.1; % cSVD threshold (0.1 for data obtained at 1.5T)
optionsOUT.deconv.cSVD.residual = 1; % 0: do not save residuals, 1: save residuals

optionsOUT.deconv.oSVD.OIthres = 0.035; % 10% threshold as in Ostergaard and Calamante
optionsOUT.deconv.oSVD.OIcounter = 1;
optionsOUT.deconv.oSVD.residual = 1; % 0: do not save residuals, 1: save residuals

% ADD PARAMETERS FOR STABLE SPLINE HERE
optionsOUT.deconv.SS.residual = 1;

% Methods to apply for perfusion calculation
if app.SVDButton.Value == 1
    optionsOUT.deconv.method = {'SVD'};
end

if app.cSVDButton.Value == 1
    optionsOUT.deconv.method = {'cSVD'};
end

if app.oSVDButton.Value == 1
    optionsOUT.deconv.method = {'oSVD'};
end


% PROPORTIONALITY CONSTANTS
optionsOUT.par.kh = 1;
optionsOUT.par.rho = 1;
optionsOUT.par.kvoi = 1;


% Progress bar settings
optionsOUT.step_cvb = 1;
optionsOUT.step_cbv_lc = 1;
optionsOUT.step_cbf = 1;

