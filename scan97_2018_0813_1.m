% this scripts summarizes the parameters to analyze scan 97 from 201_0813_1
% file

% temperatures in C
TCV = [800]; % Equilibrium Thermocouple (new heater)

% power in Watts
PWV = [10];

% Filenames
specfilenameM = ['2018_0813_1'];

% scans number (see lab book)
SCNstrM = ['097']; 

% center of larger ROI in pixels
XCENV = [93];% del
YCENV = [138]; % nu

% Width of the larger ROI
XWIDV = [25];
% Pilatus detector on sevchex arm (X is del and Y is nu)
YWIDV = [45];


% Offset of two of the ROIS around the Positions of CTRs
ymax = [8];

% Min and max on time range for delta-time average
tminv = [500];
tmaxv = [4000];
