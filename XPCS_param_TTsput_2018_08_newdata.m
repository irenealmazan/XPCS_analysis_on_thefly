%%%% This scripts initializes the data file name, the scans numbers, the
%%%% temperatures, the power at which we measure for new data taken in
%%%% August 2018


%% Section 0: 

% flag:
flag_equil_or_growth = 'equilibrium'; % choose among 'equilibrium', 'growth' or 'power_series'

pilatus_flag = 'pilatus4';% choose between 'pilatus' (beamtime March 2018) and 'pilatus4'





%% Section 1

% temperatures in C
TCV = [800]; % Equilibrium Thermocouple (new heater)

% power in Watts
PWV = [0];

% Filenames
%specfilenameM = ['2018_0809_5'];
specfilenameM = ['2018_0810_1'];

% name of the tifs for reading the data
% imname = []; %specfilenameM;   % ???? CT if imname is pased to 


% scans number (see lab book)
%SCNstrM = ['015']; % Equilibrium
SCNstrM = ['019']; % Equilibrium

% Center pixels in between the CTRS
XCENV = [128];% del
YCENV = [98]; % nu
XCENV = [135];% del
YCENV = [91]; % nu


% Width of the larger ROI
XWIDV = [25];
% Pilatus detector on sevchex arm (X is del and Y is nu)
YWIDV = [90];

% Area of the detector you want to consider for the analysis
CROPV = [1 194 1 200];

% Positions of CTRs for ROIS
ymax = [8];

% Min and max on time range for delta-time average in sec
tminv = [2000];
tmaxv = [5000];


%% Section 2

 %%%%%% Set of parameters to calculate the area where the 2 times correlation
            % function is calculated:

hwttr_allT = [32]; % row half width (pixels)=> box of 2*hwttr+1 pixels
hwttc_allT = [1]; % col half width (pixels)

wrq_allT = [0];
wcq_allT = [7];

% tuning the center of the reciprocal space for the CC2avg
% analysis (it can be larger than Ncs and Nrs respectively
offsetcc_allT = [0];
offsetrc_allT = [0];

% number of scans to bin together for 2-time calcs
tbin_allT = [10];


CWID_allT = [0.3]; % Parameter for integer/half-integer ML integration

 % degree of the polynomial when we fit the IInormbb reference
 % in the new way of analysis
 N_degree_allT = [1 ];


 %% Section 3


%%%%%%% fit parameters and range for the single exponential decays (CC2Ns) and the time constant vs
% q

fitrange_time_iiT = [2/5];
pin_iiT = [0 1e-2 50 0];
dp_iiT =  [[1 1 1 0]*0.0001];
qfitrange = [2.5e-5 1.5e-2]; %in nu or del depending the orientation of the detector (back wall with vs arm)
   
   
   
   
















