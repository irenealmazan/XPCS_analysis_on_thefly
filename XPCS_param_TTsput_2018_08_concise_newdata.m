%%%% This scripts initializes the data file name, the scans numbers, the
%%%% temperatures, the power at which we measure for new data taken in
%%%% August 2018


%% Section 0: 

% flag:
flag_equil_or_growth = 'equilibrium'; % choose among 'equilibrium', 'growth' or 'power_series'

pilatus_flag = 'pilatus4';% choose between 'pilatus' (beamtime March 2018) and 'pilatus4'




% Parameters initialized in
% XPCS_initialize_parameters.TTparameters_singlerun
%XCENV = [135];YCENV = [91];XWIDV = [25];YWIDV = [90];CROPV = [1 194 1 200];
XCENV = [91];YCENV = [135];XWIDV = [90];YWIDV = [25];CROPV = [1 194 1 200];
specfilenameM = ['2018_0810_1']; 
flag_equil_or_growth = 'equilibrium'; % choose among 'equilibrium', 'growth' or 'power_series'
		TCV = [800]; PWV = [0];SCNstrM = ['019']; tminv = [2000];tmaxv = [5000]; DOCU = '[???] not printed';flagrotate = [1];
%{
        TCV = [800]; PWV = [0];SCNstrM = ['037']; tminv = [1];tmaxv = [5000]; DOCU = '10um s6hgap 3%O2 30mT shutter closed';flagrotate = [1];
		TCV = [800]; PWV = [0];SCNstrM = ['038']; tminv = [1];tmaxv = [2500]; DOCU = '20um s6hgap 3%O2 30mT shutter closed; oscillations in later';flagrotate = [1];

flag_equil_or_growth = 'equilibrium'; %?? should it be growth? and what other parameters need tobe set?c
		TCV = [800]; PWV = [15];SCNstrM = ['032']; tminv = [1000];tmaxv = [5000]; DOCU = '10um s6hgap 3%O2 30mT';flagrotate = [1];
		TCV = [800]; PWV = [0];SCNstrM = ['019']; tminv = [2000];tmaxv = [5000]; DOCU = '[???] not printed';flagrotate = [1];
        TCV = [400]; PWV = [0];SCNstrM = ['063']; tminv = [1];tmaxv = [2000]; DOCU = '400 C, 10 microns hslit';flagrotate = [1];
        TCV = [400]; PWV = [0];SCNstrM = ['064']; tminv = [1];tmaxv = [2000]; DOCU = '400 C, 20 microns hslit';flagrotate = [1];
        TCV = [400]; PWV = [0];SCNstrM = ['065']; tminv = [1];tmaxv = [2000]; DOCU = '400 C, 20 microns hslit';flagrotate = [1];
		TCV = [400]; PWV = [0];SCNstrM = ['067']; tminv = [1];tmaxv = [5000]; DOCU = '10um s6hgap 3%O2 30mT shutter closed';flagrotate = [1];

% set of data to investigate the constrast:
XCENV = [91];YCENV = [135];XWIDV = [90];YWIDV = [25];CROPV = [1 194 1 200]
specfilenameM = ['2018_0810_1']; 
flag_equil_or_growth = 'equilibrium';
        TCV = [700]; PWV = [0];SCNstrM = ['078']; tminv = [1];tmaxv = [1000]; DOCU = '10um s6hgap 3%O2 30mT shutter closed';flagrotate = [1];
        TCV = [700]; PWV = [0];SCNstrM = ['079']; tminv = [1];tmaxv = [2200]; DOCU = '20um s6hgap 3%O2 30mT shutter closed';flagrotate = [1];
        TCV = [700]; PWV = [0];SCNstrM = ['080']; tminv = [1];tmaxv = [2000]; DOCU = '20um s6hgap 3%O2 30mT shutter closed';flagrotate = [1];
flag_equil_or_growth = 'growth'                  
TCV = [700]; PWV = [0];SCNstrM = ['081']; tminv = [1000];tmaxv = [2976]; DOCU = '20um s6hga 700 C growth 10 W';flagrotate = [1];
        
%{
        TCV = [750]; PWV = [0];SCNstrM = ['085']; tminv = [1];tmaxv = [2000]; DOCU = '20um s6hgap 750 C shutter close and rf off';flagrotate = [1];
        TCV = [750]; PWV = [0];SCNstrM = ['087']; tminv = [1];tmaxv = [2200]; DOCU = '20um s6hgap 750 c shutter close and rf on';flagrotate = [1];

        
flag_equil_or_growth = 'growth'        
        TCV = [750]; PWV = [0];SCNstrM = ['088']; tminv = [500];tmaxv = [2200]; DOCU = '20um s6hgap 750 C growth 10 W';flagrotate = [1];

flag_equil_or_growth = 'equilibrium';
        TCV = [800]; PWV = [0];SCNstrM = ['104']; tminv = [1];tmaxv = [2100]; DOCU = '20um s6hgap 800 C shutter closed and rf on';flagrotate = [1];

flag_equil_or_growth = 'growth'        
        TCV = [800]; PWV = [0];SCNstrM = ['106']; tminv = [500];tmaxv = [2200]; DOCU = '20um s6hgap 800 C growth 10 W';flagrotate = [1];
  %}    

XCENV = [91];YCENV = [135];XWIDV = [90];YWIDV = [25];CROPV = [1 194 1 200]       
specfilenameM = ['2018_0811_1'];   
        TCV = [850]; PWV = [0];SCNstrM = ['009']; tminv = [1];tmaxv = [10]; DOCU = '850 C test';flagrotate = [1];
        TCV = [850]; PWV = [0];SCNstrM = ['010']; tminv = [1];tmaxv = [4900]; DOCU = '20 microns shgap 850 C shutter close, rf off';flagrotate = [1];
        TCV = [850]; PWV = [0];SCNstrM = ['013']; tminv = [700];tmaxv = [2800]; DOCU = '20 microns shgap 850 C growth at 10 W';flagrotate = [1];
%}
        
%% Section 1
%{
% temperatures in C
%TCV = [800]; % Equilibrium Thermocouple (new heater)

% power in Watts
%PWV = [0];

% Filenames
%specfilenameM = ['2018_0809_5'];
%specfilenameM = ['2018_0810_1'];

% name of the tifs for reading the data
% imname = []; %specfilenameM;   % ???? CT if imname is pased to 


% scans number (see lab book)
%SCNstrM = ['015']; % Equilibrium
%SCNstrM = ['019']; % Equilibrium

% Center pixels in between the CTRS
%XCENV = [128];% del
%YCENV = [98]; % nu
%XCENV = [135];% del
%YCENV = [91]; % nu


% Width of the larger ROI
%XWIDV = [25];
% Pilatus detector on sevchex arm (X is del and Y is nu)
%YWIDV = [90];

% Area of the detector you want to consider for the analysis  
% ? CT is this 'absolute' or relative?
%CROPV = [1 194 1 200];



% Min and max on time range for delta-time average in sec
%tminv = [2000];
%tmaxv = [5000];

%}
% Positions of CTRs for ROIS   % CT note - this seems to be needed, whatever itis
ymax = [8];    % CT note - this seems to be needed, whatever itis

%% Section 2

 %%%%%% Set of parameters to calculate the area where the 2 times correlation
 %%%%%function is calculated:

hwttr_allT = [1]; % row half width (pixels)=> box of 2*hwttr+1 pixels
hwttc_allT = [32]; % col half width (pixels)

wrq_allT = [7];
wcq_allT = [0];

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
 
 % parameters to initialize the 2-D Savitzky-Golay smoothing filter
% For symmetric di or dj, use even pdi, pdj; next higher odd gives same
% answer
% these values should be moved to outer program
% need to optimize by overplotting smoothed function and data
maxpd = 2; % determines maximum degree of smoothing polynomial
iii = 5; % half-width of pixel range in del
jjj = 5; % half-width of time steps


 %% Section 3


%%%%%%% fit parameters and range for the single exponential decays (CC2Ns) and the time constant vs
% q

fitrange_time_iiT = [2/5]; % fraction of the total time range for the fit of the CC2NS
pin_iiT = [0 1e-2 50 0]; % initial parameters for the fit of the CC2NS  to a single exponential decay
dp_iiT =  [[1 1 1 0]*0.0001]; % tolerance on the fitted parameters of the CC2NS
qfitrange = [2.5e-5 1.5e-2]; % Q range for the fit of the time constants
   
   
   
   

















