% Two-time analysis of March 2018 TiO2 data
% Read and analyze dataset(s)
% March 17, 2018 GBS
% just use equil 2-time, average over a time range

 %startup; % when the detector is Pilatus (March 2018 beamtime)
 startup_pilatus4;


% General flag to decide what to do:
%0 means read again the data; 
%1 means read previously saved data;

skip = 0; 


% General flag to initialize the parameters for the 'equilibrium' or the
% 'growth' or 'power_series' data:
XPCS_param_file = 'XPCS_param_TTsput_2018_08_newdata';  %
XPCS_param_file = 'XPCS_param_TTsput_2018_08_concise_newdata';

[plotsmooth,plotall,plotorig,plotnew,ixplotmin,ixplotmax,...
    pcolorCC2,plothalf,plotfull,plotcorrorig,plotcorrnew,...
    plotcontrast,pcolorextra,posflag,fitorig,fitnew, altnorm,plotdeltat]...
    = XPCS_initialize_parameters.TTsput_flags();

[DOCU0,DOCU1,LOGFLAG,delay,XFRAC] = XPCS_initialize_parameters.TTparameters_allruns();

[TCV,nT,iiT,specfilenameM,SCNstrM,XCENV,YCENV,XWIDV,YWIDV,CROPV,...
    ymax,tminv,tmaxv,clrs, DXImin, ...
    DXImax,NRY,PWV,flag_equil_or_growth,pilatus_flag] = XPCS_initialize_parameters.TTparameters_singlerun(XPCS_param_file);

[fitrange_time_iiT,pin_iiT,dp_iiT,qfitrange ] ...
    = XPCS_initialize_parameters.TTparameters_fit_singlerun(XPCS_param_file);


[hwttr_allT,hwttc_allT,wrq_allT,wcq_allT,offsetcc_allT,offsetrc_allT,...
    tbin_allT,CWID_allT,N_degree_allT] =  XPCS_initialize_parameters.TTparameters_2timecorr_calc(XPCS_param_file);

[POSITION,PAPERPOSITION,FONTSIZE,CMAX,CLIM,XCOLlabel,YROWlabel,...
                AXISdet,DOCUclim,INFOstr] = XPCS_initialize_parameters.TTplot_parameters();


 [ImageJ,Xstepcol,SINGLE,BKG,...
 scanflag,imname,p_image,ending,POINTSUMS] = XPCS_initialize_parameters.TTsput_read_ini();

[D_ds,kvector,lambda,pixel_size,th_Bragg] = XPCS_initialize_parameters.TT_experiment();

% set the scans that you want to read
iTV = 1;%[1:numel(iiT)];

%_._._._._._._._._._._._._._._._._._._._._._._._._.._._._._._._._._._.._.._._._______________________
if ~skip
    disp('Recalculating 2-time from experimental datasets.');
    
    for iT = 1:numel(iTV)
        clear Read_Allscans;
        % read data
        
         warning('off','all');
        [Read_Allscans(iT).IIstruct] = XPCS_read_data.TTsput_read(iiT(iT),TCV,specfilenameM,SCNstrM,DOCU0,DOCU1,ImageJ,Xstepcol,BKG,scanflag,imname,p_image,ending,POINTSUMS,pilatus_flag,CROPV); % reads a dataset, makes IInormb
        warning('on','all');
        
        %%{
        %Save the data
        Singlescan_struct = Read_Allscans(iT);
        
        switch flag_equil_or_growth
            
            case 'equilibrium'
                  Singlescan_namefile = ['Scan' SCNstrM(iT,:) '_' specfilenameM(iT,:) '.mat'];
            case 'growth'
                  Singlescan_namefile = ['Scan_growth' SCNstrM(iT,:) '_' specfilenameM(iT,:) '.mat'];
                  
            case 'power_series'
                  Singlescan_namefile = ['Scan_growth' SCNstrM(iT,:) '_' specfilenameM(iT,:) '.mat'];

        end
      
        save(Singlescan_namefile,'Singlescan_struct','-v7.3');
        
        
        %}
        
        
    end
             
end

for iT = iTV
        
         SCNstrM_numer = str2num(SCNstrM(iT,:)) ;  
         switch flag_equil_or_growth
             
             case 'equilibrium'
                 Singlescan_namefile = ['Scan' SCNstrM(iT,:) '_' specfilenameM(iT,:) '.mat'];
             case 'growth'
                 Singlescan_namefile = ['Scan_growth' SCNstrM(iT,:) '_' specfilenameM(iT,:) '.mat'];
                 
             case 'power_series'
                Singlescan_namefile = ['Scan_growth' SCNstrM(iT,:) '_' specfilenameM(iT,:) '.mat'];
         end
        load(Singlescan_namefile);
        Read_Allscans(iT).IIstruct = Singlescan_struct.IIstruct;
    end

      

for iT = iTV
    % Initialize ROIs:
    [Allscans(iT).ROIS_struct] = XPCS_read_data.TTsput_prepare_ROIS(iiT(iT),XCENV,YCENV,XWIDV,YWIDV,ymax);
    
    
    % Calculate Q-range
     [Allscans(iT).Qval_struct] = XPCS_analysis.calculate_qval(XCENV(iT),YCENV(iT),[1:Read_Allscans(iT).IIstruct.Nc],[1:Read_Allscans(iT).IIstruct.Nr],D_ds,kvector,pixel_size,th_Bragg);

    
    
    %Initialize the range of the 2 times correlation functions:
    [Allscans(iT).itt_range_struct] = XPCS_analysis.prepare_subrange_for2corr(iT,Read_Allscans(iT).IIstruct.Xamount,XCENV,YCENV,XWIDV,YWIDV,tminv,tmaxv);
    
    % Bin IInorm
    [Allscans(iT).IIbin_struct] = XPCS_analysis.bin_scans_in_time(Read_Allscans(iT),Allscans(iT),iT,tbin_allT,ImageJ,POINTSUMS);
    
    % Calculate 2 time correlaiont functions
    flag_mean_or_poly = 'poly';
    [Allscans(iT).CCN2_struct,Allscans(iT).IIbin_struct] = XPCS_analysis.calc_2time_corr(Allscans(iT).IIbin_struct,iT,N_degree_allT,flag_mean_or_poly);
    
    [Allscans(iT).CCN2avg_struct] = XPCS_analysis.from_CCN2V_to_CCN2avg(Allscans(iT),iT,hwttr_allT,hwttc_allT,wrq_allT,wcq_allT,offsetcc_allT,offsetrc_allT,D_ds,kvector,pixel_size,th_Bragg);
    
    flag_row_col = 'col';
    num_col_or_row = 1;
    
    %DisplayFunctions_XPCS.display_IInormbbref(Allscans(iT).IIbin_struct,Allscans(iT).CCN2avg_struct.boxcenterrc,50+iT,AXISdet,D_ds,kvector,pixel_size,th_Bragg);

     % plot the CTR images in pixels:
    DisplayFunctions_XPCS.display_CCN2avg(Allscans(iT).CCN2avg_struct,num_col_or_row,flag_row_col,40+iT+num_col_or_row,AXISdet);
    QvalFlag = 0;
    fignum = 90;
    DisplayFunctions_XPCS.display_IInormb(Allscans(iT).IIbin_struct,Allscans(iT).IIbin_struct.IInormbb,'Log 10 binned IInormb',QvalFlag,fignum+iT,ImageJ,SINGLE,XCOLlabel,YROWlabel,AXISdet,INFOstr,D_ds,kvector,pixel_size,th_Bragg);
    DisplayFunctions_XPCS.display_grid_CCN2avg(Allscans(iT).CCN2avg_struct.ittccen,Allscans(iT).CCN2avg_struct.ittrcen,Allscans(iT).CCN2avg_struct.Ncq_Nrq(1),Allscans(iT).CCN2avg_struct.Ncq_Nrq(2),QvalFlag,iT,hwttr_allT,hwttc_allT,wrq_allT,wcq_allT,offsetcc_allT,offsetrc_allT,fignum+iT,ImageJ,D_ds,kvector,pixel_size,th_Bragg);
    
    % plot the CTR images in reciprocal space units
    QvalFlag = 1;
    fignum = 80;
    DisplayFunctions_XPCS.display_IInormb(Allscans(iT).IIbin_struct,Allscans(iT).IIbin_struct.IInormbb,'Log 10 binned IInormb',QvalFlag,fignum+iT,ImageJ,SINGLE,XCOLlabel,YROWlabel,AXISdet,INFOstr,D_ds,kvector,pixel_size,th_Bragg);
    DisplayFunctions_XPCS.display_grid_CCN2avg(Allscans(iT).CCN2avg_struct.ittccen,Allscans(iT).CCN2avg_struct.ittrcen,Allscans(iT).CCN2avg_struct.Ncq_Nrq(1),Allscans(iT).CCN2avg_struct.Ncq_Nrq(2),QvalFlag,iT,hwttr_allT,hwttc_allT,wrq_allT,wcq_allT,offsetcc_allT,offsetrc_allT,fignum+iT,ImageJ,D_ds,kvector,pixel_size,th_Bragg);
    
    
    [Allscans(iT).CCN2S_struct] = XPCS_analysis.from_CCN2avg_to_CCN2S(Allscans(iT).CCN2avg_struct);
    
    [Allscans(iT).CCN2S_struct] = XPCS_analysis.fit_CCN2S_with_leasqr(Allscans(iT).CCN2S_struct,iT,pin_iiT,dp_iiT,'FittingFunctions.CCN2single_fit',[2:1:round(fitrange_time_iiT(iT)*Allscans(iT).IIbin_struct.Ntb)],num_col_or_row,flag_row_col,200);
    
    figh = 400+iT+num_col_or_row;
    [Allscans_fit(iT).pout,Allscans_fit(iT).sigma] = DisplayFunctions_XPCS.display_fit_result(Allscans(iT).CCN2S_struct,num_col_or_row,flag_row_col,figh);
    
    DisplayFunctions_XPCS.display_CCN2S(Allscans(iT).CCN2S_struct,num_col_or_row,flag_row_col,101+iT,[2:Allscans(iT).IIbin_struct.Ntb]);
    
end


%disp('Fit the data');

 % Fit the taus at different temperatures 
pp_index = 3;
[pout_struct] = DisplayFunctions_XPCS.display_fit_result_log(Allscans(iTV),pp_index,num_col_or_row,flag_row_col ,601);

save(['fit_Scan_growth' SCNstrM(iT,:) '.mat'],'pout_struct');

% fit pout(3,:)= tau to a power law 1/x

for iT = 1:numel(pout_struct)
    figure;
   
    pout_struct(iT).pout_fit = XPCS_analysis.fit_tau_with_leasqr(pout_struct(iT),'FittingFunctions.powerlaw_tau',qfitrange);
      
    
    figure(601);
    subplot(1,numel(pout_struct),iT);
    %title(['#' SCNstrM(iT,:)]);
    hold on;
    plot(pout_struct(iT).pout_fit.pout.x,pout_struct(iT).pout_fit.pout.fitfunc,'r','LineWidth',3.0);
    
    %figure(602);
    %hold on;
    %plot(pout_struct(iT).pout_fit.pout.x,pout_struct(iT).pout_fit.pout.fitfunc,'r','LineWidth',3.0);

end

% display the velocity vs temperature:
island_size = 2e3; % in Angstroms
[vel_struct] = DisplayFunctions_XPCS.display_vel_vs_temp(TCV,pout_struct,'Tau vs temperature for equilibrium',island_size,700);


island_size = 2e3; % in Angstroms
%[vel_struct] = DisplayFunctions_XPCS.display_vel_vs_temp(PWV,pout_struct,'Tau vs power',island_size,800);

DisplayFunctions_XPCS.display_growth_vs_power(PWV(iT),vel_struct,[1:numel(PWV(iT))], '900 C',1000);


% display the velocity vs the power
%power_vect = [10 20 50];
%DisplayFunctions_XPCS.display_growth_vs_power(power_vect,vel_struct,[5 6 7],'900 C',800);

%{
for iT = 1%iTV
    TTM_plot1; 
 end
%}



