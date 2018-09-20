
classdef XPCS_analysis
    % This library contains all the functions which allow us to analyze the
    % time correlatin functions
    properties(Constant)
    end
    
    
    methods(Static)
  
         function Qval_struct = calculate_qval(XCEN_II,YCEN_II,Ncq_vect,Nrq_vect,D_ds,kvector,pixel_size,th_Bragg)
            

            Qval_struct.nu = (kvector*(Nrq_vect-YCEN_II )*pixel_size/D_ds)*1e-10; % nu direction, in 1/Angstroms            
            Qval_struct.del = (kvector*(Ncq_vect-XCEN_II)*pixel_size/D_ds)*sind(th_Bragg)*1e-10; % del direction  in 1/Angstroms
           
        end
        
        
        
        function [itt_range_struct] = prepare_subrange_for2corr(iT,Xamount,XCENV,YCENV,XWIDV,YWIDV,tminv,tmaxv)
            % subrange for 2-time correlation function calculation
            
            itt_range_struct.xx = [-XWIDV(iT):XWIDV(iT)]; % Pixel locations relative to center 
            itt_range_struct.yy = [-YWIDV(iT):YWIDV(iT)]; % Pixel locations relative to center
            itt_range_struct.ixtt = XCENV(iT) + itt_range_struct.xx;
            itt_range_struct.iytt = YCENV(iT) + itt_range_struct.yy;
            itt_range_struct.Ncs = length(itt_range_struct.xx);
            itt_range_struct.Nrs = length(itt_range_struct.yy);
            itt_range_struct.ittt = find(Xamount > tminv(iT) & Xamount < tmaxv(iT));
            itt_range_struct.Nts = length(itt_range_struct.ittt);
            
            
        end
        
      
        function [IIbin_struct] = bin_scans_in_time(Read_Singlescan_struct,Singlescan_struct,iT,tbinV,ImageJ,POINTSUMS)
            
            % Use sub-range of pixels and time
            
            iytt = Singlescan_struct.itt_range_struct.iytt;
            ixtt = Singlescan_struct.itt_range_struct.ixtt;
            ittt = Singlescan_struct.itt_range_struct.ittt;
            
          
            
            IInormbs = Read_Singlescan_struct.IIstruct.IInormb(iytt,ixtt,ittt);
            
            Nrs = size(IInormbs,1);
            Ncs = size(IInormbs,2);
            Nts = size(IInormbs,3);
            
            timestampXs = Read_Singlescan_struct.IIstruct.timestampX(ittt);
            Xamounts = Read_Singlescan_struct.IIstruct.Xamount(ittt);
            Xsteps = Read_Singlescan_struct.IIstruct.Xsteps(ittt);
            
            tbin = tbinV(iT);
            
            % First bin scans in time
            if tbin > 1
                tbin = floor(tbin);
                Ntb = floor(Nts/tbin);
                IIbin_struct.Ntb =  Ntb;
                IInormbe = reshape(IInormbs(:,:,1:Ntb*tbin),Nrs,Ncs,tbin,Ntb);
                IIbin_struct.IInormbb = squeeze(sum(IInormbe,3));
                Xamounte = reshape(Xamounts(1:Ntb*tbin),tbin,Ntb);
                Xstepse = reshape(Xsteps(1:Ntb*tbin),tbin,Ntb);
                IIbin_struct.Xamountb = squeeze(mean(Xamounte,1));
                IIbin_struct.Xstepsb = squeeze(mean(Xstepse,1));
                timestampe = reshape(timestampXs(1:Ntb*tbin),tbin,Ntb);
                IIbin_struct.timeb = squeeze(mean(timestampe,1));
                
                if ~isempty(POINTSUMS)
                    POINTSUMSB(:,1) = floor((POINTSUMS(:,1)-1+ImageJ)/tbin)+1-ImageJ;
                    POINTSUMSB(:,2) = floor((POINTSUMS(:,2)+ImageJ)/tbin) - ImageJ;
                else
                    POINTSUMSB = POINTSUMS;
                end
                
            else
                Ntb = Nts;
                IIbin_struct.Ntb = Nts;
                IIbin_struct.IInormbb = IInormbs;
                IIbin_struct.Xamountb = Xamounts';
                IIbin_struct.Xstepsb = Xsteps;
                IIbin_struct.timeb = timestampXs;
                POINTSUMSB = POINTSUMS;
            end
            
            if isempty(POINTSUMSB)
                IIbin_struct.POINTSB=[1 IIbin_struct.Ntb]-ImageJ;
            else
                IIbin_struct.POINTSB = POINTSUMSB;
            end
                        
            IIbin_struct.YROWpts = [1:Nrs] - ImageJ;
            IIbin_struct.XCOLpts = [1:Ncs] - ImageJ;
            IIbin_struct.TITLEstuct = Read_Singlescan_struct.IIstruct.TITLEstuct;
        end
        
        function [CCN2_struct,IIbin_struct] = calc_2time_corr(IIbin_struct,iT,N_degree,flag_mean_or_poly)
                     
            IInormbb =  IIbin_struct.IInormbb;
            
            Nrs = size(IInormbb,1);
            Ncs = size(IInormbb,2);
            Ntb = size(IInormbb,3);
            
            switch flag_mean_or_poly
                case 'mean'
                    IInormbb_ref = mean(IInormbb,3);
                case 'poly'
                    IInormbb_ref = FittingFunctions.fit_IInormb_intime(IIbin_struct, N_degree(iT));
                    IIbin_struct.N_degree = N_degree(iT);
            end
            
            IIbin_struct.IInormbb_ref = IInormbb_ref;
            
            
            dI = IInormbb - IInormbb_ref;
            dlnI = IInormbb./IInormbb_ref - 1;
             
            % Calc 2-time using time average mean, no ensemble
            
            %CCN2_struct.IIM2 = NaN*ones(Nrs,Ncs,Ntb,Ntb);
            CCN2_struct.CCN2 = NaN*ones(Nrs,Ncs,Ntb,Ntb);
            
            for ii = 1:Ntb
                for jj = 1:ii
                    %CCN2_struct.IIM2(:,:,ii,jj) = dI(:,:,ii).*dI(:,:,jj);
                    %CCN2_struct.IIM2(:,:,jj,ii) = CCN2_struct.IIM2(:,:,ii,jj);
                   
                    CCN2_struct.CCN2(:,:,ii,jj) = dlnI(:,:,ii).*dlnI(:,:,jj);
                    CCN2_struct.CCN2(:,:,jj,ii) = CCN2_struct.CCN2(:,:,ii,jj);
                end
            end
            %IID2 = diag(IIM2);
            %CC2 = IIM2./sqrt(IID2*IID2'); % Normalized to make diagonal unity
            CCN2_struct.TITLEstuct = IIbin_struct.TITLEstuct;
        end
        
        function CCN2avg_struct = from_CCN2V_to_CCN2avg(Single_scan_struct,iT,hwttr_allT,hwttc_allT,wrq_allT,wcq_allT,offsetcc_allT,offsetrc_allT,D_ds,kvector,pixel_size,th_Bragg)
            
            hwttr = hwttr_allT(iT);
            hwttc = hwttc_allT(iT);
            wrq = wrq_allT(iT);
            wcq = wcq_allT(iT);
            offsetcc = offsetcc_allT(iT);
            offsetrc = offsetrc_allT(iT);
            
            
            Nts = size(Single_scan_struct.CCN2_struct.CCN2,3);
            Ncs = size(Single_scan_struct.CCN2_struct.CCN2,2);
            Nrs = size(Single_scan_struct.CCN2_struct.CCN2,1);
            
            ittccen =1 + (Ncs - 1)/2 ;% index of col center
            ittrcen = 1 + (Nrs - 1)/2; % index of row center        

            CCN2V = Single_scan_struct.CCN2_struct.CCN2;
            
            Ncq = 2*wcq + 1;
            Nrq = 2*wrq + 1;
            for icq = 1:Ncq
                offttc = (icq-wcq-1)*(2*hwttc+1)+ offsetcc;               
                ittc = round(ittccen + offttc + [-hwttc:hwttc]) ;
               
                for irq = 1:Nrq
                    offttr = (irq - wrq - 1)*(2*hwttr+1)+ offsetrc;
                    ittr = round(ittrcen + offttr + [-hwttr:hwttr]);
                    
                    Qval_struct = XPCS_analysis.calculate_qval(ittccen,ittrcen,ittccen + offttc,ittrcen + offttr,D_ds,kvector,pixel_size,th_Bragg);

                   
                    CCN2avg_struct.scancq(icq).scanrq(irq).CCN2avg = squeeze(mean(mean(CCN2V(ittr,ittc,:,:),1),2));
                    CCN2avg_struct.scancq(icq).scanrq(irq).timex = Single_scan_struct.IIbin_struct.Xamountb;
                    CCN2avg_struct.scancq(icq).scanrq(irq).TITLEstr2V = Single_scan_struct.IIbin_struct.TITLEstuct.TITLEstr2;
                    CCN2avg_struct.scancq(icq).scanrq(irq).nu = Qval_struct.nu;
                    CCN2avg_struct.scancq(icq).scanrq(irq).del = Qval_struct.del;
                    CCN2avg_struct.qvector.nu(irq) = Qval_struct.nu;
                    CCN2avg_struct.boxcenterrc.offttr(irq) = ittrcen + offttr;
                end
               CCN2avg_struct.qvector.del(icq) = Qval_struct.del;
               CCN2avg_struct.boxcenterrc.offttc(icq) = ittccen + offttc;
            end
            CCN2avg_struct.Nts_Ncs_Nrs = [Nts Ncs Nrs];
            CCN2avg_struct.Ncq_Nrq = [Ncq Nrq];
            CCN2avg_struct.TITLEstruct = Single_scan_struct.CCN2_struct.TITLEstuct;
            CCN2avg_struct.ittccen = ittccen ;
            CCN2avg_struct.ittrcen = ittrcen ;
            CCN2avg_struct.hwttr = hwttr;
            CCN2avg_struct.hwttc = hwttc;
            CCN2avg_struct.wrq = wrq;
            CCN2avg_struct.wcq = wcq;
            CCN2avg_struct.offsetcc = offsetcc;
            CCN2avg_struct.offsetrc = offsetrc;
        end
                
        function [CCN2S_struct] = from_CCN2avg_to_CCN2S(CCN2avg_struct)
            % This function builds the CCN_structure where we store the
            % integrated 2-times correlation function. The inputs are:
            %   Nrq: the number of pixels/number of q's
            %   CCN_struct: the structure containing the integrated 2-times correlation function
            %   TITLEstr2V: the title string which specifies the name of
            %   the scans in iiT
            
            
            
            Ndt = CCN2avg_struct.Nts_Ncs_Nrs(1);
            idt = 1:Ndt;
            
            
            Ncq = CCN2avg_struct.Ncq_Nrq(1);
            Nrq = CCN2avg_struct.Ncq_Nrq(2);
            
            for icq = 1:Ncq
                
                CCNdt = NaN*ones(Nrq,Ndt);
                
                for irq = 1:Nrq
                    time_scan = CCN2avg_struct.scancq(icq).scanrq(irq).timex - CCN2avg_struct.scancq(icq).scanrq(irq).timex(1); % Delta in time from region start
                    CCN2S = CCN2avg_struct.scancq(icq).scanrq(irq).CCN2avg(idt,idt);
                    for ii = idt
                        CCNdt(irq,ii) = mean(diag(CCN2S,ii-1));
                    end
                    CCN2S_struct.scancq(icq).scanrq(irq).CCNdtV = CCNdt(irq,:);
                    CCN2S_struct.scancq(icq).scanrq(irq).time_1D = time_scan;
                    CCN2S_struct.scancq(icq).scanrq(irq).nu = CCN2avg_struct.scancq(icq).scanrq(irq).nu;
                    CCN2S_struct.scancq(icq).scanrq(irq).del = CCN2avg_struct.scancq(icq).scanrq(irq).del;
                end
            end
            CCN2S_struct.Ncq_Nrq = CCN2avg_struct.Ncq_Nrq;
            CCN2S_struct.Ndt = Ndt;
            CCN2S_struct.ittccen= CCN2avg_struct.ittccen;
            CCN2S_struct.ittrcen =  CCN2avg_struct.ittccen ;
            
            CCN2S_struct.TITLEstruct = CCN2avg_struct.TITLEstruct;
        end
        
        function [CCN2Sfit_struct] =fit_CCN2S(CCN2S_struct,fitfunc_param_str,fitfunc_str,param_legend,fit_range,fig_plot_fit_results)
            
            % initialize fitting results
            a1 = zeros(numel(CCN2S_struct),numel(CCN2S_struct(1).scanq),1);
            b1 = zeros(numel(CCN2S_struct),numel(CCN2S_struct(1).scanq),1);
            
          
            figh = figure(fig_plot_fit_results);
            clf;
            
            for iT = 1:numel(CCN2S_struct)
                for irq = 1:numel(CCN2S_struct(iT).scanq)
                    p_lower =[0 0 0]; 
                    p_upper =[7E-1 100 1e-2]; 
                    p_start = [6e-2 50 1e-3];
                    
                    [CCN2Sfit_struct(iT).scanq(irq).CCNdtV_fit] = FittingFunctions.fit_2time_corr(CCN2S_struct(iT).scanq(irq),fitfunc_param_str,param_legend,fitfunc_str,fit_range,p_lower,p_upper,p_start);                              

                    %DisplayFunctions_XPCS.display_fit_result(CCN_struct,iT,irq,figh);
                    
                end
               
             CCN2Sfit_struct(iT).TITLEstruct = CCN2S_struct.TITLEstruct;                        
            end
        end
        
        function [CCN2S_struct] =fit_CCN2S_with_leasqr(CCN2S_struct,iT,pin_iiT,dp_iiT,fitfunc_str,fit_range,indexq,flag_row_or_col,fig_plot_fit_results)
            
            figure(fig_plot_fit_results);
            clf;
            
            switch flag_row_or_col
                
                case 'row'
                    qtolook = [ 1:CCN2S_struct.Ncq_Nrq(2)];
                case 'col'
                    qtolook = [ 1:CCN2S_struct.Ncq_Nrq(1)];
                otherwise
                    disp('please select rows or cols')
                    return;
            end
            
            for iq = qtolook
                
                
                switch flag_row_or_col
                    
                    case 'row'
                        CCN2S_struct_singlescan = CCN2S_struct.scancq(indexq).scanrq(iq);
                        
                    case 'col'
                        CCN2S_struct_singlescan = CCN2S_struct.scancq(iq).scanrq(indexq);
                        
                    otherwise
                        disp('please select rows or cols')
                        return;
                end
                
                pin = pin_iiT(iT,:);
                dp =  dp_iiT(iT,:);
                w = ones(length(fit_range),1);
                
                [ fitres] = FittingFunctions.fit_2time_corr_with_leasqr(CCN2S_struct_singlescan,fitfunc_str,fit_range,pin,dp, w);
                
                %                     CCN2Sfit_struct(iT).TITLEstruct = CCN2S_struct(iT).TITLEstruct;
                %
                switch flag_row_or_col
                    case 'row'
                        CCN2S_struct.scancq(indexq).scanrq(iq).CCNdtV_fit =  fitres;
                    case 'col'
                        CCN2S_struct.scancq(iq).scanrq(indexq).CCNdtV_fit =  fitres;
                    otherwise
                        disp('please select rows or cols')
                        return;
                end
                
                
                
            end
            
            %DisplayFunctions_XPCS.display_fit_result(CCN_struct,iT,figh);
            
            
        end
        
        function [pout_struct] = fit_tau_with_leasqr(pout_struct,fitfunc_str,qrange)
                      
            
            for iT = 1:numel(pout_struct)  
                
                    pout_struct_fitrange.qvector = pout_struct(iT).qvector(pout_struct(iT).fitrange);
                    sigma = pout_struct(iT).sigma(pout_struct(iT).fitrange);
                    pout_struct_fitrange.tau = pout_struct(iT).tau(pout_struct(iT).fitrange);
                    [~,low_index] = find(pout_struct_fitrange.qvector>qrange(1));
                    fit_range = low_index;
                
                    pin = [1];
                    dp =  [1]'*0.0001;
                    w = 1./sigma;
                    
                    [ pout_struct(iT).pout] = FittingFunctions.fit_tau_with_leasqr(pout_struct_fitrange,fitfunc_str,fit_range,pin,dp, w);
                
                   
            end
        end
        
    end
    
end