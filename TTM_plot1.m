% This scripts display the set of intensities taken during a time scan,
% their integrated intensity in the different ROIS

Qvalflag = 0;
[Allscans(iT).HS_struct] = DisplayFunctions_XPCS.display_IInormb(Read_Allscans(iT).IIstruct,Read_Allscans(iT).IIstruct.IInormb,'Linear',Qvalflag,13,ImageJ,SINGLE,XCOLlabel,YROWlabel,AXISdet,INFOstr);


[Allscans(iT).HS_struct] = DisplayFunctions_XPCS.display_IInormb_with_rois(Read_Allscans(iT),Allscans(iT).ROIS_struct,Read_Allscans(iT).IIstruct.IInormb(:,:,Allscans(iT).itt_range_struct.ittt),'Log 10',iT+13,ImageJ,SINGLE,XCOLlabel,YROWlabel,AXISdet,INFOstr,CLIM);

[Allscans(iT).HS_struct] = DisplayFunctions_XPCS.show_rois_only(Read_Allscans(iT),Allscans(iT),iT+14,XCOLlabel,YROWlabel,AXISdet);

SoXFLAG = 1;SoYFLAG = 0;LOGFLAG = 1;
[Allscans(iT).HS_sumimag_struct] = DisplayFunctions_XPCS.make_summed_images(Read_Allscans(iT),Allscans(iT).ROIS_struct,SoXFLAG,SoYFLAG,LOGFLAG,iT+15,ImageJ,XCOLlabel,YROWlabel,AXISdet);
[Allscans(iT).HS_sumimag_struct] = DisplayFunctions_XPCS.display_sumpoints_1D(Read_Allscans(iT),Allscans(iT),SoXFLAG,SoYFLAG,iT+16,ImageJ,XCOLlabel,YROWlabel,AXISdet);

SoXFLAG = 0;SoYFLAG = 1;LOGFLAG = 1;
[Allscans(iT).HS_sumimag_struct] = DisplayFunctions_XPCS.make_summed_images(Read_Allscans(iT),Allscans(iT).ROIS_struct,SoXFLAG,SoYFLAG,LOGFLAG,iT+17,ImageJ,XCOLlabel,YROWlabel,AXISdet);
[Allscans(iT).HS_sumimag_struct] = DisplayFunctions_XPCS.display_sumpoints_1D(Read_Allscans(iT),Allscans(iT),SoXFLAG,SoYFLAG,iT+18,ImageJ,XCOLlabel,YROWlabel,AXISdet);

DisplayFunctions_XPCS.plot_summed_images(Read_Allscans(iT),Allscans(iT).ROIS_struct,(Read_Allscans(iT).IIstruct.IInormb),iT+19,ImageJ,CLIM,XCOLlabel,YROWlabel,AXISdet,INFOstr,POINTSUMS);

ROIS_to_plot = 1;%[2,3,4];
[Allscans(iT).ROIS_struct.sum] = DisplayFunctions_XPCS.plot_ROIs_as_counters(Read_Allscans(iT),Allscans(iT),ROIS_to_plot,iT+20) ;