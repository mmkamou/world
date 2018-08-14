;MMKAMOU SCRIPT
;plot wind from ECMWF's output data 
;at different level 
;using streamline
;wind barb
;ploting wind intensity
;overlay this all on a same map


begin
  f = addfile("ecmwf_frcst_March18_00h.nc","r")
  u    = f->u(14,5,:,:)   ; read in example data [2D only here]
  v    = f->v(14,5,:,:)
t  = f->time(14)
print(t)
; =================================================;
; create plots
; =================================================;
wtype = "png"
  u = u*1.94386                    ; Convert wind into knots
  v = v*1.94386
;copy_VarMeta(u,w)
;---Change the metadata
  u@units       = "kts"
  v@units       = "kts"
 w = wind_speed(u,v)
 printVarSummary(w)
name_file  = "ecmwf_wind_00hh_1000hpa"
  wks = gsn_open_wks("png",name_file)
;Set common resources for all plots
gsn_define_colormap(wks,"wh-bl-gr-ye-re")
  vec_res                = True
  vec_res@gsnFrame       = False
  vec_res@gsnDraw        = False
  vec_res@gsnLeftString  = ""
  vec_res@gsnRightString = ""
  vec_res@gsnMaximize    = True	
  stres                  = vec_res
  tf_res                 = vec_res
  mpres                  = vec_res
  mpres@mpOutlineOn      = True

  ;---Make sure vectors are drawn in "predraw" phase.
 ; vec_res@vcVectorDrawOrder        = "Predraw" 
;---Temperature filled contour plot
  tf_res                      = vec_res
  tf_res@cnFillOn             = True  
  tf_res@cnLevelSelectionMode = "ExplicitLevels" ;AutomaticLevels"
  tf_res@cnLevels             = ispan(10,30,3)
  tf_res@lbLabelFontHeightF   = 0.015
  tf_res@lbOrientation        = "Horizontal"; "Vertical"
  tf_res@pmLabelBarOrthogonalPosF = -0.005
  tf_res@lbTitleString            = "Wind Intensity (KTS)"
  tf_res@lbTitleOffsetF           = -0.25
  tf_res@lbTitleFontHeightF       = 0.02
  tf_res@lbLabelFontHeightF       = 0.015
  tf_res@lbLabelAutoStride        = False
;---Make sure contours are drawn in "predraw" phase.
 ; tf_res@cnFillDrawOrder      = "Predraw"
;---Control appearance of map.
  mpres                        = vec_res
  mpres@mpProjection           = "LambertEqualArea"
  mpres@mpLabelsOn             = True
  mpres@mpPerimOn              = True
  mpres@mpGridAndLimbOn        = False
  mpres@mpFillOn               = True
  mpres@mpOutlineOn            = True
 ; mpres@mpOutlineDrawOrder     = "PostDraw"
 ; mpres@mpFillDrawOrder        = "Predraw"
  mpres@mpOceanFillColor       = "lightskyblue1"
  mpres@mpLandFillColor        = "gray"
; Create, but don't draw, a vector, streamline, contour
; and map plot. Stride the vector and contour data to
; thin it a little.
;
  mpres@mpLimitMode            = "LatLon" 
  mpres@mpMinLatF              =   3.
  mpres@mpMaxLatF              =   28.
  mpres@mpMinLonF              = -26.
  mpres@mpMaxLonF              =  14.
  mpres@mpGeophysicalLineThicknessF=2 
 stres@stLineThicknessF      = 2.0         ; default is 1.0
 title = "1000hpa_Wind_over_west Africa"
 mpres@tiMainString           = title 
 mpres@gsnLeftString  = "00H"
 mpres@gsnRightString = "15 March 2018"
 stres@stLineColor = "NavyBlue"
  vcid = gsn_csm_vector(wks,u,v,vec_res)
  stid = gsn_csm_streamline(wks,u,v,stres)
  cnid = gsn_csm_contour(wks,w,tf_res)
  mpid = gsn_csm_map(wks,mpres)
;---Overlay contour, streamline, and vector plots on the map plot.
  overlay(mpid,cnid)
  overlay(mpid,stid)
  overlay(mpid,vcid)
  lnres                  = True
  lnres@gsLineColor      = "black"
  lnres@gsLineThicknessF = 3.0
  lines_id = gsn_add_shapefile_polylines(wks,mpid,"../shape/countries.shp",lnres) 
 ;draw(mpid)
; frame(wks)
  maximize_output(wks,True)
end
