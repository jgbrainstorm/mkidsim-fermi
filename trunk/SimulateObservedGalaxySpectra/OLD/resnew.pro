;Zoubian lines included
;line width taken into consideration
;OII flux split-up
;Jeff's line widths being used instead of Zoubian's
;Jeffs' noise adding approach
;OII split-up only if line width allows for it

exposure=3600.*120.  ;exposure time in secs  
detection1=5.    ;S/N required for line detection
detection2=3.
detection3=2.47    ;S/N required for line detection
ndetect1=2L
ndetect2=3L 
ndetect3=4L      ;# of detections required in detection2 criterion
ndetect4=1L
maglim=25.3
specmin=3800.
specmax=13000.
typecut=8L   ;no ellipticals
noisechoice=2L  ;0: 0 deg. zenith - 1:60 deg. zenith - 2: average
withly=0L    ;0: no Ly-alpha. 1: with Ly-alpha
corr=3. ;Line flux correction factor
oiisplit=10000. ;Width of vel. dispersion in Ang. after which OII doublet can't be resolved
;======================================================================
;input/output file-naming
;======================================================================

if withly EQ 1 then begin
;h_cat = 'ssr_exp_'+strtrim(exposure,1)+'det1_'+strtrim(detection1,1)+'nd1_'+strtrim(ndetect1,1)+'det2_'+strtrim(detection2,1)+'nd2_'+strtrim(ndetect2,1)+'nc_'+strtrim(noisechoice,1)+'.dat'
h_cat = 'ssr_exp_'+strtrim(exposure,1)+'.dat'
endif else if withly EQ 0 then begin
;h_cat = 'ssr_exp_'+strtrim(exposure,1)+'det1_'+strtrim(detection1,1)+'nd1_'+strtrim(ndetect1,1)+'det2_'+strtrim(detection2,1)+'nd2_'+strtrim(ndetect2,1)+'nc_'+strtrim(noisechoice,1)+'.noly.dat'
h_cat = 'ssr_exp_'+strtrim(exposure,1)+'.noly.dat'
endif
sn_cat = 'sn_exp_'+strtrim(exposure,1)+'nc_'+strtrim(noisechoice,1)+'.dat'



if (noisechoice EQ 0L) then begin
   f_cat = 'pfssndata/noise_best.dat'
   f_cat2 = 'pfsdata/total_best.dat'
endif else if (noisechoice EQ 1L) then begin
   f_cat = 'pfssndata/noise_worst.dat'
   f_cat2 = 'pfsdata/total_worst.dat'
endif else if (noisechoice EQ 2L) then begin
   f_cat = 'pfssndata/noise_avg.dat'
   f_cat2 = 'pfsdata/total_avg.dat'
endif

;=======================================================================================
; Survey specs
;=======================================================================================


hp=6.62607e-27                  ;ergs.secs
c=2.9979e18                     ;Ang/secs
pi=3.1415926535897
Area=pi*(8.2/2)^2         ;telescope area

Areacm=Area*10000.              ;converting from m^2 to cm^2
Nres=1.7      ;# of pixels per resolution element (aka dispersion)


;============================================
;READ NOISE FILE
;============================================
    d_cat = './'
    f_in_cat = d_cat+'/'+f_cat
    format = "D,D"
    readcol,f_in_cat $
      , lam $
      , sn $
      , no $
      , noiseatm $
      , noisetherm $
      , noiseread $
      , format=format $
      ,/silent
    
    
;============================================
;
;============================================
  
NL=N_elements(lam)

;============================================
;READ TOTAL TRANSMISSION  FILE
;============================================
    d_cat = './'
    f_in_cat2 = d_cat+'/'+f_cat2
    format = "D,D"
    readcol,f_in_cat2 $
      , abswave $
      , absatm $
      , format=format $
      ,/silent
    
    Nabs=n_elements(abswave)


  
;============================================
;interpolate absatm to same grid as atm
;============================================
Atmtrans=DBLARR(NL)
Atmtrans[*]=interpol(absatm,abswave,lam)
 


;8888888888888888888888888888888888888888888888888888888888888888888888888
;8888888888888888888888888888888888888888888888888888888888888888888888888
;=========================================================================
; Single-line flux/Angstrom
;=========================================================================
;8888888888888888888888888888888888888888888888888888888888888888888888888
;8888888888888888888888888888888888888888888888888888888888888888888888888


;=======================================================================
;Read - Cosmos sims line fluxes
;=======================================================================


in=mrdfits('../../newlinesims/CMC081211_all.fits',1,hdr) 

ct=0L
cte=0L
Na=N_elements(in.z)
for i=0L,Na-1L do begin
   if (in[i].type EQ 1L) then begin
      ct++
      if (in[i]._Mod GT typecut AND in[i].RAN_I_SUBARU LT maglim) then begin
         cte++
      endif
   endif 
endfor

Ns=cte
Nse=cte

;;
 nemis=11L+1L    
 Lflux=DBLARR(nemis,Ns)
 Ldisp=DBLARR(nemis,Ns)
 Llam=DBLARR(nemis,Ns)
 Roii=DBLARR(Ns)
 red=DBLARR(Ns)
 stype=DBLARR(Ns)
 Imag=DBLARR(Ns)

ct=0L
for i=0L,Na-1L do begin
   if (in[i].type EQ 1L AND (in[i]._Mod GT typecut AND in[i].RAN_I_SUBARU LT maglim)) then begin
      
      red[ct]=in[i].z
      Imag[ct]=in[i].RAN_I_SUBARU
;Roii[ct]=in[i].ROII
      Ratoii=1.3
      
      if withly EQ 1 then begin
         Llam[0,ct]=in[i].LAMBDA_LY
      endif else if withly EQ 0 then begin
         Llam[0,ct]=0
      endif
      
;=====================================
;OII needs special attention
;======================================
;; Ldisp[1,ct]=in[i].DISP_OII
;; Ldisp[11,ct]=in[i].DISP_OII
 

      ifoii=0L
      if in[i].DISP_OII LT oiisplit*(1.+in[i].z) then begin
         ifoii=1L
      endif

      if ifoii EQ 1L then begin   ;decide whether the doublet will be resolved or not.
         Llam[1,ct]=in[i].LAMBDA_OII -0.98*(1.+in[i].z)
         Llam[11,ct]=in[i].LAMBDA_OII +1.83*(1.+in[i].z) 
         Lflux[1,ct]=in[i].FLUX_OII/(1.+Ratoii)
         Lflux[11,ct]=in[i].FLUX_OII/(1.+1./Ratoii)
      endif else begin
         Llam[1,ct]=in[i].LAMBDA_OII
         Llam[11,ct]=0.
         Lflux[1,ct]=in[i].FLUX_OII
         Lflux[11,ct]=0.
      endelse

;=====================================
;The other lines
;======================================
      
      Llam[2,ct]=in[i].LAMBDA_HB
      Llam[3,ct]=in[i].LAMBDA_OIIIA
      Llam[4,ct]=in[i].LAMBDA_OIIIB
      Llam[5,ct]=in[i].LAMBDA_HA
      Llam[6,ct]=in[i].LAMBDA_HD
      Llam[7,ct]=in[i].LAMBDA_HG
      Llam[8,ct]=in[i].LAMBDA_NII
      Llam[9,ct]=in[i].LAMBDA_SIIA
      Llam[10,ct]=in[i].LAMBDA_SIIB

      ;; Lflux[0,ct]=in[i].FLUX_LY
      ;; Lflux[2,ct]=in[i].FLUX_HB 
      ;; Lflux[3,ct]=in[i].FLUX_OIIIA
      ;; Lflux[4,ct]=in[i].FLUX_OIIIB
      ;; Lflux[5,ct]=in[i].FLUX_HA
      ;; Lflux[6,ct]=in[i].FLUX_HD
      ;; Lflux[7,ct]=in[i].FLUX_HG
      ;; Lflux[8,ct]=in[i].FLUX_NII
      ;; Lflux[9,ct]=in[i].FLUX_SIIA
      ;; Lflux[10,ct]=in[i].FLUX_SIIB

      Lflux[0,ct]=in[i].FLUX_LY/corr
      Lflux[2,ct]=in[i].FLUX_HB /corr
      Lflux[3,ct]=in[i].FLUX_OIIIA/corr
      Lflux[4,ct]=in[i].FLUX_OIIIB/corr
      Lflux[5,ct]=in[i].FLUX_HA/corr
      Lflux[6,ct]=in[i].FLUX_HD/corr
      Lflux[7,ct]=in[i].FLUX_HG/corr
      Lflux[8,ct]=(in[i].FLUX_NII < in[i].FLUX_HA*10^(-0.5))/corr
      Lflux[9,ct]=in[i].FLUX_SIIA/corr
      Lflux[10,ct]=in[i].FLUX_SIIB/corr

      
;; Ldisp[0,ct]=in[i].DISP_LY

;; Ldisp[2,ct]=in[i].DISP_HB 
;; Ldisp[3,ct]=in[i].DISP_OIIIA
;; Ldisp[4,ct]=in[i].DISP_OIIIB
;; Ldisp[5,ct]=in[i].DISP_HA
;; Ldisp[6,ct]=in[i].DISP_HD
;; Ldisp[7,ct]=in[i].DISP_HG
;; Ldisp[8,ct]=in[i].DISP_NII
;; Ldisp[9,ct]=in[i].DISP_SIIA
;; Ldisp[10,ct]=in[i].DISP_SIIB
      

      
      Ldisp[*,ct]=Llam[*,ct]*(65./3.0d5)       ;Line widths in mocks seem suspicious

      if ifoii EQ 1L then begin
         Ldisp[1,ct]=in[i].DISP_OII
      endif

      ct++
   endif
endfor

;=======================================================================
;Convert fluxes to counts 
;=======================================================================



Ltrans=dblarr(nemis,Ns)
Lcounts=DBLARR(nemis,Ns)

for i=0,nemis-1L do begin
temp=Llam[i,*]
Ltrans[i,*]=interpol(Atmtrans,lam,temp[*])
Lcounts[i,*]=Lflux[i,*]*Llam[i,*]*Areacm*exposure*Ltrans[i,*]/(hp*c)
endfor



;=======================================================================
;Find pixel containing line center
;=======================================================================
pix=DBLARR(NL)
pixt=LONARR(nemis,Ns)
pwidth=DBLARR(nemis,Ns)
for i=0L,NL-1L do begin
   pix[i]=double(i)
endfor
for i=0,nemis-1L do begin
   temp=Llam[i,*]
   pixt[i,*]=long(interpol(pix,lam,temp[*]))
  
endfor
for i=0,nemis-1L do begin
 pwidth[i,*]=lam[pixt[i,*]]-lam[pixt[i,*]-1L]
endfor

for i=0L,nemis-1L do begin
   for j=0L,Ns-1L do begin
      if pwidth[i,j] GT 1 then pwidth[i,j]=0.8  ;fix gap in pixel width.
   endfor
endfor




;8888888888888888888888888888888888888888888888888888888888888888888888888
;=========================================================================
;Calculate S/N
;=========================================================================
;8888888888888888888888888888888888888888888888888888888888888888888888888



noisetot=dblarr(NL)
noisepart=dblarr(NL)
darknoise=dblarr(NL)
rdnoise=dblarr(NL)
Vatm=dblarr(NL)

nhour=3 ;# of readouts per hour
darknoise[*]=noisetherm[*]*noisetherm[*]*exposure/1000.
rdnoise[*]=noiseread[*]*noiseread[*]*nhour*exposure/3600.
Vatm[*]=noiseatm[*]*noiseatm[*]*exposure/1000./0.55  ;0.55 is Jeff's suggested correction


noisepart[*]= darknoise[*] +rdnoise[*] + Vatm[*]

npart=dblarr(nemis,Ns)
ntot=dblarr(nemis,Ns)
stot=dblarr(nemis,Ns)
ston=dblarr(nemis,Ns)


for i=0L,nemis-1L do begin
temp=Llam[i,*]
npart[i,*]=interpol(noisepart,lam,temp[*])
endfor

;=========================================================================
;calculate effective emission-line width
;its a combination of intrinsic width + resolution
;=========================================================================

sig_tot=DBLARR(Nemis,Ns) ;total line width = pix_width (*) intrinsic width (int_width)

offs=10L
Noff=2*offs+1L



for i=0L,Nemis-1L do begin
   sig_tot[i,*]=sqrt((Nres*pwidth[i,*])^2 + Ldisp[i,*]^2 )  ;effective line width
   signal=dblarr(Noff,Ns)
   frac_in_pix=dblarr(Noff,Ns)
   noise_sq=dblarr(Noff,Ns)
   for j=0L,Noff-1L do begin
      sigtemp=dblarr(Ns)
      noisetemp=dblarr(Ns)
      sntemp=dblarr(Ns)
      place=j-Noff+offs+1L

      xf=(lam[pixt[i,*]+place]-Llam[i,*])/(sig_tot[i,*]*sqrt(2.))
      xi=(lam[pixt[i,*]+place+1L]-Llam[i,*])/(sig_tot[i,*]*sqrt(2.))
      fraction=0.5*(erfc(xf)-erfc(xi))

      sigtemp = Lcounts[i,*]*fraction
      signal[j,*]=sigtemp

      noisetemp[*] = noisepart[pixt[i,*]+place]
      frac_in_pix[j,*] = fraction    ; fraction of all signal in the j'th pixel for each spectrum
 
      noise_sq [j,*] = sigtemp + noisetemp ; variance in the j'th pixel for each spectrum

     
   endfor
   

   ;for k = 0L, ns-1 do begin 
    ;  frac_in_pix[*,k] = frac_in_pix[*,k] / total(frac_in_pix[*,k]) ; this normalizes the fractions to sum to 1
 ;  endfor

   sigtot = total(signal,1)    ; this will total the array over the first axis, i.e. summing over the 5 pixels for each spectrum

   noisetot = 1/sqrt( total(  frac_in_pix^2/noise_sq , 1) )
   ston[i,*] = sigtot/noisetot

;   if i EQ 4 then stop

endfor


cdetect1=LONARR(Ns)
cdetect2=LONARR(Ns)
cdetect3=LONARR(Ns)
cdetect1[*]=0L
cdetect2[*]=0L
cdetect3[*]=0L

ctot=0L
ctot1=0L
ctot2=0L
ctot3=0L
cmlim=0L

for i=0L, Ns-1L do begin
   if Imag[i] LT maglim then begin
      cmlim++
      for j=0L, nemis-1L do begin
         if (ston[j,i] GE detection3 AND Llam[j,i] GT specmin AND Llam[j,i] LT specmax) then begin
            cdetect3[i]++
            if (ston[j,i] GE detection2) then begin
               cdetect2[i]++
               if (ston[j,i] GE detection1) then begin
                  cdetect1[i]++
               endif
            endif
         endif
      endfor
   endif
   if ((cdetect1[i] GE ndetect1) OR (cdetect2[i] GE ndetect2) OR (cdetect2[i] GE ndetect3) OR ((cdetect1[i] GE ndetect4) AND (cdetect3[i] GE ndetect2))) then begin
      ctot++
   endif
   if cdetect2[i] GE ndetect2 then begin
      ctot2++
   endif
   if cdetect1[i] GE ndetect1 then begin
      ctot1++
   endif
   if cdetect3[i] GE ndetect3 then begin
      ctot3++
   endif
endfor  


print,'Fraction detected', double(ctot)/double(cmlim), double(ctot1)/double(cmlim), double(ctot2)/double(cmlim), double(ctot3)/double(cmlim)

;====================================================
;histogram stuff
;====================================================
zmin=0.
zmax=6.
imin=20.
imax=maglim
nzbins=40
nibins=40


hz=dblarr(nzbins)
hi=dblarr(nibins)
hzdet=dblarr(nzbins)
hidet=dblarr(nibins)
hztot=dblarr(nzbins)
hitot=dblarr(nibins)
zdetfrac=dblarr(nzbins)
idetfrac=dblarr(nibins)

hzdet[*]=0
hidet[*]=0
hztot[*]=0
hitot[*]=0

delz=(zmax-zmin)/double(nzbins)
deli=(imax-imin)/double(nibins)

for i=0L,nzbins-1 do begin
hz[i]=zmin+i*delz
endfor
for i=0L,nibins-1 do begin
hi[i]=imin+i*deli
endfor

countz=0L

for i=0L,Ns-1 do begin
   for j=0L,nzbins-1 do begin
      if ((red[i] GE hz[j]) AND (red[i] LT (hz[j]+delz))  AND Imag[i] LT maglim ) then begin
         hztot[j]++
         countz++
         if ((cdetect1[i] GE ndetect1) OR (cdetect2[i] GE ndetect2) OR (cdetect2[i] GE ndetect3) OR ((cdetect1[i] GE ndetect4) AND (cdetect3[i] GE ndetect2))) then begin
            hzdet[j]++            
         endif
      endif
   endfor
endfor

zdetfrac[*]=hzdet[*]/hztot[*]

;==================================================
counti=0L

for i=0L,Ns-1 do begin
   for j=0L,nibins-1 do begin
      if ((Imag[i] GE hi[j]) AND (Imag[i] LT (hi[j]+deli))) then begin
         hitot[j]++
         counti++
         if ((cdetect1[i] GE ndetect1) OR (cdetect2[i] GE ndetect2) OR (cdetect2[i] GE ndetect3) OR ((cdetect1[i] GE ndetect4) AND (cdetect3[i] GE ndetect2))) then begin
            hidet[j]++            
         endif
      endif
   endfor
endfor
idetfrac[*]=hidet[*]/hitot[*]
;=========================================================
;plot,hz,zdetfrac
;plot,hi,idetfrac


openw,1,sn_cat 
for i=0L,Ns-1 do begin
printf,1,format='(14F14.6, 2I)',red[i],Imag[i],ston[0,i],ston[1,i],ston[2,i],ston[3,i],ston[4,i],ston[5,i],ston[6,i],ston[7,i],ston[8,i],ston[9,i],ston[10,i],ston[11,i],cdetect1[i],cdetect2[i]
endfor
close,1

openw,2,h_cat  
for i=0L,nzbins-1 do begin
printf,2,format='(8F14.6 , 2I)',hz[i],hzdet[i],hztot[i],zdetfrac[i],hi[i],hidet[i],hitot[i],idetfrac[i],countz,counti
endfor
close,2




print,'FINISHED'
;end








print,'FINISHED'
end
