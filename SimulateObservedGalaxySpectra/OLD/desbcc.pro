pro desbcc,expnum,resnum,in_0,in_f ; 

; expnum = exposure time in seconds

suf=string('.fit')
expint=long(expnum)
resint=long(resnum)
explab=string(expint)
explab=strtrim(explab,1)
reslab=string(resint)
reslab=strtrim(reslab,1)

i_0=long(in_0)
i_f=long(in_f)

indx_0=string(i_0)
indx_f=string(i_f)

indx_0=strtrim(indx_0,1)
indx_f=strtrim(indx_f,1)

fnames=strarr(330)
files=findfile('/nfs/slac/g/ki/ki11/des/mbusha/catalogs/BCC1/Catalog/Aard*truth_no_photo*.fit')
Nfile=n_elements(files)
print,Nfile


;read a few of the input files
for indx=i_0,i_f do begin

filename=files[indx]
print,filename
infile=mrdfits(filename,1, hdr) 

id=infile.id
Nrows=n_elements(id)

print,"N_elements = ",Nrows

;=======================================================================
; Calculate  coeffs using kcorrect.
;=======================================================================

Nr=1L
;analysis is only performed for a subsample of the galaxies in each
;file

for jj=0L, Nr-1 do begin
	
	Nskip=392L*4  ;1 out of Nskip objects are used
	;Nskip=50L
	Nshort=long(((double(Nrows)/double(Nskip))+1.)/double(Nr))
	print,'Nshort = ',Nshort,"==================================="
	Ni=jj*Nshort
	print,"FILENUMBER = ",indx+1L,"out of ",i_f-i_0,"jj = ",jj,"   Ni = ",Ni,"(Ni+Nshort-1)*Nskip =" ,(Ni+Nshort-1)*Nskip ;Ni=0L

	tmagtemp=FLTARR(5,Nshort)
	tmagerrtemp=FLTARR(5,Nshort)
	zspectemp=FLTARR(Nshort)
	
	listname="/nfs/slac/g/ki/ki11/des/carlos/bcclist."+indx_0+"."+indx_f+".dat"
	listname=strtrim(listname,1)

	openw,1,listname,/append  ;CAREFUL! Output always appended, so remember to delete old files	 
	count =0L
	for i=Ni,Ni+Nshort-1L do begin
		ii=i-Ni

		tmagtemp[*,ii]=infile[i*Nskip].tmag[*]
		tmagerrtemp[*,ii]=0.01;infile[i].tmagerr[i] ;needed for kcorrect call (set to low value since error free mags are used).
		zspectemp[ii]=infile[i*Nskip].z


;print ascii table with photometric properties of galaxies for which
;spectrum will be calculated
		printf,1, format='(I,I, 18F11.6,I)',infile[i*Nskip].id,i*Nskip,infile[i*Nskip].ra,infile[i*Nskip].dec,zspectemp[ii],tmagtemp[0,ii],tmagtemp[1,ii],tmagtemp[2,ii],tmagtemp[3,ii],tmagtemp[4,ii],infile[i*Nskip].omagerr[*],infile[i*Nskip].omag[*],indx
	endfor

	close,1


	print,'BEGIN KCORRECT'
	
	kcorrect,tmagtemp,tmagerrtemp,zspectemp,kcorrect,rmatrix=rmatrix,zvals=zvals,coeffs=coeffs,vname=vname,mass=mass,omega0=0.25,omegal0=0.75,chi2=chi2,mtol=mtol,filterlist=['DES_g.par','DES_r.par','DES_i.par','DES_z.par','DES_Y.par'],/magnitude 
	

	print,'DONE WITH KCORRECT'
	
;=================================================================================
;; Calculate SED's based on coeffs.
;=================================================================================
	
	
	k_reconstruct_spec,coeffs,loglam,mass=mass,flux,vdisp=300,vname=vname
	;outflux=k_smooth(loglam,flux,1219/2.35)   ;smooth galaxy spectra to survey resolution
	;flux=outflux
	
	NL=10000
								;Nshort=Nshort
	lam=DBLARR(NL)
	lam[*]=10^loglam[*]
	
;============================================
;READ SKY BACKGROUND FILE
;============================================
	d_cat = './'
	f_cat = 'skybg_50_10.dat'
	f_in_cat = d_cat+'/'+f_cat
	format = "D,D"
	readcol,f_in_cat $
	  , wave $
	  , temp $
	  , format=format $
	  ,/silent
	
	N=n_elements(wave)

	atm=DBLARR(N)
	
;=======================================================================================
;;Convert sky photons to photons/Ang.
;in skybg_50_10.dat the units are photons/s/nm/arcsec^2/m^2
;=======================================================================================
	pi=3.1415926535897
	length=5.				   ;arcsecs
	width=0.5				   ;arcsecs
	aperture=length*width

	exposure=expnum

;exposure=2400.00			;wide VVDS
;exposure=16200.			 ;deep VVDS
;exposure=40000.			 ;ultra-deep VVDS (I made up the exposure time)

	Area=pi*4^2
	Areacm=Area*10000.		  ;converting from m^2 to cm^2
	worsening=1.
	atm[*]=worsening*exposure*temp[*]*aperture*Area/10. ;because 10ang=1 nm
	wave[*]=10.*wave[*]
;=======================================================================================
;;READ TRANSMISSION FILE
;=======================================================================================
	
	d_cat = './'
	f_cat = 'VIMOS_total_transmission.dat'
	f_in_cat = d_cat+'/'+f_cat
	format = "D,D"
	readcol,f_in_cat $
	  , wavetrans $
	  , trans $
	  , format=format $
	  ,/silent
	
	Ntrans=n_elements(wavetrans)
	wavetrans[*]=10.*wavetrans[*]
	trans=0.01*trans
	
;============================================
;READ ATMOSPHERIC ABSORPTION  FILE
;============================================
	d_cat = './'
	f_cat = 'palomarextinct.dat'
	f_in_cat = d_cat+'/'+f_cat
	format = "D,D"
	readcol,f_in_cat $
	  , abswave $
	  , abstemp $
	  , format=format $
	  ,/silent
	
	Nabs=n_elements(abswave)
	absatm=DBLARR(Nabs)
	
	airmass=1.3
	altitude=2635
	absatm[*]=abstemp[*]*airmass*exp((1700.-altitude)/7000.)
	for i=0L,Nabs-1 do begin
		if (absatm[i] GT 1.) then begin
			absatm[i]=1.
		endif
		
	endfor

;====================================
;smooth spectra to VVDS resolution
	logwave=alog10(wave)
	;outatm=k_smooth(logwave,atm,1320)
	outatm=k_smooth(logwave,atm,300)
	atm=outatm

;I think I shouldn't smooth over atmospheric absorption since 
;smoothing only happens after the photons have been removed from atmosphere.
	;logabs=alog10(abswave) ;
	;outabs=k_smooth(logabs,absatm,1320)
	;absatm=outabs
;=======================================================================================
;;Convert sed flux to photon counts.
;sed flux is in ergs/cm^2/s/Ang.
;photon counts is in photons/Ang.
;=======================================================================================
	sedphot=DBLARR(NL,Nshort)
	sedflux=DBLARR(NL,Nshort)
	hp=6.62607e-27			  ;ergs.secs
	c=2.9979e18				 ;Ang/secs
	
	for i=0L,Nshort-1 do begin
		flux[*,i]=flux[*,i]/(1.+zspectemp[i]) ;add redshift dimming
		;flux[*,i]=flux[*,i]/(1.+zspectemp[i])^2 ;add redshift dimming
		sedphot[*,i]=exposure*flux[*,i]*lam[*]*Areacm/(hp*c) ;approximate - correct is to integrate 
		sedflux[*,i]=exposure*flux[*,i]*Areacm ; 
	endfor
;=======================================================================================
;;Redshift galaxy sed's. 
;;rest-frame flux assumes galaxy is 10pc away.
;doesn't matter when using coeffs calculated from
;kcorrect, since the code already applies the 1/r^2 dimming.
;=======================================================================================
	lamz=DBLARR(NL,Nshort)
	for i=0L,Nshort-1 do begin
		lamz[*,i]=lam[*]*(1.+zspectemp[i])
	endfor
;=======================================================================================
;;Interpolate wavelengths onto grid with resolution
;set by VVDS resolution R=230=lam/delta_lam.
;1) make grid 
;2) interpolate counts/fluxes onto grid
;3) Note that delta_lam will affect number of photons counted.
;that is, remember to multiply by delta_lam.
;=======================================================================================
	;Nwav=557
	;R=230.*3.		  ;The 3 converts from resolution to dispersion. 
	;dellam=7.14

	;Nwav=557*5.
	;R=230.*3.*5.		  ;The 3 converts from resolution to dispersion. 
	;dellam=7.14/5.

	Nwav=557*resnum
	R=230.*3.*resnum		  ;The 3 converts from resolution to dispersion. 
	dellam=7.14/resnum		;R is not being used. Only dellam matters now.

	print,"EXPOSURE = ",expnum,"NWAV = ",Nwav,"resnum = ",resnum

	lam0=5500.
	lamvvds=DBLARR(Nwav)
	binwidth=DBLARR(Nwav)
	Vsed=DBLARR(Nwav,Nshort)
	Vflux=DBLARR(Nwav,Nshort)
	Vtrans=DBLARR(Nwav)
	Vabsatm=DBLARR(Nwav)
	lamvvds[0]=lam0
	
	for i=1L,Nwav-1 do begin
		;lamvvds[i]=lamvvds[i-1]*(1.+1./R)
		lamvvds[i]=lamvvds[i-1]+dellam
	endfor
	for i=0L,Nwav-2 do begin
		lamvvds[i]=lamvvds[i]+(lamvvds[i+1]-lamvvds[i])/2.
		binwidth[i]=(lamvvds[i+1]-lamvvds[i])
	endfor
;=======================================================================================
;;Bin spectra to survey's resolutions by integrating over pixels
;to convert flux density to flux. 
;????????????????????
;ISSUE:
;Not sure if this should be done before or after redshifting templates
;??????????????????????
;=======================================================================================
	Vatm= k_binspec(wave,atm,lamvvds)
	
	for i=0L,Nshort-1 do begin
		Vsed[*,i]=k_binspec(lamz[*,i],sedphot[*,i],lamvvds)
		Vflux[*,i]=k_binspec(lamz[*,i],sedflux[*,i],lamvvds)
	endfor
	
;=======================================================================================
;;Multply fluxes by telescope transmission
;=======================================================================================
	Vtrans[*]=interpol(trans,wavetrans,lamvvds[*])
	Vabsatm[*]=interpol(absatm,abswave,lamvvds[*])
	Vatm[*]=Vatm[*]*Vtrans[*]*(1.-Vabsatm[*])
	for i=0L,Nshort-1 do begin
		Vsed[*,i]=Vsed[*,i]*Vtrans[*]*(1.-Vabsatm[*])
	endfor
;=======================================================================================
;;Calculate S/N
;=======================================================================================
	Vatmerr=DBLARR(Nwav)
	Vatmerr=sqrt(Vatm)
	signoise=DBLARR(Nwav,Nshort)
	readnoise=5
	for i=0L,Nshort-1 do begin
		Vatmerr=sqrt(Vatm+Vsed[*,i]+readnoise*readnoise)
		signoise[*,i]=Vsed[*,i]/(Vatmerr[*])
	endfor
	Vfluxerr=DBLARR(Nwav,Nshort)
	Vfluxerr=Vflux/signoise
;=======================================================================================
;;Generate observed spectra
;=======================================================================================
	
	Vfluxobs=DBLARR(Nwav,Nshort)
	gasdev=DBLARR(Nwav)
	
	
	
	for i=0L,Nshort-1 do begin
		
		for j=0L,Nwav-1 do begin
			gasdev[j]=(RANDOMN(seed,/normal,/double))
		endfor
		
		Vfluxobs[*,i]=Vflux[*,i]+gasdev[*]*Vfluxerr[*,i]
		if (expnum GT 1000000) then begin
			Vfluxobs[*,i]=Vflux[*,i]
		endif

	endfor
	tempname=string(indx)
	tempname=strtrim(tempname,1)
	ext="."+tempname+".fits"
	exta="."+tempname+".dat"

	for i=0,Nshort-1 do begin
		folder="carlos/ascii.bcc16200/"
		folder=strtrim(folder,1)
		name0=string(infile[(i+Ni)*Nskip].id)
		name0=strtrim(name0,1)
	   
		name=folder+name0
		name=strtrim(name,1)
		name2=name+exta
		name2=strtrim(name2,1) 
 

		openw,4,name2
		for j=0,Nwav-1 do begin
		printf,4,format='(3E)',lamvvds[j],Vfluxobs[j,i],Vflux[j,i]	
		endfor
		close,4

	endfor



	for i=0,Nshort-1 do begin
;name=string(id[i+Ni])
		folder="/nfs/slac/g/ki/ki11/des/carlos/bcc."
		folder=folder+explab+"."+reslab+"/"
		;folder="temp/"
		;folder="vvdsdeep.des.skip.hires/"
		;folder="vvdsdeep.des.skip.noerr/"
		folder=strtrim(folder,1)
		
								;name0=string(i+Ni)
								;name0=string(id[i+ni])
		;name0=string(infile[i+ni].id)
		name0=string(infile[(i+Ni)*Nskip].id)
		name0=strtrim(name0,1)
	   
		name=folder+name0
		name=strtrim(name,1)
		name2=name+ext
		name2=strtrim(name2,1) 
		aa=Vfluxobs[*,i]
		mwrfits,aa,name2,/create
		
		h = headfits(name2,exten=0)
		sxaddpar, h, 'CRVAL1', 5500
		sxaddpar, h, 'CRPIX1', 1
		sxaddpar, h, 'CTYPE1', 'LINEAR'
		sxaddpar, h, 'CD1_1', 7.14 ;1.428 ;7.14
								;sxaddpar, h, 'NAXIS', 557  
		
		modfits,name2,0,h,exten=0
								;print,"NAME = ",name2
	endfor
;#############################################################
;Create catalog with atmospheric lines removed
;#############################################################
;comment out everything below if you don't want to remove the
;lines at this stage.
	

		
	for i=0L,Nshort-1 do begin
		for j=0L,Nwav-1 do begin
			
			if ((j EQ 10) OR (j EQ 111) OR (j EQ 120) OR (j EQ 294) OR (j EQ 465)) then begin 
				dx1=1./3.
				Vfluxobs[j,i]=Vfluxobs[j+2,i]*dx1+Vfluxobs[j-1,i]*(1.-dx1)
			endif
			if ((j EQ 11) OR (j EQ 112) OR (j EQ 121) OR (j EQ 296) OR (j EQ 466)) then begin 
				dx1=2./3.
				Vfluxobs[j,i]=Vfluxobs[j+1,i]*dx1+Vfluxobs[j-2,i]*(1.-dx1)
			endif
			if (j EQ 295) then begin 
				dx1=1./2.
				Vfluxobs[j,i]=Vfluxobs[j+1,i]*dx1+Vfluxobs[j-2,i]*(1.-dx1)
			endif
			
			
		endfor
	endfor
	
	ext=".fits"
	exta=".dat"

; ;	 for i=0,Nshort-1 do begin
; ;		 folder="ascii.smt/"
; ;		 folder=strtrim(folder,1)
; ;		 name0=string(infile[(i+Ni)*Nskip].id)
; ;		 name0=strtrim(name0,1)
	   
; ;		 name=folder+name0
; ;		 name=strtrim(name,1)
; ;		 name2=name+exta
; ;		 name2=strtrim(name2,1) 
 

; ;		 openw,4,name2
; ;		 for j=0,Nwav-1 do begin
; ;		 printf,4,format='(3E)',lamvvds[j],Vfluxobs[j,i],Vflux[j,i]	
; ;		 endfor
; ;		 close,4

; ;	 endfor



	for i=0,Nshort-1 do begin

		folder="/nfs/slac/g/ki/ki11/des/carlos/bcc.atm."
		folder=folder+explab+"."+reslab+"/"
		folder=strtrim(folder,1)
		

		name0=string(infile[(i+Ni)*Nskip].id)
		name0=strtrim(name0,1)
	   
		name=folder+name0
		name=strtrim(name,1)
		name2=name+ext
		name2=strtrim(name2,1) 
		aa=Vfluxobs[*,i]
		mwrfits,aa,name2,/create
		
		h = headfits(name2,exten=0)
		sxaddpar, h, 'CRVAL1', 5500
		sxaddpar, h, 'CRPIX1', 1
		sxaddpar, h, 'CTYPE1', 'LINEAR'
		sxaddpar, h, 'CD1_1', 7.14 ;1.428 ;7.14
								;sxaddpar, h, 'NAXIS', 557  
		
		modfits,name2,0,h,exten=0
								;print,"NAME = ",name2
	endfor


endfor
;endif
endfor ;indx for loop

print,'FINISHED'
end
