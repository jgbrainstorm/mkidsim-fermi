

pro cull_galaxy_table, inputfile_base, inputdir=inputdir

    if ~keyword_set(inputdir) then inputdir = '../../data/'

    ; set full file name of input file
    inputfile = inputdir + inputfile_base + '.fit'
    print, '***'
    print, 'inputfile: ', inputfile

    ; read input file
    gal = mrdfits(inputfile,1)

    ; count galaxies in input catalog
    ngal = n_elements(gal)

    ;set up new array of structures from input catalog
    gal_new = replicate({z:0.0D, id:0L, coeffs:fltarr(5)},ngal)

    ; copy values from original to new subset catalog
    gal_new.z       = gal.z
    gal_new.id      = gal.id
    gal_new.coeffs  = gal.coeffs

    ; set output file name
    outputfile = inputdir + inputfile_base + '_tag_subset.fit'
    print, 'outputfile: ', outputfile
    print, '***'
    print

    ; write output file
    mwrfits,gal_new,outputfile,/create




end


pro cull_galaxy_table_main


    inputfile_base_list = [ 'Aardvark_v0.5d_truth_des_masked.114',$
                            'Aardvark_v0.5d_truth_des_masked.115',$
                            'Aardvark_v0.5d_truth_des_masked.147',$
                            'Aardvark_v0.5d_truth_des_masked.86' ]

    for i=0,n_elements(inputfile_base_list)-1 do cull_galaxy_table, inputfile_base_list[i]
       


end
