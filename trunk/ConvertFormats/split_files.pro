; to split incoming galaxy truth files


pro split_function, inputfile_base, num_new_files=num_new_files

    if ~keyword_set(num_new_files) then num_new_files = 4

    ; read in galaxies
    gal  = mrdfits(inputfile, 1)

    ; count the number of galaxies
    ngal = n_elements(gal)

    ; divide by the number of files we want
    ngal_new = ceil(float(ngal)/float(num_new_files))

    ;;; START LOOP
    start = 0
    finish = 0
    for i=0,num_new_files-1 do begin

        ; set last index of new sub-set of galaxies
        ;   -- use the basic number given from dividing by the number of files wanted added to the starting index
        ;   or
        ;   -- use the largest number in the original catalog
        finish = (start + ngal_new) < (ngal-1)

        ; get new galaxies using these subsets
        gal_new = gal[start:finish]

        ; set new starting point for galaxy indices
        start = finish + 1


    endfor
    ;;; END LOOP

    return

end


    

