#
# I have altered this routine to be compatible with zero-offset
# arrays
#
def k_locate(xx, n, x, j):

	#int ascnd;

	jl = -1
	ju = n
	ascnd=(xx[n-1] > xx[0])             # !!! boolean in python; will it work with the comparison 4 lines below?
    # start while
    while (ju-jl) > 1 :
		jm=(ju+jl) >> 1                 # !!! python has the same operator; does work exactly the same?
        if (x > xx[jm]) == ascnd :
			jl=jm
        else :
			ju=jm
    # end while

	j=jl

