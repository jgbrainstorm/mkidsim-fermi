NR_END  = 1
#define FREE_ARG char*


# allocate a float vector with subscript range v[nl..nh] */      # not sure there's an equivalent in python
#float *k_vector(long nl, long nh)
def k_vector(nl, nh):
	float *v;

	v=(float *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(float)));
	if (!v) {
		fprintf(stderr,"allocation failure in vector()");
		exit(1);
	}
	return v-nl+NR_END;
}

void k_free_vector(float *v, long nl, long nh)
/* free a float vector allocated with vector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}
