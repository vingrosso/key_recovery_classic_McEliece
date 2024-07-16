/*
  This file is for syndrome computation
*/

#include "synd.h"

#include "params.h"
#include "root.h"
#include "simpleserial.h"

#include <stdio.h>



/* input: Goppa polynomial f, support L, received word r */
/* output: out, the syndrome of length 2t */
void synd(gf *out, gf *f, gf *L, unsigned char *r)
{
	int i, j;
	gf e, e_inv, c;
	
	for (j = 0; j < 2*SYS_T; j++)
		out[j] = 0;
	
	#ifdef STREAM_MODE
	for (i = 0; i < 102; i++)
	{
		c = (r[i/8] >> (i%8)) & 1;
		trigger_high();
		e = eval(f, L[i]);

		e_inv = gf_inv(gf_mul(e,e));

		for (j = 0; j < 2*SYS_T ; j++)
		{
			out[j] = gf_add(out[j], gf_mul(e_inv, c));
			e_inv = gf_mul(e_inv, L[i]);
		}
	}	
	#else
	for (i = 0; i < 1; i++)
	{
		c = (r[i/8] >> (i%8)) & 1;
		trigger_high();
		e = eval(f, L[i]);

		e_inv = gf_inv(gf_mul(e,e));

		for (j = 0; j < 28 ; j++)
		{
			out[j] = gf_add(out[j], gf_mul(e_inv, c));
			e_inv = gf_mul(e_inv, L[i]);
		}
	}
	#endif
}