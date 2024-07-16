#include "hal.h"
#include "simpleserial.h"

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

#include "util.h"
#include "synd.h"
#include "params.h"
#include "root.h"

unsigned char r[ SYS_N/8 ] = {0};

gf g[ SYS_T+1 ];
gf L[ SYS_N ];

gf s[ SYS_T*2 ];

uint8_t store_received_word(uint8_t* received_word, uint8_t len)
{
  /* The 34-byte pk array comprises the following fields:
     - 1 byte   : offset
     - 1 byte   : number of values
     - 32 bytes : values to store in the received word r
     At most 34 bytes can be sent, expect serial problems otherwise
  */
  uint8_t OK[32] = {0};
  uint8_t offset = received_word[0];
  uint8_t nb_values = received_word[1];
  for (int val = 0; val<nb_values; val++) {
    r[32*offset+val] = received_word[2+val];
    OK[val] = r[32*offset+val];
  }
  #ifdef STREAM_MODE
  simpleserial_put('r', 32, OK);
  #else 
  simpleserial_put('r', 1, OK);
  #endif
  return 0x00;
}

uint8_t store_goppa_polynomial(uint8_t* goppa_polynomial, uint8_t len)
{
  /* The 34-byte pk array comprises the following fields:
     - 1 byte   : offset
     - 1 byte   : number of values
     - 32 bytes : values to store in the goppa polynomial g
     At most 34 bytes can be sent, expect serial problems otherwise
  */
  uint8_t OK[32] = {0};
  uint8_t offset = goppa_polynomial[0];
  uint8_t nb_values = goppa_polynomial[1]/sizeof(gf);
  for (int val = 0; val<nb_values; val++) {
  	g[16*offset+val] = load_gf(goppa_polynomial+2*(val+1));
  	store_gf(OK + 2*val, g[16*offset+val]);
  }
  simpleserial_put('g', 32, OK);
  return 0x00;
}

uint8_t store_support(uint8_t* support, uint8_t len)
{
  /* The 34-byte pk array comprises the following fields:
     - 1 byte   : offset
     - 1 byte   : number of values
     - 32 bytes : values to store in the support L
     At most 34 bytes can be sent, expect serial problems otherwise
  */
  uint8_t OK[32] = {0};
  uint8_t offset = support[0];
  uint8_t nb_values = support[1]/sizeof(gf);
  for (int val = 0; val<nb_values; val++) {
  	L[16*offset+val] = load_gf(support+2*(val+1));
  	store_gf(OK + 2*val, L[16*offset+val]);
  }
  #ifdef STREAM_MODE
  simpleserial_put('L', 32, OK);
  #else
  simpleserial_put('L', 2, OK);
  #endif
  return 0x00;
}

uint8_t compute_synd(uint8_t* msg, uint8_t len)
{
  synd(s, g, L, r);
  trigger_low();
  #ifdef CHECK_SYND
  for (int offset = 0; offset<(((SYS_T*4)+7)/32); offset++)
  {
    simpleserial_put('s', 32, (uint8_t*)(s+16*offset));
  }
  #else
  for (int offset = 0; offset<1; offset++)
  {
  	simpleserial_put('s', 1, (uint8_t*)(s+16*offset));
  }
  #endif
  return 0x00;
}

int main(void)
{
	platform_init();
	init_uart();
	trigger_setup();

	simpleserial_init();
  #ifdef STREAM_MODE
  simpleserial_addcmd('r', 2+32, store_received_word); 	                //Store the received word r
  #else
  simpleserial_addcmd('r', 2+1, store_received_word);                   //Store the first byte of the received word r
  #endif
	simpleserial_addcmd('g', 2+32, store_goppa_polynomial); 		          //Store the goppa polynomial g
  #ifdef STREAM_MODE
	simpleserial_addcmd('L', 2+32, store_support);	                      //Store the support L
  #else
	simpleserial_addcmd('L', 2+2, store_support);                         //Store the first element of the support L
  #endif
	simpleserial_addcmd('s', 0, compute_synd);              		          //Compute the syndrome s
	while(1)
		simpleserial_get();
}