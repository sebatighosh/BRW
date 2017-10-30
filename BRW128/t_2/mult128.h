/*
###############################################################################
# BRW128 developers and authors:                                              #                                                                           
# Sebati Ghosh,        		Indian Statistical Institute                  #
# Palash Sarkar,       		Indian Statistical Institute                  #
###############################################################################
#                                                                             #
###############################################################################
#                                                                             #
# Copyright (c) 2017, Sebati Ghosh, Palash Sarkar,                            #
#                                                                             #
# Permission to use this code for BRW128 is granted.                          #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are      #
# met:                                                                        #
#                                                                             #
# * Redistributions of source code must retain the above copyright notice,    #
#   this list of conditions and the following disclaimer.                     #
#                                                                             #
# * Redistributions in binary form must reproduce the above copyright         #
#   notice, this list of conditions and the following disclaimer in the       #
#   documentation and/or other materials provided with the distribution.      #
#                                                                             #
# * The names of the contributors may not be used to endorse or promote       #
# products derived from this software without specific prior written          #
# permission.                                                                 #
#                                                                             #
###############################################################################
#                                                                             #
###############################################################################
# THIS SOFTWARE IS PROVIDED BY THE AUTHORS ""AS IS"" AND ANY                  #
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE           #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR          #
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE CONTRIBUTORS BE LIABLE        #
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 		      #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT 	      #
# OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR          	      #
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF      #
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING        #
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS          #
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                #
###############################################################################
*/



#ifndef MULT_H_
#define MULT_H_

#define MAX1 15

	#define XOR(var1,var2,var3){\
		var3 = _mm_xor_si128(var1,var2);}

	#define XORNEW(var1,var2,var3){\
		var3[0] = _mm_xor_si128(var1[0],var2[0]);\
		var3[1] = _mm_xor_si128(var1[1],var2[1]);}

	#define ASSIGN(var1,var2){\
		var1 = var2;}
	#define ASSIGNNEW(var1,var2){\
		var1[0] = var2[0];\
		var1[1] = var2[1];}


	
	
	__m128i keyPow[MAX1];
	
	
	
	#define schoolbook(aa,bb,cc,xx)  ({\
	__m128i r,s;\
	cc[0] =_mm_clmulepi64_si128( aa, bb, 0x00 );\
	cc[1] =_mm_clmulepi64_si128( aa, bb, 0x11 );\
	r = _mm_clmulepi64_si128( aa, bb, 0x01 );\
	s = _mm_clmulepi64_si128( aa, bb, 0x10 );\
	xx = _mm_xor_si128(r,s);\
	/*p = _mm_slli_si128 (r, 8);\
	s = _mm_srli_si128 (r, 8);\
	cc[0] = _mm_xor_si128(cc[0],p);\
	cc[1] = _mm_xor_si128(cc[1],s);*/})


	#define reductionbymult(cc,xx,dd)   ({\
	__m128i r,s,p;\
	p = _mm_slli_si128 (xx, 8);\
	s = _mm_srli_si128 (xx, 8);\
	cc[0] = _mm_xor_si128(cc[0],p);\
	cc[1] = _mm_xor_si128(cc[1],s);\
    	r = _mm_clmulepi64_si128( cc[1], POLY, 0x01 );\
	p = _mm_srli_si128 (r, 8);\
	p = _mm_xor_si128(p,cc[1]);\
	s =_mm_clmulepi64_si128( p, POLY, 0x00 );\
	p = _mm_slli_si128 (r, 8);\
	dd = _mm_xor_si128(s,p);\
	dd = _mm_xor_si128(dd,cc[0]);})

	

	#define mult(aa,bb,dd) ({__m128i cc[2],xx;\
	schoolbook(aa,bb,cc,xx);\
	reductionbymult(cc,xx,dd);})


	








void partialone(__m128i a[],__m128i dig[2], __m128i *mid)
	{
		ASSIGN(*dig,a[0]);
		dig[1] = _mm_setzero_si128 ();
		*mid = _mm_setzero_si128 ();
		
	}


void partialtwo(__m128i a[],__m128i dig[2], __m128i *mid)
	{
  		__m128i POLY;
		
		POLY = _mm_set_epi32(0x0,0x0,0x0,0x00000087);
		schoolbook(a[0],keyPow[0],dig,*mid);
	        XOR(a[1],dig[0],dig[0]);
		//printf("in");
		//printf("%x",_mm_extract_epi8 (dig[0], 0));
		//reductionbymult(multpart1,xx,*dig);


	}



void partialthree(__m128i a[],__m128i dig[2], __m128i *mid)
	{
  		__m128i POLY;
		__m128i btmp1,btmp2;
		POLY = _mm_set_epi32(0x0,0x0,0x0,0x00000087);
		XOR(a[0],keyPow[0],btmp1);
		XOR(a[1],keyPow[1],btmp2);
		schoolbook(btmp1,btmp2,dig,*mid);
		XOR(a[2],dig[0],dig[0]);
		//reductionbymult(multpart1,xx,*dig);
	}

	


	


#endif
