/*
###############################################################################
# BRW256 developers and authors:                                              #                                                                           
# Sebati Ghosh,        		Indian Statistical Institute                  #
# Palash Sarkar,       		Indian Statistical Institute                  #
###############################################################################
#                                                                             #
###############################################################################
#                                                                             #
# Copyright (c) 2017, Sebati Ghosh, Palash Sarkar,                            #
#                                                                             #
# Permission to use this code for BRW256 is granted.                          #
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


#ifndef _MULT256_H_
#define _MULT256_H_

#define MAX1 13


	#define XOR(var1,var2,var3){\
		var3[0] = _mm_xor_si128(var1[0],var2[0]);\
		var3[1] = _mm_xor_si128(var1[1],var2[1]);}
	
	
	#define XORNEW(var1,var2,var3){\
		var3[0] = _mm_xor_si128(var1[0],var2[0]);\
		var3[1] = _mm_xor_si128(var1[1],var2[1]);\
		var3[2] = _mm_xor_si128(var1[2],var2[2]);\
		var3[3] = _mm_xor_si128(var1[3],var2[3]);}

	#define ASSIGN(var1,var2){\
		var1[0] = var2[0];\
		var1[1] = var2[1];}


	#define ASSIGNNEW(var1,var2){\
		var1[0] = var2[0];\
		var1[1] = var2[1];\
		var1[2] = var2[2];\
		var1[3] = var2[3];}


	__m128i keyPow[MAX1][2];
	
	



	#define karatsuba(aa,bb,cc)  ({__m128i l1[2],l2[2],l3[2],jnew1[2],j2[2],j3[2],xx[2],ee[2];\
							   l1[0] = _mm_srli_si128 (aa[0], 8);\
							   l1[0] = _mm_xor_si128(l1[0],aa[0]);\
							   jnew1[0] = _mm_srli_si128 (bb[0], 8);\
							   jnew1[0] = _mm_xor_si128(jnew1[0],bb[0]);\
							   l2[0] = _mm_srli_si128 (aa[1], 8);\
							   l2[0] = _mm_xor_si128(l2[0],aa[1]);\
							   j2[0] = _mm_srli_si128 (bb[1], 8);\
							   j2[0] = _mm_xor_si128(j2[0],bb[1]);\
							   xx[0] = _mm_xor_si128(aa[0],aa[1]);\
							   xx[1] = _mm_xor_si128(bb[0],bb[1]);\
							   l3[0] = _mm_srli_si128 (xx[0], 8);\
							   l3[0] = _mm_xor_si128(l3[0],xx[0]);\
							   j3[0] = _mm_srli_si128 (xx[1], 8);\
							   j3[0] = _mm_xor_si128(j3[0],xx[1]);\
							   cc[0] =_mm_clmulepi64_si128( aa[0], bb[0], 0x00 );\
							   cc[1] =_mm_clmulepi64_si128( aa[0], bb[0], 0x11 );\
							   l1[0] = _mm_clmulepi64_si128( l1[0], jnew1[0], 0x00 );\
							   cc[2] =_mm_clmulepi64_si128( aa[1], bb[1], 0x00 );\
							   cc[3] =_mm_clmulepi64_si128( aa[1], bb[1], 0x11 );\
							   l2[0] = _mm_clmulepi64_si128( l2[0], j2[0], 0x00 );\
							   ee[0] =_mm_clmulepi64_si128( xx[0], xx[1], 0x00 );\
							   ee[1] =_mm_clmulepi64_si128( xx[0], xx[1], 0x11 );\
							   l3[0] = _mm_clmulepi64_si128( l3[0], j3[0], 0x00 );\
							   jnew1[0] = _mm_xor_si128(cc[0],cc[1]);\
							   /*cc =_mm256_loadu2_m128i (__m128i const* hiaddr, __m128i const* loaddr)*/\
							   l1[0] = _mm_xor_si128(l1[0],jnew1[0]);\
							   jnew1[0] = _mm_slli_si128 (l1[0], 8);\
							   l1[0] = _mm_srli_si128 (l1[0], 8);\
							   cc[0] = _mm_xor_si128(cc[0],jnew1[0]);\
							   cc[1] = _mm_xor_si128(cc[1],l1[0]);\
							   j2[0] = _mm_xor_si128(cc[2],cc[3]);\
							   /*cc =_mm256_loadu2_m128i (__m128i const* hiaddr, __m128i const* loaddr)*/\
							   l2[0] = _mm_xor_si128(l2[0],j2[0]);\
							   j2[0] = _mm_slli_si128 (l2[0], 8);\
							   l2[0] = _mm_srli_si128 (l2[0], 8);\
							   cc[2] = _mm_xor_si128(cc[2],j2[0]);\
							   cc[3] = _mm_xor_si128(cc[3],l2[0]);\
							   j3[0] = _mm_xor_si128(ee[0],ee[1]);\
							   /*cc =_mm256_loadu2_m128i (__m128i const* hiaddr, __m128i const* loaddr)*/\
							   l3[0] = _mm_xor_si128(l3[0],j3[0]);\
							   j3[0] = _mm_slli_si128 (l3[0], 8);\
							   l3[0] = _mm_srli_si128 (l3[0], 8);\
							   ee[0] = _mm_xor_si128(ee[0],j3[0]);\
							   ee[1] = _mm_xor_si128(ee[1],l3[0]);\
							   ee[0] = _mm_xor_si128(ee[0],cc[0]);\
							   ee[0] = _mm_xor_si128(ee[0],cc[2]);\
							   ee[1] = _mm_xor_si128(ee[1],cc[1]);\
							   ee[1] = _mm_xor_si128(ee[1],cc[3]);\
							   cc[1] = _mm_xor_si128(ee[0],cc[1]);\
							   cc[2] = _mm_xor_si128(ee[1],cc[2]);})






	


	#define reductionbymult(cc, dd) ({__m128i l1[2],jnew1[2],j2[2],j3[2];\
	jnew1[0] = _mm_clmulepi64_si128( cc[3], g1, 0x01 );\
	j2[0] = _mm_clmulepi64_si128( cc[3], g1, 0x00 );\
	j3[0] = _mm_slli_si128 (jnew1[0], 8);\
	jnew1[0] = _mm_srli_si128 (jnew1[0], 8);\
	jnew1[0] = _mm_xor_si128(cc[2],jnew1[0]);\
	j2[0] = _mm_xor_si128(j2[0],j3[0]);\
	j3[0] = _mm_clmulepi64_si128( g1, jnew1[0], 0x00 );\
	jnew1[0] = _mm_clmulepi64_si128( g1, jnew1[0], 0x10 );\
	l1[0] = _mm_srli_si128 (jnew1[0], 8);\
	jnew1[0] = _mm_slli_si128 (jnew1[0], 8);\
	dd[0] = _mm_xor_si128(j3[0],jnew1[0]);\
	dd[0] = _mm_xor_si128(cc[0],dd[0]);\
	dd[1] = _mm_xor_si128(l1[0],j2[0]);\
	dd[1] = _mm_xor_si128(dd[1],cc[1]);})

	

	#define mult(aa,bb,yy) ({__m128i cc[4];\
			karatsuba(aa,bb,cc);\
			reductionbymult(cc,yy);})


	







void partialone(__m128i a[][2],__m128i dig[4])
	{
		//ASSIGN(dig,a[0]);
		
		dig[0] = a[0][0];
		dig[1] = a[0][1];
		dig[2] = dig[3] = _mm_setzero_si128 ();
	}

void partialtwo(__m128i a[][2],__m128i dig[4])
	{
  		__m128i multpart1[4],g1;
		g1 = _mm_set_epi32(0x0,0x0,0x0,0x00000425);
		karatsuba(a[0],keyPow[0],dig);
		XOR(a[1],dig,dig);
		//dig[2] = dig[3] = _mm_setzero_si128 ();

	}

void partialthree(__m128i a[][2],__m128i dig[4])
	{
  		__m128i multpart1[4],g1;
		__m128i btmp1[2],btmp2[2];
		g1 = _mm_set_epi32(0x0,0x0,0x0,0x00000425);
		XOR(a[0],keyPow[0],btmp1);
		XOR(a[1],keyPow[1],btmp2);
		karatsuba(btmp1,btmp2,dig);
		XOR(a[2],dig,dig);
		//dig[2] = dig[3] = _mm_setzero_si128 ();
		//reductionbymult(multpart1,dig);


	}

void partialfour(__m128i a[][2],__m128i dig[4])
	{
  		__m128i multpart1[4],g1;
		__m128i btmp1[2],btmp2[2];
		g1 = _mm_set_epi32(0x0,0x0,0x0,0x00000425);
		XOR(a[0],keyPow[0],btmp1);
		XOR(a[1],keyPow[1],btmp2);
		karatsuba(btmp1,btmp2,multpart1);
		XOR(a[2],multpart1,multpart1);
		reductionbymult(multpart1,btmp1);
		XOR(a[3],keyPow[2],btmp2);
		karatsuba(btmp1,btmp2,dig);

	}

void partialfive(__m128i a[][2],__m128i dig[4])
	{
  		
		__m128i multpart1[4],g1;
		__m128i btmp1[2],btmp2[2];
		g1 = _mm_set_epi32(0x0,0x0,0x0,0x00000425);
		XOR(a[0],keyPow[0],btmp1);
		XOR(a[1],keyPow[1],btmp2);
		karatsuba(btmp1,btmp2,multpart1);
		XOR(a[2],multpart1,multpart1);
		reductionbymult(multpart1,btmp1);
		XOR(a[3],keyPow[2],btmp2);
		karatsuba(btmp1,btmp2,dig);
		XOR(a[4],dig,dig);
		//reductionbymult(multpart1,dig);

	}


void partialsix(__m128i a[][2],__m128i dig[4])
	{
  		
		__m128i multpart1[4],multpart2[4],g1;
		__m128i btmp1[2],btmp2[2];
		g1 = _mm_set_epi32(0x0,0x0,0x0,0x00000425);
		XOR(a[0],keyPow[0],btmp1);
		XOR(a[1],keyPow[1],btmp2);
		karatsuba(btmp1,btmp2,multpart1);
		XOR(a[2],multpart1,multpart1);
		reductionbymult(multpart1,btmp1);
		XOR(a[3],keyPow[2],btmp2);
		
		karatsuba(btmp1,btmp2,multpart1);
		karatsuba(a[4],keyPow[0],dig);
		XOR(a[5],dig,dig);
		XORNEW(multpart1,dig,dig);
		//reductionbymult(multpart1,dig);
	}



void partialseven(__m128i a[][2],__m128i dig[4])
	{
  		
		__m128i multpart1[4],multpart2[4],g1;
		__m128i btmp1[2],btmp2[2],btmp9[2],btmp10[2];
		g1 = _mm_set_epi32(0x0,0x0,0x0,0x00000425);				
		XOR(a[0],keyPow[0],btmp1);
		XOR(a[1],keyPow[1],btmp2);
		karatsuba(btmp1,btmp2,multpart1);
		XOR(a[2],multpart1,multpart1);
		reductionbymult(multpart1,btmp1);
		XOR(a[3],keyPow[2],btmp2);
		XOR(a[4],keyPow[0],btmp9);
		XOR(a[5],keyPow[1],btmp10);
		
		karatsuba(btmp1,btmp2,multpart1);
		karatsuba(btmp9,btmp10,dig);
		XOR(a[6],dig,dig);
		XORNEW(multpart1,dig,dig);
		//reductionbymult(multpart1,dig);

	}

#endif
