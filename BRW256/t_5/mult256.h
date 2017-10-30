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


	



#define evalbrwseq(aa,dig) ({__m128i new1[4], new2[4], new26[4], new10[4];\
__m128i btmp1[2],btmp2[2],btmp9[2],btmp10[2],btmp17[2],btmp18[2],btmp25[2],btmp26[2];\
	XOR(aa[0],keyPow[0],btmp1);\
	XOR(aa[1],keyPow[1],btmp2);\
	/*XOR(aa[4],key,btmp5);\
	XOR(aa[5],keyPow[1],btmp6);\
	/*XOR(aa[7],keyPow[3],btmp8);*/\
	XOR(aa[8],keyPow[0],btmp9);\
	XOR(aa[9],keyPow[1],btmp10);\
	XOR(aa[16],keyPow[0],btmp17);\
	XOR(aa[17],keyPow[1],btmp18);\
	karatsuba(btmp2,btmp1,new1);\
	karatsuba(btmp10,btmp9,new2);\
	karatsuba(btmp18,btmp17,new26);\
	XOR(aa[2],new1,new1);\
	XOR(aa[10],new2,new2);\
	XOR(aa[18],new26,new26);\
	reductionbymult(new1,btmp2);\
	reductionbymult(new2,btmp10);\
	reductionbymult(new26,btmp18);\
	XOR(aa[4],keyPow[0],btmp1);\
	XOR(aa[5],keyPow[1],btmp9);\
	XOR(aa[3],keyPow[2],btmp17);\
	XOR(aa[24],keyPow[0],btmp25);\
	XOR(aa[25],keyPow[1],btmp26);\
	karatsuba(btmp1,btmp9,new1);\
	karatsuba(btmp2,btmp17,new2);\
	karatsuba(btmp26,btmp25,new26);\
	XOR(aa[26],new26,new26);\
	reductionbymult(new26,btmp26);\
	XOR(aa[6],new1,new1);\
	XORNEW(new1,new2,new1);\
	/*XOR(new1[1],new2[1],new1[1]);*/\
	reductionbymult(new1,btmp1);\
	XOR(aa[7],keyPow[3],btmp9);\
	XOR(aa[11],keyPow[2],btmp2);\
	XOR(aa[12],keyPow[0],btmp17);\
	XOR(aa[13],keyPow[1],btmp25);\
	/*XOR(aa[11],keyPow[2],btmp12);*/\
	karatsuba(btmp1,btmp9,new1);\
	karatsuba(btmp2,btmp10,new2);\
	karatsuba(btmp25,btmp17,new26);\
	XOR(aa[14],new26,new26);\
	XORNEW(new2,new26,new26);\
	/*XOR(new2[1],new26[1],new26[1]);*/\
	XORNEW(new1,new26,new26);\
	/*XOR(new1[1],new26[1],new26[1]);*/\
	reductionbymult(new26,btmp25);\
	XOR(aa[15],keyPow[4],btmp1);\
	/*XOR(aa[16],key,btmp17);\
	XOR(aa[17],keyPow[1],btmp18);*/\
	XOR(aa[19],keyPow[2],btmp2);\
	XOR(aa[20],keyPow[0],btmp9);\
	XOR(aa[21],keyPow[1],btmp17);\
	karatsuba(btmp1,btmp25,new1);\
	karatsuba(btmp2,btmp18,new2);\
	karatsuba(btmp9,btmp17,new26);\
	XOR(aa[22],new26,new26);\
	XORNEW(new2,new26,new2);\
	/*XOR(new2[1],new26[1],new2[1]);*/\
	reductionbymult(new2,btmp2);\
	XOR(aa[23],keyPow[3],btmp9);\
	XOR(aa[27],keyPow[2],btmp17);\
	XOR(aa[28],keyPow[0],btmp1);\
	XOR(aa[29],keyPow[1],btmp10);\
	karatsuba(btmp2,btmp9,new2);\
	karatsuba(btmp26,btmp17,new26);\
	karatsuba(btmp10,btmp1,new10);\
	XOR(aa[30],new10,new10);\
	XORNEW(new10,new26,new10);\
	/*XOR(new10[1],new26[1],new10[1]);*/\
	XORNEW(new10,new2,new10);\
	/*XOR(new10[1],new2[1],new10[1]);*/\
	XORNEW(new10,new1,dig);\
	/*XOR(new10[1],new1[1],new10[1]);*/\
	/*reductionbymult(new10, dig);*/\
	})







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

void partialeight(__m128i a[][2],__m128i dig[4])
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
		karatsuba(btmp9,btmp10,multpart2);
		XOR(a[6],multpart2,multpart2);
		XORNEW(multpart1,multpart2,multpart1);
		reductionbymult(multpart1,btmp1);
		XOR(a[7],keyPow[3],btmp2);
		karatsuba(btmp1,btmp2,dig);

	}


void partialnine(__m128i a[][2],__m128i dig[4])
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
		karatsuba(btmp9,btmp10,multpart2);
		XOR(a[6],multpart2,multpart2);
		XORNEW(multpart1,multpart2,multpart1);
		reductionbymult(multpart1,btmp1);
		XOR(a[7],keyPow[3],btmp2);
		karatsuba(btmp1,btmp2,dig);
		XOR(a[8],dig,dig);
		
		//reductionbymult(multpart1,dig);
	}


void partialten(__m128i a[][2],__m128i dig[4])
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
		karatsuba(btmp9,btmp10,multpart2);
		XOR(a[6],multpart2,multpart2);
		XORNEW(multpart1,multpart2,multpart1);
		reductionbymult(multpart1,btmp1);
		XOR(a[7],keyPow[3],btmp2);
		
		karatsuba(btmp1,btmp2,multpart1);
		karatsuba(a[8],keyPow[0],multpart2);
		XOR(a[9],multpart2,multpart2);
		XORNEW(multpart1,multpart2,dig);
		//reductionbymult(multpart1,dig);
	}

void partialeleven(__m128i a[][2],__m128i dig[4])
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
		karatsuba(btmp9,btmp10,multpart2);
		XOR(a[6],multpart2,multpart2);
		XORNEW(multpart1,multpart2,multpart1);
		reductionbymult(multpart1,btmp1);
		XOR(a[7],keyPow[3],btmp2);
		XOR(a[8],keyPow[0],btmp9);
		XOR(a[9],keyPow[1],btmp10);
		
		karatsuba(btmp1,btmp2,multpart1);
		karatsuba(btmp9,btmp10,multpart2);
		XOR(a[10],multpart2,multpart2);
		XORNEW(multpart1,multpart2,dig);
		//reductionbymult(multpart1,dig);

	}


void partialtwelve(__m128i a[][2],__m128i dig[4])
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
		karatsuba(btmp9,btmp10,multpart2);
		XOR(a[6],multpart2,multpart2);
		XORNEW(multpart1,multpart2,multpart1);
		reductionbymult(multpart1,btmp1);
		XOR(a[7],keyPow[3],btmp2);
		XOR(a[8],keyPow[0],btmp9);
		XOR(a[9],keyPow[1],btmp10);
	
		karatsuba(btmp1,btmp2,multpart1);
		karatsuba(btmp9,btmp10,multpart2);
		XOR(a[10],multpart2,multpart2);
		reductionbymult(multpart2,btmp9);
		XOR(a[11],keyPow[2],btmp10);
		karatsuba(btmp9,btmp10,multpart2);
		XORNEW(multpart1,multpart2,dig);
		//reductionbymult(multpart1,dig);
	}

void partialthirteen(__m128i a[][2],__m128i dig[4])
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
		karatsuba(btmp9,btmp10,multpart2);
		XOR(a[6],multpart2,multpart2);
		XORNEW(multpart1,multpart2,multpart1);
		reductionbymult(multpart1,btmp1);
		XOR(a[7],keyPow[3],btmp2);
		XOR(a[8],keyPow[0],btmp9);
		XOR(a[9],keyPow[1],btmp10);
		
		karatsuba(btmp1,btmp2,multpart1);
		karatsuba(btmp9,btmp10,multpart2);
		XOR(a[10],multpart2,multpart2);
		reductionbymult(multpart2,btmp9);
		XOR(a[11],keyPow[2],btmp10);
		karatsuba(btmp9,btmp10,multpart2);
		XOR(a[12],multpart2,multpart2);
		XORNEW(multpart1,multpart2,dig);
		//reductionbymult(multpart1,dig);
	}

void partialfourteen(__m128i a[][2],__m128i dig[4])
	{
  		
		__m128i multpart1[4],multpart2[4],multpart3[4],g1;
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
		karatsuba(btmp9,btmp10,multpart2);
		XOR(a[6],multpart2,multpart2);
		XORNEW(multpart1,multpart2,multpart1);
		reductionbymult(multpart1,btmp1);
		XOR(a[7],keyPow[3],btmp2);
		XOR(a[8],keyPow[0],btmp9);
		XOR(a[9],keyPow[1],btmp10);
		
		karatsuba(btmp1,btmp2,multpart1);
		karatsuba(btmp9,btmp10,multpart2);
		XOR(a[10],multpart2,multpart2);
		reductionbymult(multpart2,btmp9);
		XOR(a[11],keyPow[2],btmp10);
		
		karatsuba(btmp9,btmp10,multpart2);
		karatsuba(a[12],keyPow[0],multpart3);
		XOR(a[13],multpart3,multpart3);
		XORNEW(multpart1,multpart2,multpart1);
		XORNEW(multpart1,multpart3,dig);
		//reductionbymult(multpart1,dig);

	}

void partialfifteen(__m128i a[][2],__m128i dig[4])
	{
  		
		__m128i multpart1[4],multpart2[4],multpart3[4],g1;
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
		karatsuba(btmp9,btmp10,multpart2);
		XOR(a[6],multpart2,multpart2);
		XORNEW(multpart1,multpart2,multpart1);
		reductionbymult(multpart1,btmp1);
		XOR(a[7],keyPow[3],btmp2);
		XOR(a[8],keyPow[0],btmp9);
		XOR(a[9],keyPow[1],btmp10);
		
		karatsuba(btmp1,btmp2,multpart1);
		karatsuba(btmp9,btmp10,multpart2);
		XOR(a[10],multpart2,multpart2);
		reductionbymult(multpart2,btmp9);
		XOR(a[11],keyPow[2],btmp10);
		XOR(a[12],keyPow[0],btmp1);
		XOR(a[13],keyPow[1],btmp2);
		
		karatsuba(btmp9,btmp10,multpart2);
		karatsuba(btmp1,btmp2,multpart3);
		XOR(a[14],multpart3,multpart3);
		XORNEW(multpart1,multpart2,multpart1);
		XORNEW(multpart1,multpart3,dig);
		//reductionbymult(multpart1,dig);

		
	}


void partialrest(__m128i a[][2],int noOfBlock,__m128i dig[4], __m128i dig1[4])
	{
		
		__m128i tmpdig[4],g1, tmp[2];
		__m128i btmp1[2];
		g1 = _mm_set_epi32(0x0,0x0,0x0,0x00000425);

		
		partialfifteen(a,dig);
		reductionbymult(dig,tmp);
		
		XOR(a[15],keyPow[4],btmp1);  		
		
		karatsuba(tmp,btmp1,dig);
		//reductionbymult(tmpdig,dig);
		if(noOfBlock == 16)
			return;
		if(noOfBlock == 17)
		{
			XOR(a[16],dig,dig);
			//reductionbymult(tmpdig,dig);
			return;

		}
		switch(noOfBlock)
		{
		/*case 17:
			partialone(a[16],key,dig1);
			break;*/
		case 18:
			partialtwo(a[16],dig1);
			break;
		case 19:
			partialthree(a[16],dig1);
			break;
		case 20:
			partialfour(a[16],dig1);
			break;
		case 21:
			partialfive(a[16],dig1);
			break;
		case 22:
			partialsix(a[16],dig1);
			break;
		case 23:
			partialseven(a[16],dig1);
			break;
		case 24:
			partialeight(a[16],dig1);
			break;
		case 25:
			partialnine(a[16],dig1);
			break;
		case 26:
			partialten(a[16],dig1);
			break;
		case 27:
			partialeleven(a[16],dig1);
			break;
		case 28:
			partialtwelve(a[16],dig1);
			break;
		case 29:
			partialthirteen(a[16],dig1);
			break;
		case 30:
			partialfourteen(a[16],dig1);
			break;
		}
		XORNEW(dig, dig1, dig);
	}






#endif
