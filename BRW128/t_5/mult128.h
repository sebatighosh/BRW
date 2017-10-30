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


	




#define evalbrwseq(aa,cc,yy) ({__m128i new1[2], new2[2], new26[2], new10[2];\
	__m128i btmp1,btmp2,btmp9,btmp10,btmp17,btmp18,btmp25,btmp26,xx,xx1,xx2,xx3;\
/* new19[2], new23[2], new27[2], new31[2], new32[2], new34[2], new33[2], new39[2], new35[2], new40[2], new42[2], sum32[2], sum37[2], sum38[2], sum41[2], sum42[2];*/\
	XOR(aa[0],keyPow[0],btmp1);\
	XOR(aa[1],keyPow[1],btmp2);\
	/*XOR(aa[4],keyPow[0],btmp5);\
	XOR(aa[5],keyPow[1],btmp6);\
	/*XOR(aa[7],keyPow[3],btmp8);*/\
	XOR(aa[8],keyPow[0],btmp9);\
	XOR(aa[9],keyPow[1],btmp10);\
	XOR(aa[16],keyPow[0],btmp17);\
	XOR(aa[17],keyPow[1],btmp18);\
	schoolbook(btmp2,btmp1,new1,xx);\
	schoolbook(btmp10,btmp9,new2,xx1);\
	schoolbook(btmp18,btmp17,new26,xx2);\
	XOR(aa[2],new1[0],new1[0]);\
	XOR(aa[10],new2[0],new2[0]);\
	XOR(aa[18],new26[0],new26[0]);\
	reductionbymult(new1,xx,btmp2);\
	reductionbymult(new2,xx1,btmp10);\
	reductionbymult(new26,xx2,btmp18);\
	XOR(aa[4],keyPow[0],btmp1);\
	XOR(aa[5],keyPow[1],btmp9);\
	XOR(aa[3],keyPow[2],btmp17);\
	XOR(aa[24],keyPow[0],btmp25);\
	XOR(aa[25],keyPow[1],btmp26);\
	schoolbook(btmp1,btmp9,new1,xx);\
	schoolbook(btmp2,btmp17,new2,xx1);\
	schoolbook(btmp26,btmp25,new26,xx2);\
	XOR(aa[26],new26[0],new26[0]);\
	reductionbymult(new26,xx2,btmp26);\
	XOR(aa[6],new1[0],new1[0]);\
	/*XOR(new1[0],new2[0],new1[0]);\
	XOR(new1[1],new2[1],new1[1]);*/\
	XORNEW(new1,new2,new1);\
	XOR(xx,xx1,xx);\
	reductionbymult(new1,xx,btmp1);\
	XOR(aa[7],keyPow[3],btmp9);\
	XOR(aa[11],keyPow[2],btmp2);\
	XOR(aa[12],keyPow[0],btmp17);\
	XOR(aa[13],keyPow[1],btmp25);\
	/*XOR(aa[11],keyPow[2],btmp12);*/\
	schoolbook(btmp1,btmp9,new1,xx);\
	schoolbook(btmp2,btmp10,new2,xx1);\
	schoolbook(btmp25,btmp17,new26,xx2);\
	XOR(aa[14],new26[0],new26[0]);\
	/*XOR(new2[0],new26[0],new26[0]);\
	XOR(new2[1],new26[1],new26[1]);*/\
	XORNEW(new2,new26,new26);\
	XOR(xx1,xx2,xx2);\
	/*XOR(new1[0],new26[0],new26[0]);\
	XOR(new1[1],new26[1],new26[1]);*/\
	XORNEW(new1,new26,new26);\
	XOR(xx,xx2,xx2);\
	reductionbymult(new26,xx2,btmp25);\
	XOR(aa[15],keyPow[4],btmp1);\
	/*XOR(aa[16],keyPow[0],btmp17);\
	XOR(aa[17],keyPow[1],btmp18);*/\
	XOR(aa[19],keyPow[2],btmp2);\
	XOR(aa[20],keyPow[0],btmp9);\
	XOR(aa[21],keyPow[1],btmp17);\
	schoolbook(btmp1,btmp25,new1,xx);\
	schoolbook(btmp2,btmp18,new2,xx1);\
	schoolbook(btmp9,btmp17,new26,xx2);\
	XOR(aa[22],new26[0],new26[0]);\
	/*XOR(new2[0],new26[0],new2[0]);\
	XOR(new2[1],new26[1],new2[1]);*/\
	XORNEW(new2,new26,new2);\
	XOR(xx1,xx2,xx1);\
	reductionbymult(new2,xx1,btmp2);\
	XOR(aa[23],keyPow[3],btmp9);\
	XOR(aa[27],keyPow[2],btmp17);\
	XOR(aa[28],keyPow[0],btmp1);\
	XOR(aa[29],keyPow[1],btmp10);\
	schoolbook(btmp2,btmp9,new2,xx1);\
	schoolbook(btmp26,btmp17,new26,xx2);\
	schoolbook(btmp10,btmp1,new10,xx3);\
	XOR(aa[30],new10[0],new10[0]);\
	/*XOR(new10[0],new26[0],new10[0]);\
	XOR(new10[1],new26[1],new10[1]);*/\
	XORNEW(new10,new26,new10);\
	XOR(xx3,xx2,xx3);\
	/*XOR(new10[0],new2[0],new10[0]);\
	XOR(new10[1],new2[1],new10[1]);*/\
	XORNEW(new10,new2,new10);\
	XOR(xx3,xx1,xx3);\
	/*XOR(new10[0],new1[0],new10[0]);\
	XOR(new10[1],new1[1],new10[1]);*/\
	XORNEW(new10,new1,cc);\
	XOR(xx3,xx,yy);\
	/*printf("in-cc");\
	printf("%x",_mm_extract_epi8 (cc[0], 0));\
	printf("in-yy");\
	printf("%x",_mm_extract_epi8 (yy, 0));*/\
	/*reductionbymult(new10, xx3, dig);*/})



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

	

void partialfour(__m128i a[],__m128i dig[2], __m128i *mid)
	{
  		__m128i POLY;
		__m128i btmp1,btmp2;
		POLY = _mm_set_epi32(0x0,0x0,0x0,0x00000087);
		XOR(a[0],keyPow[0],btmp1);
		XOR(a[1],keyPow[1],btmp2);
		schoolbook(btmp1,btmp2,dig,*mid);
		XOR(a[2],dig[0],dig[0]);
		reductionbymult(dig,*mid,btmp1);
		XOR(a[3],keyPow[2],btmp2);
		schoolbook(btmp1,btmp2,dig,*mid);

	}



void partialfive(__m128i a[],__m128i dig[2], __m128i *mid)
	{
  		
		__m128i multpart1[2], POLY;
		__m128i btmp1,btmp2,btmp25,btmp26,xx;
		POLY = _mm_set_epi32(0x0,0x0,0x0,0x00000087);
		XOR(a[0],keyPow[0],btmp1);
		XOR(a[1],keyPow[1],btmp2);
		schoolbook(btmp1,btmp2,multpart1,xx);
		XOR(a[2],multpart1[0],multpart1[0]);
		reductionbymult(multpart1,xx,btmp1);
		XOR(a[3],keyPow[2],btmp2);
		schoolbook(btmp1,btmp2,dig,*mid);
		XOR(a[4],dig[0],dig[0]);
		//reductionbymult(multpart1,xx,*dig);

	}


  		
void partialsix(__m128i a[],__m128i dig[2], __m128i *mid)
	{
  		
		__m128i multpart1[2],multpart2[2],POLY;
		__m128i btmp1,btmp2,xx,xx1;
		POLY = _mm_set_epi32(0x0,0x0,0x0,0x00000087);
		XOR(a[0],keyPow[0],btmp1);
		XOR(a[1],keyPow[1],btmp2);
		schoolbook(btmp1,btmp2,multpart1,xx);
		XOR(a[2],multpart1[0],multpart1[0]);
		reductionbymult(multpart1,xx,btmp1);
		XOR(a[3],keyPow[2],btmp2);
		schoolbook(btmp1,btmp2,multpart1,xx);
		schoolbook(a[4],keyPow[0],multpart2,xx1);
		XOR(a[5],multpart2[0],multpart2[0]);
		XORNEW(multpart1,multpart2,dig);
		XOR(xx,xx1,*mid);
		//reductionbymult(multpart1,xx,*dig);
	}


  		

void partialseven(__m128i a[],__m128i dig[2], __m128i *mid)
	{
  		
		__m128i multpart1[2],multpart2[2],POLY;
		__m128i btmp1,btmp2,btmp9,btmp10,xx,xx1;
		POLY = _mm_set_epi32(0x0,0x0,0x0,0x00000087);
		XOR(a[0],keyPow[0],btmp1);
		XOR(a[1],keyPow[1],btmp2);
		schoolbook(btmp1,btmp2,multpart1,xx);
		XOR(a[2],multpart1[0],multpart1[0]);
		reductionbymult(multpart1,xx,btmp1);
		XOR(a[3],keyPow[2],btmp2);
		XOR(a[4],keyPow[0],btmp9);
		XOR(a[5],keyPow[1],btmp10);
		schoolbook(btmp1,btmp2,multpart1,xx);
		schoolbook(btmp9,btmp10,multpart2,xx1);
		XOR(a[6],multpart2[0],multpart2[0]);
		XORNEW(multpart1,multpart2,dig);
		XOR(xx,xx1,*mid);
		//reductionbymult(multpart1,xx,*dig);

	}


void partialeight(__m128i a[],__m128i dig[2], __m128i *mid)
	{
  		
		__m128i multpart1[2],multpart2[2],POLY;
		__m128i btmp1,btmp2,btmp9,btmp10,xx,xx1;
		POLY = _mm_set_epi32(0x0,0x0,0x0,0x00000087);
		XOR(a[0],keyPow[0],btmp1);
		XOR(a[1],keyPow[1],btmp2);
		schoolbook(btmp1,btmp2,multpart1,xx);
		XOR(a[2],multpart1[0],multpart1[0]);
		reductionbymult(multpart1,xx,btmp1);
		XOR(a[3],keyPow[2],btmp2);
		XOR(a[4],keyPow[0],btmp9);
		XOR(a[5],keyPow[1],btmp10);
		schoolbook(btmp1,btmp2,multpart1,xx);
		schoolbook(btmp9,btmp10,multpart2,xx1);
		XOR(a[6],multpart2[0],multpart2[0]);
		XORNEW(multpart1,multpart2,multpart1);
		XOR(xx,xx1,xx);
		reductionbymult(multpart1,xx,btmp1);
		XOR(a[7],keyPow[3],btmp2);
		schoolbook(btmp1,btmp2,dig,*mid);

	}


void partialnine(__m128i a[],__m128i dig[2], __m128i *mid)
	{
  		
		__m128i multpart1[2],multpart2[2],POLY;
		__m128i btmp1,btmp2,btmp9,btmp10,xx,xx1;
		POLY = _mm_set_epi32(0x0,0x0,0x0,0x00000087);
		XOR(a[0],keyPow[0],btmp1);
		XOR(a[1],keyPow[1],btmp2);
		schoolbook(btmp1,btmp2,multpart1,xx);
		XOR(a[2],multpart1[0],multpart1[0]);
		reductionbymult(multpart1,xx,btmp1);
		XOR(a[3],keyPow[2],btmp2);
		XOR(a[4],keyPow[0],btmp9);
		XOR(a[5],keyPow[1],btmp10);
		schoolbook(btmp1,btmp2,multpart1,xx);
		schoolbook(btmp9,btmp10,multpart2,xx1);
		XOR(a[6],multpart2[0],multpart2[0]);
		XORNEW(multpart1,multpart2,multpart1);
		XOR(xx,xx1,xx);
		reductionbymult(multpart1,xx,btmp1);
		XOR(a[7],keyPow[3],btmp2);
		schoolbook(btmp1,btmp2,dig,*mid);
		XOR(a[8],dig[0],dig[0]);
		//reductionbymult(multpart1,xx,*dig);
	}


void partialten(__m128i a[],__m128i dig[2], __m128i *mid)
	{
  		
		__m128i multpart1[2],multpart2[2],POLY;
		__m128i btmp1,btmp2,btmp9,btmp10,xx,xx1;
		POLY = _mm_set_epi32(0x0,0x0,0x0,0x00000087);
		XOR(a[0],keyPow[0],btmp1);
		XOR(a[1],keyPow[1],btmp2);
		schoolbook(btmp1,btmp2,multpart1,xx);
		XOR(a[2],multpart1[0],multpart1[0]);
		reductionbymult(multpart1,xx,btmp1);
		XOR(a[3],keyPow[2],btmp2);
		XOR(a[4],keyPow[0],btmp9);
		XOR(a[5],keyPow[1],btmp10);
		schoolbook(btmp1,btmp2,multpart1,xx);
		schoolbook(btmp9,btmp10,multpart2,xx1);
		XOR(a[6],multpart2[0],multpart2[0]);
		XORNEW(multpart1,multpart2,multpart1);
		XOR(xx,xx1,xx);
		reductionbymult(multpart1,xx,btmp1);
		XOR(a[7],keyPow[3],btmp2);
		schoolbook(btmp1,btmp2,multpart1,xx);
		schoolbook(a[8],keyPow[0],multpart2,xx1);
		XOR(a[9],multpart2[0],multpart2[0]);
		XORNEW(multpart1,multpart2,dig);
		XOR(xx,xx1,*mid);
		//reductionbymult(multpart1,xx,*dig);
	}

void partialeleven(__m128i a[],__m128i dig[2], __m128i *mid)
	{
  		
		__m128i multpart1[2],multpart2[2],POLY;
		__m128i btmp1,btmp2,btmp9,btmp10,xx,xx1;
		POLY = _mm_set_epi32(0x0,0x0,0x0,0x00000087);
		XOR(a[0],keyPow[0],btmp1);
		XOR(a[1],keyPow[1],btmp2);
		schoolbook(btmp1,btmp2,multpart1,xx);
		XOR(a[2],multpart1[0],multpart1[0]);
		reductionbymult(multpart1,xx,btmp1);
		XOR(a[3],keyPow[2],btmp2);
		XOR(a[4],keyPow[0],btmp9);
		XOR(a[5],keyPow[1],btmp10);
		schoolbook(btmp1,btmp2,multpart1,xx);
		schoolbook(btmp9,btmp10,multpart2,xx1);
		XOR(a[6],multpart2[0],multpart2[0]);
		XORNEW(multpart1,multpart2,multpart1);
		XOR(xx,xx1,xx);
		reductionbymult(multpart1,xx,btmp1);
		XOR(a[7],keyPow[3],btmp2);
		XOR(a[8],keyPow[0],btmp9);
		XOR(a[9],keyPow[1],btmp10);
		schoolbook(btmp1,btmp2,multpart1,xx);
		schoolbook(btmp9,btmp10,multpart2,xx1);
		XOR(a[10],multpart2[0],multpart2[0]);
		XORNEW(multpart1,multpart2,dig);
		XOR(xx,xx1,*mid);
		//reductionbymult(multpart1,xx,*dig);

	}


void partialtwelve(__m128i a[],__m128i dig[2], __m128i *mid)
	{
  		
		__m128i multpart1[2],multpart2[2],POLY;
		__m128i btmp1,btmp2,btmp9,btmp10,xx,xx1;
		POLY = _mm_set_epi32(0x0,0x0,0x0,0x00000087);
		XOR(a[0],keyPow[0],btmp1);
		XOR(a[1],keyPow[1],btmp2);
		schoolbook(btmp1,btmp2,multpart1,xx);
		XOR(a[2],multpart1[0],multpart1[0]);
		reductionbymult(multpart1,xx,btmp1);
		XOR(a[3],keyPow[2],btmp2);
		XOR(a[4],keyPow[0],btmp9);
		XOR(a[5],keyPow[1],btmp10);
		schoolbook(btmp1,btmp2,multpart1,xx);
		schoolbook(btmp9,btmp10,multpart2,xx1);
		XOR(a[6],multpart2[0],multpart2[0]);
		XORNEW(multpart1,multpart2,multpart1);
		XOR(xx,xx1,xx);
		reductionbymult(multpart1,xx,btmp1);
		XOR(a[7],keyPow[3],btmp2);
		XOR(a[8],keyPow[0],btmp9);
		XOR(a[9],keyPow[1],btmp10);
		schoolbook(btmp1,btmp2,multpart1,xx);
		schoolbook(btmp9,btmp10,multpart2,xx1);
		XOR(a[10],multpart2[0],multpart2[0]);
		reductionbymult(multpart2,xx1,btmp9);
		XOR(a[11],keyPow[2],btmp10);
		schoolbook(btmp9,btmp10,multpart2,xx1);
		XORNEW(multpart1,multpart2,dig);
		XOR(xx,xx1,*mid);
		//reductionbymult(multpart1,xx,*dig);
	}

void partialthirteen(__m128i a[],__m128i dig[2], __m128i *mid)
	{
  		
		__m128i multpart1[2],multpart2[2],POLY;
		__m128i btmp1,btmp2,btmp9,btmp10,xx,xx1;
		POLY = _mm_set_epi32(0x0,0x0,0x0,0x00000087);
		XOR(a[0],keyPow[0],btmp1);
		XOR(a[1],keyPow[1],btmp2);
		schoolbook(btmp1,btmp2,multpart1,xx);
		XOR(a[2],multpart1[0],multpart1[0]);
		reductionbymult(multpart1,xx,btmp1);
		XOR(a[3],keyPow[2],btmp2);
		XOR(a[4],keyPow[0],btmp9);
		XOR(a[5],keyPow[1],btmp10);
		schoolbook(btmp1,btmp2,multpart1,xx);
		schoolbook(btmp9,btmp10,multpart2,xx1);
		XOR(a[6],multpart2[0],multpart2[0]);
		XORNEW(multpart1,multpart2,multpart1);
		XOR(xx,xx1,xx);
		reductionbymult(multpart1,xx,btmp1);
		XOR(a[7],keyPow[3],btmp2);
		XOR(a[8],keyPow[0],btmp9);
		XOR(a[9],keyPow[1],btmp10);
		schoolbook(btmp1,btmp2,multpart1,xx);
		schoolbook(btmp9,btmp10,multpart2,xx1);
		XOR(a[10],multpart2[0],multpart2[0]);
		reductionbymult(multpart2,xx1,btmp9);
		XOR(a[11],keyPow[2],btmp10);
		schoolbook(btmp9,btmp10,multpart2,xx1);
		XOR(a[12],multpart2[0],multpart2[0]);
		XORNEW(multpart1,multpart2,dig);
		XOR(xx,xx1,*mid);
		//reductionbymult(multpart1,xx,*dig);
	}

void partialfourteen(__m128i a[],__m128i dig[2], __m128i *mid)
	{
  		
		__m128i multpart1[2],multpart2[2],multpart3[2],POLY;
		__m128i btmp1,btmp2,btmp9,btmp10,xx,xx1,xx2;
		POLY = _mm_set_epi32(0x0,0x0,0x0,0x00000087);
		XOR(a[0],keyPow[0],btmp1);
		XOR(a[1],keyPow[1],btmp2);
		schoolbook(btmp1,btmp2,multpart1,xx);
		XOR(a[2],multpart1[0],multpart1[0]);
		reductionbymult(multpart1,xx,btmp1);
		XOR(a[3],keyPow[2],btmp2);
		XOR(a[4],keyPow[0],btmp9);
		XOR(a[5],keyPow[1],btmp10);
		schoolbook(btmp1,btmp2,multpart1,xx);
		schoolbook(btmp9,btmp10,multpart2,xx1);
		XOR(a[6],multpart2[0],multpart2[0]);
		XORNEW(multpart1,multpart2,multpart1);
		XOR(xx,xx1,xx);
		reductionbymult(multpart1,xx,btmp1);
		XOR(a[7],keyPow[3],btmp2);
		XOR(a[8],keyPow[0],btmp9);
		XOR(a[9],keyPow[1],btmp10);
		schoolbook(btmp1,btmp2,multpart1,xx);
		schoolbook(btmp9,btmp10,multpart2,xx1);
		XOR(a[10],multpart2[0],multpart2[0]);
		reductionbymult(multpart2,xx1,btmp9);
		XOR(a[11],keyPow[2],btmp10);
		schoolbook(btmp9,btmp10,multpart2,xx1);
		schoolbook(a[12],keyPow[0],multpart3,xx2);
		XOR(a[13],multpart3[0],multpart3[0]);
		XORNEW(multpart1,multpart2,multpart1);
		XOR(xx,xx1,xx);
		XORNEW(multpart1,multpart3,dig);
		XOR(xx,xx2,*mid);
		//reductionbymult(multpart1,xx,*dig);

	}

void partialfifteen(__m128i a[],__m128i dig[2], __m128i *mid)
	{
  		
		__m128i multpart1[2],multpart2[2],multpart3[2],POLY;
		__m128i btmp1,btmp2,btmp9,btmp10,xx,xx1,xx2;
		POLY = _mm_set_epi32(0x0,0x0,0x0,0x00000087);
		XOR(a[0],keyPow[0],btmp1);
		XOR(a[1],keyPow[1],btmp2);
		schoolbook(btmp1,btmp2,multpart1,xx);
		XOR(a[2],multpart1[0],multpart1[0]);
		reductionbymult(multpart1,xx,btmp1);
		XOR(a[3],keyPow[2],btmp2);
		XOR(a[4],keyPow[0],btmp9);
		XOR(a[5],keyPow[1],btmp10);
		schoolbook(btmp1,btmp2,multpart1,xx);
		schoolbook(btmp9,btmp10,multpart2,xx1);
		XOR(a[6],multpart2[0],multpart2[0]);
		XORNEW(multpart1,multpart2,multpart1);
		XOR(xx,xx1,xx);
		reductionbymult(multpart1,xx,btmp1);
		XOR(a[7],keyPow[3],btmp2);
		XOR(a[8],keyPow[0],btmp9);
		XOR(a[9],keyPow[1],btmp10);
		schoolbook(btmp1,btmp2,multpart1,xx);
		schoolbook(btmp9,btmp10,multpart2,xx1);
		XOR(a[10],multpart2[0],multpart2[0]);
		reductionbymult(multpart2,xx1,btmp9);
		XOR(a[11],keyPow[2],btmp10);
		XOR(a[12],keyPow[0],btmp1);
		XOR(a[13],keyPow[1],btmp2);
		schoolbook(btmp9,btmp10,multpart2,xx1);
		schoolbook(btmp1,btmp2,multpart3,xx2);
		XOR(a[14],multpart3[0],multpart3[0]);
		XORNEW(multpart1,multpart2,multpart1);
		XOR(xx,xx1,xx);
		XORNEW(multpart1,multpart3,dig);
		XOR(xx,xx2,*mid);
		//reductionbymult(multpart1,xx,*dig);

		
	}


		
void partialrest(__m128i a[],int noOfBlock,__m128i dig[2], __m128i  *mid1)
	{
		__m128i dig3[2],tmpdig,POLY, tmpMid, btmp1;
		POLY = _mm_set_epi32(0x0,0x0,0x0,0x00000087);
		
		partialfifteen(a,dig,mid1);
		reductionbymult(dig, *mid1, tmpdig);
		
		
		XOR(a[15],keyPow[4],btmp1);  		
		schoolbook(tmpdig,btmp1,dig,*mid1);
		//reductionbymult(tmpdig,xx,*dig);

		
		if(noOfBlock == 16)
			return;
		if(noOfBlock == 17)
		{
			XOR(a[16],dig[0],dig[0]);
			//reductionbymult(tmpdig,xx,*dig);
			return;

		}
		switch(noOfBlock)
		{
		
		case 18:
			partialtwo(a+16,dig3,&tmpMid);
			break;
		case 19:
			partialthree(a+16,dig3,&tmpMid);
			break;
		case 20:
			partialfour(a+16,dig3,&tmpMid);
			break;
		case 21:
			partialfive(a+16,dig3,&tmpMid);;
			break;
		case 22:
			partialsix(a+16,dig3,&tmpMid);
			break;
		case 23:
			partialseven(a+16,dig3,&tmpMid);
			break;
		case 24:
			partialeight(a+16,dig3,&tmpMid);
			break;
		case 25:
			partialnine(a+16,dig3,&tmpMid);
			break;
		case 26:
			partialten(a+16,dig3,&tmpMid);
			break;
		case 27:
			partialeleven(a+16,dig3,&tmpMid);
			break;
		case 28:
			partialtwelve(a+16,dig3,&tmpMid);
			break;
		case 29:
			partialthirteen(a+16,dig3,&tmpMid);
			break;
		case 30:
			partialfourteen(a+16,dig3,&tmpMid);
			break;
		}
		XORNEW(dig, dig3, dig);
		XOR(*mid1, tmpMid, *mid1);
	}
	


#endif
