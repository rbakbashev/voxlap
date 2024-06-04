.686
.XMM

USEZBUFFER EQU 1       ;To disable, put ; in front of line
LVSID EQU 10           ;log2(VSID) - used for mip-mapping index adjustment
VSID EQU (1 SHL LVSID) ;should match VSID in VOXLAP5.H (adjust LVSID, not this)
LOGPREC EQU (8+12)

EXTRN _gi : dword
EXTRN _gpixy : dword
EXTRN _gixy : dword     ;long[2]
EXTRN _gpz : dword      ;long[2]
EXTRN _gdz : dword      ;long[2]
EXTRN _gxmip : dword
EXTRN _gxmax : dword
EXTRN _gcsub : dword    ;long[4]
EXTRN _gylookup : dword ;long[256+4+128+4+...]
EXTRN _gmipnum : dword

EXTRN _sptr : dword

EXTRN _skyoff : dword   ;Memory offset to start of longitude line
EXTRN _skyxsiz : dword  ;Size of longitude line
EXTRN _skylat : dword   ;long[_skyxsiz] : latitude's unit dir. vector

CODE SEGMENT PUBLIC USE32 'CODE'
ASSUME cs:CODE,ds:CODE

PUBLIC _v5_asm_dep_unlock ;Data Execution Prevention unlock (works under XP2 SP2)
_v5_asm_dep_unlock:
	EXTRN __imp__VirtualProtect@16:NEAR
	sub esp, 4
	push esp
	push 40h ;PAGE_EXECUTE_READWRITE
	push offset _dep_protect_end - offset _v5_asm_dep_unlock
	push offset _v5_asm_dep_unlock
	call dword ptr __imp__VirtualProtect@16
	add esp, 4
	ret

PUBLIC _cfasm, _skycast
ALIGN 16
_cfasm db 256*32 dup(0)
w8bmask0 dq 000ff00ff00ff00ffh
w8bmask1 dq 000f000f000f000f0h
w8bmask2 dq 000e000e000e000e0h
;gyadd dq ((-1) SHL (LOGPREC-16))
mmask dq 0ffff0000ffff0000h
_skycast dq 0
gylookoff dd 0
ngxmax dd 0
ce dd 0
espbak dd 0

gylut  dd _gylookup
		 dd _gylookup+(4*1+256)*4
		 dd _gylookup+(4*2+384)*4
		 dd _gylookup+(4*3+448)*4
		 dd _gylookup+(4*4+480)*4
		 dd _gylookup+(4*5+496)*4
		 dd _gylookup+(4*6+504)*4
		 dd _gylookup+(4*7+508)*4
		 dd _gylookup+(4*8+510)*4

gxmipk dd ((1 SHL (LVSID-0))-1)*2
		 dd ((1 SHL (LVSID-1))-1)*2
		 dd ((1 SHL (LVSID-2))-1)*2
		 dd ((1 SHL (LVSID-3))-1)*2
		 dd ((1 SHL (LVSID-4))-1)*2
		 dd ((1 SHL (LVSID-5))-1)*2
		 dd ((1 SHL (LVSID-6))-1)*2
		 dd ((1 SHL (LVSID-7))-1)*2
		 dd ((1 SHL (LVSID-8))-1)*2

gymipk dd ((1 SHL (LVSID-0))-1) SHL (LVSID+2)
		 dd ((1 SHL (LVSID-1))-1) SHL (LVSID+1)
		 dd ((1 SHL (LVSID-2))-1) SHL (LVSID  )
		 dd ((1 SHL (LVSID-3))-1) SHL (LVSID-1)
		 dd ((1 SHL (LVSID-4))-1) SHL (LVSID-2)
		 dd ((1 SHL (LVSID-5))-1) SHL (LVSID-3)
		 dd ((1 SHL (LVSID-6))-1) SHL (LVSID-4)
		 dd ((1 SHL (LVSID-7))-1) SHL (LVSID-5)
		 dd ((1 SHL (LVSID-8))-1) SHL (LVSID-6)

gamipk dd _sptr
		 dd _sptr+(VSID*VSID)*4
		 dd _sptr+(VSID*VSID + (VSID*VSID SHR 2))*4
		 dd _sptr+(VSID*VSID + (VSID*VSID SHR 2) + (VSID*VSID SHR 4))*4
		 dd _sptr+(VSID*VSID + (VSID*VSID SHR 2) + (VSID*VSID SHR 4) + (VSID*VSID SHR 6))*4
		 dd _sptr+(VSID*VSID + (VSID*VSID SHR 2) + (VSID*VSID SHR 4) + (VSID*VSID SHR 6) + (VSID*VSID SHR 8))*4
		 dd _sptr+(VSID*VSID + (VSID*VSID SHR 2) + (VSID*VSID SHR 4) + (VSID*VSID SHR 6) + (VSID*VSID SHR 8) + (VSID*VSID SHR 10))*4
		 dd _sptr+(VSID*VSID + (VSID*VSID SHR 2) + (VSID*VSID SHR 4) + (VSID*VSID SHR 6) + (VSID*VSID SHR 8) + (VSID*VSID SHR 10) + (VSID*VSID SHR 12))*4
		 dd _sptr+(VSID*VSID + (VSID*VSID SHR 2) + (VSID*VSID SHR 4) + (VSID*VSID SHR 6) + (VSID*VSID SHR 8) + (VSID*VSID SHR 10) + (VSID*VSID SHR 12) + (VSID*VSID SHR 14))*4

gmipcnt db 0
ALIGN 16

	;THE INNER LOOP:
	;#ifdef CPU <= PENTIUM II
	;
	;   movd mm3, _gylookup[ecx*4] ;mm3: [ 0   0   0  -gy]
	;   por mm3, mm6               ;mm3: [ogx  0   gx -gy]
	;      or:
	;   paddd mm3, gyadd           ;where: gyadd: dq ((-1) SHL (LOGPREC-16))
	;
	;   ...
	;
	;   movq mm7, mm0           ;mm7: [cy0.... cx0....]
	;   psrad mm7, 16           ;mm7: [----cy0 ----cx0]
	;   packssdw mm7, mm7       ;mm7: [cy0 cx0 cy0 cx0]
	;   pmaddwd mm7, mm3           ;mm7: [ 0   0  -decide]
	;   movd eax, mm7
	;   test eax, eax
	;   j ?
	;   ...
	;   paddd mm0, _gi
	;
	;#else
	;      ;Do this only when gx/ogx changes
	;   movd mm3, ogx                   ;mm3: [ 0   0  ogx  0 ]
	;      or:
	;   pshufw mm3, mm3, 0e8h            ;mm3: [ gx ogx ogx  0 ]
	;
	;      ;Do this only when ecx/edx changes
	;   pinsrw mm3, _gylookup[ecx*2], 0 ;mm3: [ 0   0  ogx -gy]
	;      or:
	;   paddd mm3, gyadd                ;where: gyadd: dq (1 SHL LOGPREC)
	;
	;   ...
	;
	;   pshufw mm7, mm0, 0ddh      ;mm7: [cy0 cx0 cy0 cx0]
	;   pmaddwd mm7, mm3           ;mm7: [ 0   0  -decide]
	;   movd eax, mm7
	;   test eax, eax
	;   j ?
	;   ...
	;   paddd mm0, _gi
	;
	;#endif


	;   Register allocation:
	;eax: [.temp1.]     mm0: [cy0.... cx0....]
	;ebx: [.temp2.]     mm1: [cy1.... cx1....]
	;ecx: [.....z0]     mm2: [    temp!!!    ]   //gi[1].. gi[0]..]
	;edx: [.....z1]     mm3: [     temp      ]
	;esi: [..ixy..]     mm4: [??????? csub...]
	;edi: [..v[]..]     mm5: [??????? coltemp]
	;ebp: [...bakj]     mm6: [gx. 0.. ogx 0..]
	;esp: [..c->..]     mm7: [     temp      ]
PUBLIC _grouscanasm ;Visual C entry point (passes parameters by stack)
_grouscanasm:
	mov eax, [esp+4]
	push ebx   ;Visual C's _cdecl requires EBX,ESI,EDI,EBP to be preserved
	push esi
	push edi
	push ebp
	mov dword ptr espbak, esp

	mov edi, eax

		;cfasm:   0-2047  (extra memory for stack)
		;      2048-4095  c and ce always sit in this range ((esp = c) <= ce)
		;      4096-6143  This is where memory for cfasm is actually stored!
		;      6144-8191  (memory never used - this seems unnecessary?)
	mov esp, offset _cfasm[2048]
	mov eax, offset _cfasm[4096]
	mov ecx, [eax+8]
	mov edx, [eax+12]
	movq mm0, [eax+16]
	movq mm1, [eax+24]
	mov dword ptr ce, esp

	mov gylookoff, offset _gylookup
	mov byte ptr gmipcnt, 0

	mov ebp, _gxmax
	cmp byte ptr _gmipnum, 1
	jle short skipngxmax0
	cmp ebp, _gxmip
	jle short skipngxmax0
	mov ebp, _gxmip
skipngxmax0:
	mov ngxmax, ebp

	mov ebp, _gpz[4]
	sub ebp, _gpz[0]
	shr ebp, 31
	movd mm6, _gpz[ebp*4]        ;update gx in mm6
	pand mm6, qword ptr mmask
	mov eax, _gdz[ebp*4]
	add _gpz[ebp*4], eax

	mov esi, _gpixy
	cmp edi, [esi]
	je drawflor
	jmp drawceil

drawfwall:
	movzx eax, byte ptr [edi+1]
	cmp eax, edx
	jge drawcwall
	mov ebx, [esp+4+2048]
loop0:
	neg eax
	add eax, edx
	dec edx
	punpcklbw mm5, [edi+eax*4]
	mov eax, gylookoff
	movd mm3, [eax+edx*4] ;mm3: [ 0   0   0  -gy]
	psubusb mm5, mm4
	pshufw mm2, mm5, 0ffh
	pmulhuw mm5, mm2
	psrlw mm5, 7
	packuswb mm5, mm5
ifdef USEZBUFFER
	punpckldq mm5, mm6         ;Stuff ogx into hi part of color for Z buffer
endif
	por mm3, mm6               ;mm3: [ gx  0  ogx -gy]
loop1: ;if (dmulrethigh(gylookup[edx*4],c->cx1,c->cy1,ogx) >= 0) jmp endloop1
	pshufw mm7, mm1, 0ddh      ;mm7: [cy1 cx1 cy1 cx1]
	pmaddwd mm7, mm3           ;mm7: [ 0   0  -decide]
	movd eax, mm7
	test eax, eax              ;if (cy1*ogx ? gy*cx1)
	jle endloop1
	psubd mm1, qword ptr _gi
ifdef USEZBUFFER
	movntq [ebx], mm5
	sub ebx, 8
else
	movd [ebx], mm5
	sub ebx, 4
endif
	cmp ebx, [esp+2048]
	jnb loop1
	jmp predeletez
endloop1:
	movzx eax, byte ptr [edi+1]
	cmp eax, edx
	jne loop0
	mov [esp+4+2048], ebx

drawcwall:
	cmp edi, [esi]
	mov edx, eax
	je predrawflor

	movzx eax, byte ptr [edi+3]
	cmp eax, ecx
	jle predrawceil
	mov ebx, [esp+2048]
loop2:
	neg eax
	add eax, ecx
	inc ecx
	punpcklbw mm5, [edi+eax*4]
	mov eax, gylookoff
	movd mm3, [eax+ecx*4] ;mm3: [ 0   0   0  -gy]
	psubusb mm5, mm4
	pshufw mm2, mm5, 0ffh
	pmulhuw mm5, mm2
	psrlw mm5, 7
	packuswb mm5, mm5
ifdef USEZBUFFER
	punpckldq mm5, mm6         ;Stuff ogx into hi part of color for Z buffer
endif
	por mm3, mm6               ;mm3: [ gx  0  ogx -gy]
loop3: ;if (dmulrethigh(gylookup[ecx*4],c->cx0,c->cy0,ogx) < 0) jmp endloop3
	pshufw mm7, mm0, 0ddh      ;mm7: [cy0 cx0 cy0 cx0]
	pmaddwd mm7, mm3           ;mm7: [ 0   0  -decide]
	movd eax, mm7
	test eax, eax              ;if (cy0*ogx ? gy*cx0)
	jg endloop3
	paddd mm0, qword ptr _gi
ifdef USEZBUFFER
	movntq [ebx], mm5
	add ebx, 8
else
	movd [ebx], mm5
	add ebx, 4
endif
	cmp ebx, [esp+4+2048]
	jna loop3
	jmp predeletez
endloop3:
	movzx eax, byte ptr [edi+3]
	cmp eax, ecx
	jne loop2
	mov [esp+2048], ebx

predrawceil:
	mov ecx, eax
	pshufw mm6, mm6, 04eh       ;swap hi & lo of mm6
drawceil: ;if (dmulrethigh(gylookup[ecx*4],c->cx0,c->cy0,gx) < 0) jmp drawflor
	mov eax, gylookoff
	movd mm3, [eax+ecx*4] ;mm3: [ 0   0   0  -gy]
	por mm3, mm6               ;mm3: [ogx  0   gx -gy]
drawceilloop:
	pshufw mm7, mm0, 0ddh      ;mm7: [cy0 cx0 cy0 cx0]
	pmaddwd mm7, mm3           ;mm7: [ 0   0  -decide]
	movd eax, mm7
	test eax, eax              ;if (cy0*gx ? gy*cx0)
	jg drawflor
	paddd mm0, qword ptr _gi
	mov eax, [esp+2048]

	punpcklbw mm5, [edi-4]
	psubusb mm5, qword ptr _gcsub[16]
	pshufw mm2, mm5, 0ffh
	pmulhuw mm5, mm2
	psrlw mm5, 7
	packuswb mm5, mm5
ifdef USEZBUFFER
	punpckldq mm5, mm6         ;Stuff gx into hi part of color for Z buffer
	movntq [eax], mm5
	add eax, 8
else
	movd [eax], mm5
	add eax, 4
endif
	mov [esp+2048], eax
	cmp eax, [esp+4+2048]
	jna drawceilloop
	jmp deletez

predrawflor:
	pshufw mm6, mm6, 04eh       ;swap hi & lo of mm6
drawflor: ;if (dmulrethigh(gylookup[edx*4],c->cx1,c->cy1,gx) >= 0) jmp enddrawflor
	mov eax, gylookoff
	movd mm3, [eax+edx*4] ;mm3: [ 0   0   0  -gy]
	por mm3, mm6               ;mm3: [ogx  0   gx -gy]
drawflorloop:
	pshufw mm7, mm1, 0ddh      ;mm7: [cy1 cx1 cy1 cx1]
	pmaddwd mm7, mm3           ;mm7: [ 0   0  -decide]
	movd eax, mm7
	test eax, eax              ;if (cy1*gx ? gy*cx1)
	jle enddrawflor
	psubd mm1, qword ptr _gi
	mov eax, [esp+4+2048]

	punpcklbw mm5, [edi+4]
	psubusb mm5, qword ptr _gcsub[24]
	pshufw mm2, mm5, 0ffh
	pmulhuw mm5, mm2
	psrlw mm5, 7
	packuswb mm5, mm5
ifdef USEZBUFFER
	punpckldq mm5, mm6         ;Stuff gx into hi part of color for Z buffer
	movntq [eax], mm5
	sub eax, 8
else
	movd [eax], mm5   ;(Used to page fault here)
	sub eax, 4
endif
	mov [esp+4+2048], eax
	cmp eax, [esp+2048]
	jnb drawflorloop
	jmp deletez

enddrawflor:
	mov ebx, esp
afterdelete:
	sub esp, 32
	cmp esp, offset _cfasm[2048]
	jae skipixy

	movq mm4, qword ptr _gcsub[ebp*8]
	add esi, _gixy[ebp*4]
	mov ebp, _gpz[4]
	mov edi, [esi]
	sub ebp, _gpz[0]
	shr ebp, 31
	mov eax, _gpz[ebp*4]
	movd mm7, eax
	punpckldq mm6, mm7
	pand mm6, qword ptr mmask
	cmp eax, ngxmax
	ja remiporend
	add eax, _gdz[ebp*4]
	mov _gpz[ebp*4], eax
	mov esp, ce
	jmp skipixy2

skipixy:
	pshufw mm6, mm6, 04eh       ;swap hi & lo of mm6
skipixy2:
	cmp ebx, esp
	je skipixy3
	add ebx, 2048
	mov [ebx+8], ecx
	mov [ebx+12], edx
	movq [ebx+16], mm0
	movq [ebx+24], mm1
	lea ebx, [esp+2048]
	mov ecx, [ebx+8]
	mov edx, [ebx+12]
	movq mm0, [ebx+16]
	movq mm1, [ebx+24]
skipixy3:

		;Find highest intersecting vbuf slab
	cmp byte ptr [edi], 0
	je drawfwall
	mov ebx, gylookoff
	jmp intoslabloop
findslabloop:
	lea edi, [edi+eax*4]
	cmp byte ptr [edi], 0
	je drawfwall
intoslabloop:
	movzx eax, byte ptr [edi+2]
		;if (dmulrethigh(gylookup[[edi+2]*4+4],c->cx0,c->cy0,ogx) >= 0)
		;   jmp findslabloopbreak
	movd mm3, [ebx+eax*4+4]    ;mm3: [ 0   0   0  -gy]
	por mm3, mm6               ;mm3: [ gx  0  ogx -gy]
	pshufw mm7, mm0, 0ddh      ;mm7: [cy0 cx0 cy0 cx0]
	pmaddwd mm7, mm3           ;mm7: [ 0   0  -decide]
	movd eax, mm7
	test eax, eax              ;if (cy0*ogx ? ?y*cx0)

	movzx eax, byte ptr [edi]
	jg findslabloop

		;If next slab ALSO intersects, split _cfasm!
		;if (dmulrethigh(v[v[0]*4+3],c->cx1,c->cy1,ogx) >= 0) jmp drawfwall
	movzx eax, byte ptr [3+edi+eax*4]
	movd mm3, [ebx+eax*4]      ;mm3: [ 0   0   0  -gy]
	por mm3, mm6               ;mm3: [ gx  0  ogx -gy]
	pshufw mm7, mm1, 0ddh      ;mm7: [cy1 cx1 cy1 cx1]
	pmaddwd mm7, mm3           ;mm7: [ 0   0  -decide]
	movd eax, mm7
	test eax, eax              ;if (cy1*ogx ? ?y*cx1)
	jle drawfwall


		;Make sure everything is in memory at this point
	lea eax, [esp+2048]
	mov [eax+8], ecx
	mov [eax+12], edx
	movq [eax+16], mm0
	movq [eax+24], mm1

		;(ecx and edx are free registers at this point)

	mov edx, [eax+4]             ;col = (long)c->i1;
	movzx eax, byte ptr [edi+2]  ;dax = c->cx1; day = c->cy1;
	movd mm3, [ebx+eax*4+4]      ;mm3: [ 0   0   0  -gy]
	por mm3, mm6                 ;mm3: [ gx  0  ogx -gy]

		;WARNING: NEW CODE!!!!!!!
prebegsearchi16:
	movq mm7, qword ptr _gi
	pslld mm7, 4
	movq mm5, mm1
	psubd mm5, mm7             ;mm7: [day.... dax....]
	pshufw mm7, mm5, 0ddh      ;mm7: [day dax day dax]
	pmaddwd mm7, mm3           ;mm7: [ 0   0  -decide]
	movd eax, mm7
	test eax, eax              ;if (day*ogx ? gy*dax)
	jle begsearchi
	movq mm1, mm5
ifdef USEZBUFFER
	sub edx, 16 SHL 3
else
	sub edx, 16 SHL 2
endif
	jmp prebegsearchi16

	jmp begsearchi
		;while (dmulrethigh(gylookup[v[2]+1],dax,day,ogx) < 0)
prebegsearchi:
ifdef USEZBUFFER
	sub edx, 4 SHL 1             ;col -= 8;
else
	sub edx, 4                   ;col -= 4;
endif
	psubd mm1, qword ptr _gi     ;dax -= gi[0]; day -= gi[1];
begsearchi:
	pshufw mm7, mm1, 0ddh      ;mm7: [day dax day dax]
	pmaddwd mm7, mm3           ;mm7: [ 0   0  -decide]
	movd eax, mm7
	test eax, eax              ;if (day*ogx ? gy*dax)
	jg prebegsearchi

	mov eax, ce            ;ce++;
	add eax, 32

	cmp eax, offset _cfasm[4096] ;VERY BAD!!! - Interrupt would overwrite data!
	ja retsub                    ;Just in case, return early to prevent lockup.

	mov dword ptr ce, eax
	cmp eax, esp           ;for(c2=ce;c2>c;c2--)   //(c2 = eax)
	jbe skipinsertloop
beginsertloop:
	movq mm5, [eax-32+24+2048]  ;c2[0] = c2[-1];
	movq mm7, [eax-32+16+2048]
	movq [eax+24+2048], mm5
	movq [eax+16+2048], mm7
	movq mm5, [eax-32+8+2048]
	movq mm7, [eax-32+0+2048]
	movq [eax+8+2048], mm5
	movq [eax+0+2048], mm7
	sub eax, 32
	cmp eax, esp
	ja beginsertloop
skipinsertloop:

	movzx eax, byte ptr [edi]
	movq mm7, mm1              ;c[1].cx1 = dax; c[1].cy1 = day;
	paddd mm7, qword ptr _gi
	movzx eax, byte ptr [3+edi+eax*4]
	mov [esp+32+4+2048], edx        ;c[1].i1 = (long *)col;
ifdef USEZBUFFER
	add edx, 8                      ;c[0].i0 = (long *)(col+(4<<1));
else
	add edx, 4                      ;c[0].i0 = (long *)(col+4);
endif
	mov [esp+2048], edx
	mov edx, eax               ;c[1].z1 = c[0].z0 = v[(v[0]<<2)+3];
	mov [esp+8+2048], eax
	movq [esp+16+2048], mm7         ;c[0].cx0 = dax+gi[0]; c[0].cy0 = day+gi[1];
	add esp, 32                ;c++;
	jmp drawfwall

remiporend:
	mov al, gmipcnt
	inc al
	cmp al, byte ptr _gmipnum
	jge startsky
	mov gmipcnt, al

	sub esi, offset _sptr

	mov eax, esi
	shl eax, 29
	xor eax, _gixy[0]
	mov eax, _gdz[0]
	js short skipbladd0
	add _gpz[0], eax
skipbladd0:
	add eax, eax
	jno short skipremip0
	mov _gpz[0], 7fffffffh
	xor eax, eax
skipremip0:
	mov _gdz[0], eax

	mov [ebx+8+2048], ecx ;this is the official place to backup ecx

	mov eax, esi
	mov cl, byte ptr gmipcnt
	add cl, 31-1-2-LVSID
	shl eax, cl
	xor eax, _gixy[4]
	mov eax, _gdz[4]
	js short skipbladd1
	add _gpz[4], eax
skipbladd1:
	add eax, eax
	jno short skipremip1
	mov _gpz[4], 7fffffffh
	xor eax, eax
skipremip1:
	mov _gdz[4], eax

	shr esi, 2
	mov eax, esi
	movzx ecx, byte ptr gmipcnt
	and esi, gxmipk[ecx*4] ;mask for x (1:1024->512, etc...)
	and eax, gymipk[ecx*4] ;mask for y (1:1024->512, etc...)
	lea esi, [eax+esi*2]
	add esi, gamipk[ecx*4] ;add offset (1:sptr+1024*1024*4, etc...)

	movzx eax, byte ptr gmipcnt
	mov eax, gylut[eax*4]
	mov gylookoff, eax

	sar _gixy[4], 1

	mov eax, offset _cfasm[2048]
startremip0:
	shr dword ptr [eax+8+2048], 1
	inc dword ptr [eax+12+2048]
	shr dword ptr [eax+12+2048], 1
	add eax, 32
	cmp eax, ce
	jbe short startremip0

	mov eax, ngxmax
	cmp eax, _gxmax
	jae short startsky
	add eax, eax
	jo skipngxmax1 ;Make sure it doesn't overflow to negative!
	cmp eax, _gxmax
	jl short skipngxmax2
skipngxmax1:
	mov eax, _gxmax
skipngxmax2:
	mov ngxmax, eax

		;register fix-ups after here:
	mov ecx, [ebx+8+2048] ;this is the official place to restore ecx
	shr ecx, 1
	inc edx
	shr edx, 1

		;this makes grid transition clean
	mov ebp, _gpz[4]
	sub ebp, _gpz[0]
	shr ebp, 31
	mov eax, _gpz[ebp*4]
	add eax, _gdz[ebp*4]
	mov _gpz[ebp*4], eax
	mov edi, [esi]

	mov esp, ce
	jmp skipixy2

startsky:
	mov esp, offset _cfasm[2048]
	cmp esp, ce
	ja retsub
	mov esi, _skyoff
	test esi, esi
	jnz short prestartskyloop

;Sky not loaded, so fill with black ------------------------------------------
endprebegloop:
	movq mm5, _skycast
	mov eax, [esp+2048]
	mov ebx, [esp+4+2048]
	cmp eax, ebx
	ja short endnextloop
endbegloop:
ifdef USEZBUFFER
	movntq [eax], mm5
	add eax, 8
else
	movd [eax], mm5
	add eax, 4
endif
	cmp eax, ebx
	jbe short endbegloop
endnextloop:
	add esp, 32
	cmp esp, ce
	jbe short endprebegloop
	jmp short retsub

;Sky loaded: do texture mapping ----------------------------------------------

prestartskyloop:
	movq qword ptr [ebx+24+2048], mm1  ;Hack to make sure [cy0,cx0] is in memory for sky

	mov esi, _skyoff
	mov ecx, _skylat
	movd mm5, dword ptr _skycast[4]
	mov edi, _skyxsiz
startskyloop:
	mov eax, [esp+2048]
	mov ebx, [esp+4+2048]
	cmp eax, ebx
	ja short endskyslab
	movq mm1, [esp+24+2048]    ;mm1: [cy1.... cx1....]
preskysearch:
	psubd mm1, qword ptr _gi
skysearch:
	pshufw mm7, mm1, 0ddh      ;mm7: [cy1 cx1 cy1 cx1]
	movd mm3, [ecx+edi*4]      ;mm3: [       xvi -yvi]
	pmaddwd mm7, mm3           ;mm7: [ 0   0  -decide]
	movd edx, mm7
	sar edx, 31
	lea edi, [edi+edx]
	jnz short skysearch        ;if (cy1*xvi ? -yvi*cx1)

	movd mm6, dword ptr [esi+edi*4]
ifdef USEZBUFFER
	punpckldq mm6, mm5
	movntq [ebx], mm6
	sub ebx, 8
else
	movd [ebx], mm6
	sub ebx, 4
endif
	cmp eax, ebx
	jbe short preskysearch
endskyslab:
	add esp, 32
	cmp esp, ce
	jbe short startskyloop

;-----------------------------------------------------------------------

retsub:
	emms
	mov esp, dword ptr espbak
	pop ebp    ;Visual C's _cdecl requires EBX,ESI,EDI,EBP to be preserved
	pop edi
	pop esi
	pop ebx
	ret

predeletez:
	pshufw mm6, mm6, 04eh       ;swap hi & lo of mm6
deletez:
	mov ebx, ce
	sub ebx, 32
	cmp ebx, offset _cfasm[2048]
	jb retsub          ;nothing to fill - skip remiporend stuff!
	mov dword ptr ce, ebx

	add ebx, 32

	cmp esp, ebx       ;while (eax <= ce)
	jae afterdelete
	mov eax, esp
deleteloop:
	movq mm5, [eax+32+0+2048]
	movq mm7, [eax+32+8+2048]
	movq [eax+0+2048], mm5
	movq [eax+8+2048], mm7
	movq mm5, [eax+32+16+2048]
	movq mm7, [eax+32+24+2048]
	movq [eax+16+2048], mm5
	movq [eax+24+2048], mm7
	add eax, 32
	cmp eax, ebx
	jb deleteloop
	jmp afterdelete

;----------------------------------------------------------------------------

PUBLIC _opti4asm

ALIGN 16
_opti4asm dd 5*4 dup(0)        ;NOTE: this used by ?render

_dep_protect_end:
CODE ENDS
END
