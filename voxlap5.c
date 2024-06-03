#include <fcntl.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

#ifdef __GNUC__
#include <stdint.h>
#define INT_PTR intptr_t
#define UINT_PTR uintptr_t
#endif

#if !defined(_WIN32) && !defined(__DOS__)
#include <unistd.h>
#include <dirent.h>
typedef long long __int64;
static __inline int filelength (int h)
{
	struct stat st;
	if (fstat(h,&st) < 0) return(-1);
	return(st.st_size);
}
#define _fileno fileno
#else
#include <io.h>
#endif

#if !defined(max)
#define max(a,b) (((a) > (b)) ? (a) : (b))
#endif
#if !defined(min)
#define min(a,b) (((a) < (b)) ? (a) : (b))
#endif

#if defined(__GNUC__)
#define _inline inline
#endif

#define VOXLAP5

#define MAXXDIM 1024
#define MAXYDIM 768
#define PI 3.141592653589793
#define VSID 1024   //Maximum .VXL dimensions in both x & y direction
#define MAXZDIM 256 //Maximum .VXL dimensions in z direction (height)

#pragma pack(push,1)

typedef struct { long x, y, z; } lpoint3d;
typedef struct { float x, y, z; } point3d;
typedef struct { float x, y, z, z2; } point4d;
typedef struct { double x, y, z; } dpoint3d;

	//Sprite structures:
typedef struct { long col; unsigned short z; char vis, dir; } kv6voxtype;

typedef struct kv6data
{
	long leng, xsiz, ysiz, zsiz;
	float xpiv, ypiv, zpiv;
	unsigned long numvoxs;
	long namoff;
	kv6data *lowermip;
	kv6voxtype *vox;      //numvoxs*sizeof(kv6voxtype)
	unsigned long *xlen;  //xsiz*sizeof(long)
	unsigned short *ylen; //xsiz*ysiz*sizeof(short)
} kv6data;

#pragma pack(pop)

#define MAXFRM 1024 //MUST be even number for alignment!

	//Voxlap5 shared global variables:
#ifndef VOXLAP5
extern
#endif
struct
{
	//------------------------ DATA coming from VOXLAP5 ------------------------

		//Clipmove hit point info (use this after calling clipmove):
	double clipmaxcr; //clipmove always calls findmaxcr even with no movement
	dpoint3d cliphit[3];
	long cliphitnum;

		//Bounding box written by last set* VXL writing call
	long minx, miny, minz, maxx, maxy, maxz;

		//Falling voxels shared data:
	long flstnum;

		//Total count of solid voxels in .VXL map (included unexposed voxels)
	long globalmass;

	//------------------------ DATA provided to VOXLAP5 ------------------------

		//Opticast variables:
	long anginc, sideshademode, mipscandist, maxscandist, vxlmipuse, fogcol;

		//Drawsprite variables:
	long kv6col;

		//Map modification function data:
	long curcol, currad, curhei;
	float curpow;

		//Procedural texture function data:
	long (*colfunc)(lpoint3d *);
	long cen, amount, *pic, bpl, xsiz, ysiz, xoru, xorv, picmode;
	point3d fpico, fpicu, fpicv, fpicw;
	lpoint3d pico, picu, picv;
	float daf;

	long fallcheck;
} vx5;

	//Initialization functions:
extern long initvoxlap ();
extern void uninitvoxlap ();

	//File related functions:
extern long loadsxl (const char *, char **, char **, char **);
extern long loadvxl (const char *, dpoint3d *, dpoint3d *, dpoint3d *, dpoint3d *);

	//Screen related functions:
extern void voxsetframebuffer (long, long, long, long);
extern void setsideshades (char, char, char, char, char, char);
extern void setcamera (dpoint3d *, dpoint3d *, dpoint3d *, dpoint3d *, float, float, float);
extern void opticast ();

	//Physics helper functions:
extern void dorthorotate (double, double, double, dpoint3d *, dpoint3d *, dpoint3d *);
extern double findmaxcr (double, double, double, double);
extern void clipmove (dpoint3d *, dpoint3d *, double);
extern void estnorm (long, long, long, point3d *);

	//VXL reading functions (fast!):
extern long isvoxelsolid (long, long, long);
extern long anyvoxelsolid (long, long, long, long);
extern long anyvoxelempty (long, long, long, long);
extern long getfloorz (long, long, long);
extern long getcube (long, long, long);

	//VXL MISC functions:
extern void updatebbox (long, long, long, long, long, long, long);
extern void updatevxl ();
extern void genmipvxl (long, long, long, long);

	//Procedural texture functions:
extern long curcolfunc (lpoint3d *);
extern long floorcolfunc (lpoint3d *);
extern long jitcolfunc (lpoint3d *);
extern long manycolfunc (lpoint3d *);
extern long sphcolfunc (lpoint3d *);
extern long woodcolfunc (lpoint3d *);
extern long pngcolfunc (lpoint3d *);

	//ZIP functions:
extern int kzopen (const char *);
extern int kzread (void *, int);
extern int kzfilelength ();
extern int kztell ();
extern void kzclose ();

#include "sysmain.h"


//VOXLAP engine by Ken Silverman (http://advsys.net/ken)

#define USEZBUFFER 1

#define PREC (256*4096)
#define CMPPREC (256*4096)
#define FPREC (256*4096)
#define USEV5ASM 1
#define SCISDIST 1.0
#define GOLDRAT 0.3819660112501052 //Golden Ratio: 1 - 1/((sqrt(5)+1)/2)
#define ESTNORMRAD 2 //Specially optimized for 2: DON'T CHANGE unless testing!

#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#else
#include <stdarg.h>
#include <stdio.h>
#include <string.h>
#include <conio.h>
#include <dos.h>
#define MAX_PATH 260
#endif
#include <stdlib.h>

#define VOXSIZ VSID*VSID*128
#ifdef __cplusplus
extern "C" {
#endif
char *sptr[(VSID*VSID*4)/3];
#ifdef __cplusplus
}
#endif
static long *vbuf = 0, *vbit = 0, vbiti;
	//WARNING: loaddta uses last 2MB of vbuf; vbuf:[VOXSIZ>>2], vbit:[VOXSIZ>>7]
	//WARNING: loadpng uses last 4MB of vbuf; vbuf:[VOXSIZ>>2], vbit:[VOXSIZ>>7]

//                     ÚÄÄÄÄÄÄÄÄÂÄÄÄÄÄÄÄÄÂÄÄÄÄÄÄÄÄÂÄÄÄÄÄÄÄÄ¿
//        vbuf format: ³   0:   ³   1:   ³   2:   ³   3:   ³
//ÚÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÅÄÄÄÄÄÄÄÄÅÄÄÄÄÄÄÄÄÅÄÄÄÄÄÄÄÄÅÄÄÄÄÄÄÄÄ´
//³      First header: ³ nextptr³   z1   ³   z1c  ³  dummy ³
//³           Color 1: ³    b   ³    g   ³    r   ³ intens ³
//³           Color 2: ³    b   ³    g   ³    r   ³ intens ³
//³             ...    ³    b   ³    g   ³    r   ³ intens ³
//³           Color n: ³    b   ³    g   ³    r   ³ intens ³
//³ Additional header: ³ nextptr³   z1   ³   z1c  ³   z0   ³
//ÀÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÁÄÄÄÄÄÄÄÄÁÄÄÄÄÄÄÄÄÁÄÄÄÄÄÄÄÄÁÄÄÄÄÄÄÄÄÙ
//  nextptr: add this # <<2 to index to get to next header (0 if none)
//       z1: z floor (top of floor color list)
//      z1c: z bottom of floor color list MINUS 1! - needed to calculate
//             slab size with slng() and used as a separator for fcol/ccol
//       z0: z ceiling (bottom of ceiling color list)

	//Memory management variables:
#define MAXCSIZ 1028
char tbuf[MAXCSIZ];
long tbuf2[MAXZDIM*3];
long templongbuf[MAXZDIM];

extern long cputype; //bit25=1: SSE, bits30&31=1,1:3DNow!+

static char nullst = 0; //nullst always NULL string

#pragma pack(push,1)
	//Rendering variables:
#if (USEZBUFFER == 0)
typedef struct { long col; } castdat;
#else
typedef struct { long col, dist; } castdat;
#endif
typedef struct { castdat *i0, *i1; long z0, z1, cx0, cy0, cx1, cy1; } cftype;
typedef struct { unsigned short x, y; } uspoint2d;
typedef struct { long x, y; } lpoint2d;
typedef struct { float x, y; } point2d;
#pragma pack(pop)

#if USEV5ASM
#ifndef __cplusplus
	extern void *cfasm;
#else
	extern "C" void *cfasm;
#endif
	#define cf ((cftype *)&cfasm)
#else
	cftype cf[256];
#endif

	//Screen related variables:
static long xres, yres, bytesperline, frameplace, xres4;
long ylookup[MAXYDIM+1];

static lpoint3d glipos;
static point3d gipos, gistr, gihei, gifor;
static point3d gixs, giys, gizs, giadd;
static float gihx, gihy, gihz, gposxfrac[2], gposyfrac[2], grd;
static long gposz, giforzsgn, gstartz0, gstartz1, gixyi[2];
static char *gstartv;

long backtag, backedup = -1, bacx0, bacy0, bacx1, bacy1;
char *bacsptr[262144];

	//Flash variables
#define LOGFLASHVANG 9
static lpoint2d gfc[(1<<LOGFLASHVANG)*8];
static long gfclookup[8] = {4,7,2,5,0,3,6,1}, flashcnt = 0;
__int64 flashbrival;

	//Norm flash variables
#define GSIZ 512  //NOTE: GSIZ should be 1<<x, and must be <= 65536
static long bbuf[GSIZ][GSIZ>>5], p2c[32], p2m[32];      //bbuf: 2.0K
static uspoint2d ffx[((GSIZ>>1)+2)*(GSIZ>>1)], *ffxptr; // ffx:16.5K
static long xbsox = -17, xbsoy, xbsof;
static __int64 xbsbuf[25*5+1]; //need few bits before&after for protection

	//Look tables for expandbitstack256:
static long xbsceil[32], xbsflor[32];

	//float detection & falling code variables...
	//WARNING: VLSTSIZ,FSTKSIZ,FLCHKSIZ can all have bounds errors! :(
#define VLSTSIZ 65536 //Theoretically should be at least: VOXSIZ\8
#define LOGHASHEAD 12
#define FSTKSIZ 8192
typedef struct { long v, b; } vlstyp;
vlstyp vlst[VLSTSIZ];
long hhead[1<<LOGHASHEAD], vlstcnt = 0x7fffffff;
lpoint3d fstk[FSTKSIZ]; //Note .z is actually used as a pointer, not z!
#define FLCHKSIZ 4096
lpoint3d flchk[FLCHKSIZ]; long flchkcnt = 0;

	//Opticast global variables:
	//radar: 320x200 requires  419560*2 bytes (area * 6.56*2)
	//radar: 400x300 requires  751836*2 bytes (area * 6.27*2)
	//radar: 640x480 requires 1917568*2 bytes (area * 6.24*2)
#define SCPITCH 256
long *radar = 0, *radarmem = 0;
#if (USEZBUFFER == 1)
static long *zbuffermem = 0, zbuffersiz = 0;
#endif
static castdat *angstart[MAXXDIM*4], *gscanptr;
#define CMPRECIPSIZ MAXXDIM+32
static float cmprecip[CMPRECIPSIZ], wx0, wy0, wx1, wy1;
static long iwx0, iwy0, iwx1, iwy1;
static point3d gcorn[4];
		 point3d ginor[4]; //Should be static, but... necessary for stupid pingball hack :/
static long lastx[max(MAXYDIM,VSID)], uurendmem[MAXXDIM*2+8], *uurend;

	//Parallaxing sky variables:
static long skypic = 0, nskypic = 0, skybpl, skyysiz, skycurlng, skycurdir;
static float skylngmul;
static point2d *skylng = 0;

#ifdef __cplusplus
extern "C" {
#endif

	//Parallaxing sky variables (accessed by assembly code)
long skyoff = 0, skyxsiz, *skylat = 0;

__int64 gi, gcsub[8] =
{
	0xff00ff00ff00ff,0xff00ff00ff00ff,0xff00ff00ff00ff,0xff00ff00ff00ff,
	0xff00ff00ff00ff,0xff00ff00ff00ff,0xff00ff00ff00ff,0xff00ff00ff00ff
};
long gylookup[512+36], gmipnum = 0; //256+4+128+4+64+4+...
long gpz[2], gdz[2], gxmip, gxmax, gixy[2], gpixy;
static long gmaxscandist;

//long reax, rebx, recx, redx, resi, redi, rebp, resp, remm[16];
void v5_asm_dep_unlock();
void grouscanasm (long);
#if (USEZBUFFER == 1)
long zbufoff;
#endif
#ifdef __cplusplus
}
#endif
#define gi0 (((long *)&gi)[0])
#define gi1 (((long *)&gi)[1])

#ifdef _MSC_VER

#pragma warning(disable:4799) //I know how to use EMMS

static _inline void fcossin (float a, float *c, float *s)
{
	_asm
	{
		fld a
		fsincos
		mov eax, c
		fstp dword ptr [eax]
		mov eax, s
		fstp dword ptr [eax]
	}
}

static _inline void dcossin (double a, double *c, double *s)
{
	_asm
	{
		fld a
		fsincos
		mov eax, c
		fstp qword ptr [eax]
		mov eax, s
		fstp qword ptr [eax]
	}
}

static _inline void ftol (float f, long *a)
{
	_asm
	{
		mov eax, a
		fld f
		fistp dword ptr [eax]
	}
}

static _inline void dtol (double d, long *a)
{
	_asm
	{
		mov eax, a
		fld qword ptr d
		fistp dword ptr [eax]
	}
}

	//WARNING: This ASM code requires >= PPRO
static _inline double dbound (double d, double dmin, double dmax)
{
	_asm
	{
		fld dmin
		fld d
		fucomi st, st(1)   ;if (d < dmin)
		fcmovb st, st(1)   ;    d = dmin;
		fld dmax
		fxch st(1)
		fucomi st, st(1)   ;if (d > dmax)
		fcmovnb st, st(1)  ;    d = dmax;
		fstp d
		fucompp
	}
	return(d);
}

static _inline long mulshr16 (long a, long d)
{
	_asm
	{
		mov eax, a
		mov edx, d
		imul edx
		shrd eax, edx, 16
	}
}

static _inline __int64 mul64 (long a, long d)
{
	_asm
	{
		mov eax, a
		imul d
	}
}

static _inline long shldiv16 (long a, long b)
{
	_asm
	{
		mov eax, a
		mov edx, eax
		shl eax, 16
		sar edx, 16
		idiv b
	}
}

static _inline long isshldiv16safe (long a, long b)
{
	_asm
	{
		mov edx, a
		test edx, edx
		js short skipneg0
		neg edx
skipneg0:
		sar edx, 14

		mov eax, b
		test eax, eax
		js short skipneg1
		neg eax
skipneg1:
			;abs((a<<16)/b) < (1<<30) ;1 extra for good luck!
			;-abs(a)>>14 > -abs(b)    ;use -abs because safe for 0x80000000
			;eax-edx < 0
		sub eax, edx
		shr eax, 31
	}
}

static _inline long umulshr32 (long a, long d)
{
	_asm
	{
		mov eax, a
		mul d
		mov eax, edx
	}
}

static _inline long scale (long a, long d, long c)
{
	_asm
	{
		mov eax, a
		imul d
		idiv c
	}
}

static _inline long dmulrethigh (long b, long c, long a, long d)
{
	_asm
	{
		mov eax, a
		imul d
		mov ecx, eax
		push edx
		mov eax, b
		imul c
		sub eax, ecx
		pop ecx
		sbb edx, ecx
		mov eax, edx
	}
}

static _inline void copybuf (void *s, void *d, long c)
{
	_asm
	{
		push esi
		push edi
		mov esi, s
		mov edi, d
		mov ecx, c
		rep movsd
		pop edi
		pop esi
	}
}

static _inline void clearbuf (void *d, long c, long a)
{
	_asm
	{
		push edi
		mov edi, d
		mov ecx, c
		mov eax, a
		rep stosd
		pop edi
	}
}

#else
#pragma message ("Compiler says it isn't Visual C.")
#endif

	//if (a < 0) return(0); else if (a > b) return(b); else return(a);
static _inline long lbound0 (long a, long b) //b MUST be >= 0
{
	if ((unsigned long)a <= b) return(a);
	return((~(a>>31))&b);
}

	//if (a < b) return(b); else if (a > c) return(c); else return(a);
static _inline long lbound (long a, long b, long c) //c MUST be >= b
{
	c -= b;
	if ((unsigned long)(a-b) <= c) return(a);
	return((((b-a)>>31)&c) + b);
}

#define LSINSIZ 8 //Must be >= 2!
static point2d usintab[(1<<LSINSIZ)+(1<<(LSINSIZ-2))];
static void ucossininit ()
{
	long i, j;
	double a, ai, s, si, m;

	j = 0; usintab[0].y = 0.0;
	i = (1<<LSINSIZ)-1;
	ai = PI*(-2)/((float)(1<<LSINSIZ)); a = ((float)(-i))*ai;
	ai *= .5; m = sin(ai)*2; s = sin(a); si = cos(a+ai)*m; m = -m*m;
	for(;i>=0;i--)
	{
		usintab[i].y = s; s += si; si += s*m; //MUCH faster than next line :)
		//usintab[i].y = sin(i*PI*2/((float)(1<<LSINSIZ)));
		usintab[i].x = (usintab[j].y-usintab[i].y)/((float)(1<<(32-LSINSIZ)));
		j = i;
	}
	for(i=(1<<(LSINSIZ-2))-1;i>=0;i--) usintab[i+(1<<LSINSIZ)] = usintab[i];
}

	//Calculates cos & sin of 32-bit unsigned long angle in ~15 clock cycles
	//  Accuracy is approximately +/-.0001
static _inline void ucossin (unsigned long a, float *cosin)
{
	float f = ((float)(a&((1<<(32-LSINSIZ))-1))); a >>= (32-LSINSIZ);
	cosin[0] = usintab[a+(1<<(LSINSIZ-2))].x*f+usintab[a+(1<<(LSINSIZ-2))].y;
	cosin[1] = usintab[a                 ].x*f+usintab[a                 ].y;
}

static long gkrand = 0;
long colorjit (long i, long jitamount)
{
	gkrand = (gkrand*27584621)+1;
	return((gkrand&jitamount)^i);
}

	//Note: ebx = 512 is no change
	//If PENTIUM III:1.Replace punpcklwd&punpckldq with: pshufw mm1, mm1, 0
	//               2.Use pmulhuw, shift by 8 & mul by 256
	//  :(  Can't mix with floating point
//#pragma aux colormul =
//   "movd mm0, eax"
//   "pxor mm1, mm1"
//   "punpcklbw mm0, mm1"
//   "psllw mm0, 7"
//   "movd mm1, ebx"
//   "punpcklwd mm1, mm1"
//   "punpckldq mm1, mm1"
//   "pmulhw mm0, mm1"
//   "packsswb mm0, mm0"
//   "movd eax, mm0"
//   parm [eax][ebx]
//   modify exact [eax]
//   value [eax]

long colormul (long i, long mulup8)
{
	long r, g, b;

	r = ((((i>>16)&255)*mulup8)>>8); if (r > 255) r = 255;
	g = ((((i>>8 )&255)*mulup8)>>8); if (g > 255) g = 255;
	b = ((((i    )&255)*mulup8)>>8); if (b > 255) b = 255;
	return((i&0xff000000)+(r<<16)+(g<<8)+b);
}

long curcolfunc (lpoint3d *p) { return(vx5.curcol); }

long floorcolfunc (lpoint3d *p)
{
	char *v;
	for(v=sptr[p->y*VSID+p->x];(p->z>v[2]) && (v[0]);v+=v[0]*4);
	return(*(long *)&v[4]);
}

long jitcolfunc (lpoint3d *p) { return(colorjit(vx5.curcol,vx5.amount)); }

static long manycolukup[64] =
{
	  0,  1,  2,  5, 10, 15, 21, 29, 37, 47, 57, 67, 79, 90,103,115,
	127,140,152,165,176,188,198,208,218,226,234,240,245,250,253,254,
	255,254,253,250,245,240,234,226,218,208,198,188,176,165,152,140,
	128,115,103, 90, 79, 67, 57, 47, 37, 29, 21, 15, 10,  5,  2,  1
};
long manycolfunc (lpoint3d *p)
{
	return((manycolukup[p->x&63]<<16)+(manycolukup[p->y&63]<<8)+manycolukup[p->z&63]+0x80000000);
}

long sphcolfunc (lpoint3d *p)
{
	long i;
	ftol(sin((p->x+p->y+p->z-vx5.cen)*vx5.daf)*-96,&i);
	return(((i+128)<<24)|(vx5.curcol&0xffffff));
}

#define WOODXSIZ 46
#define WOODYSIZ 24
#define WOODZSIZ 24
static float wx[256], wy[256], wz[256], vx[256], vy[256], vz[256];
long woodcolfunc (lpoint3d *p)
{
	float col, u, a, f, dx, dy, dz;
	long i, c, xof, yof, tx, ty, xoff;

	if (*(long *)&wx[0] == 0)
	{
		for(i=0;i<256;i++)
		{
			wx[i] = WOODXSIZ * ((float)rand()/32768.0f-.5f) * .5f;
			wy[i] = WOODXSIZ * ((float)rand()/32768.0f-.5f) * .5f;
			wz[i] = WOODXSIZ * ((float)rand()/32768.0f-.5f) * .5f;

				//UNIFORM spherical randomization (see spherand.c)
			dz = 1.0f-(float)rand()/32768.0f*.04f;
			a = (float)rand()/32768.0f*PI*2.0f; fcossin(a,&dx,&dy);
			f = sqrt(1.0f-dz*dz); dx *= f; dy *= f;
				//??z: rings,  ?z?: vertical,  z??: horizontal (nice)
			vx[i] = dz; vy[i] = fabs(dy); vz[i] = dx;
		}
	}

		//(tx&,ty&) = top-left corner of current panel
	ty = p->y - (p->y%WOODYSIZ);
	xoff = ((ty/WOODYSIZ)*(ty/WOODYSIZ)*51721 + (p->z/WOODZSIZ)*357) % WOODXSIZ;
	tx = ((p->x+xoff) - (p->x+xoff)%WOODXSIZ) - xoff;

	xof = p->x - (tx + (WOODXSIZ>>1));
	yof = p->y - (ty + (WOODYSIZ>>1));

	c = ((((tx*429 + 4695) ^ (ty*341 + 4355) ^ 13643) * 2797) & 255);
	dx = xof - wx[c];
	dy = yof - wy[c];
	dz = (p->z%WOODZSIZ) - wz[c];

		//u = distance to center of randomly oriented cylinder
	u = vx[c]*dx + vy[c]*dy + vz[c]*dz;
	u = sqrt(dx*dx + dy*dy + dz*dz - u*u);

		//ring randomness
	u += sin((float)xof*.12 + (float)yof*.15) * .5;
	u *= (sin(u)*.05 + 1);

		//Ring function: smooth saw-tooth wave
	col = sin(u*2)*24;
	col *= pow(1.f-vx[c],.3f);

		//Thin shaded borders
	if ((p->x-tx == 0) || (p->y-ty == 0)) col -= 5;
	if ((p->x-tx == WOODXSIZ-1) || (p->y-ty == WOODYSIZ-1)) col -= 3;

	//f = col+c*.12+72; i = ftolp3(&f);
	  ftol(col+c*.12f+72.0f,&i);

	return(colormul(vx5.curcol,i<<1));
}

long gxsizcache = 0, gysizcache = 0;
long pngcolfunc (lpoint3d *p)
{
	long x, y, z, u, v;
	float fx, fy, fz, rx, ry, rz;

	if (!vx5.pic) return(vx5.curcol);
	switch(vx5.picmode)
	{
		case 0:
			x = p->x-vx5.pico.x; y = p->y-vx5.pico.y; z = p->z-vx5.pico.z;
			u = (((x&vx5.picu.x) + (y&vx5.picu.y) + (z&vx5.picu.z))^vx5.xoru);
			v = (((x&vx5.picv.x) + (y&vx5.picv.y) + (z&vx5.picv.z))^vx5.xorv);
			break;
		case 1: case 2:
			fx = (float)p->x-vx5.fpico.x;
			fy = (float)p->y-vx5.fpico.y;
			fz = (float)p->z-vx5.fpico.z;
			rx = vx5.fpicu.x*fx + vx5.fpicu.y*fy + vx5.fpicu.z*fz;
			ry = vx5.fpicv.x*fx + vx5.fpicv.y*fy + vx5.fpicv.z*fz;
			rz = vx5.fpicw.x*fx + vx5.fpicw.y*fy + vx5.fpicw.z*fz;
			ftol(atan2(ry,rx)*vx5.xoru/(PI*2),&u);
			if (vx5.picmode == 1) ftol(rz,&v);
			else ftol((atan2(rz,sqrt(rx*rx+ry*ry))/PI+.5)*vx5.ysiz,&v);
			break;
		default: //case 3:
			fx = (float)p->x-vx5.fpico.x;
			fy = (float)p->y-vx5.fpico.y;
			fz = (float)p->z-vx5.fpico.z;
			ftol(vx5.fpicu.x*fx + vx5.fpicu.y*fy + vx5.fpicu.z*fz,&u);
			ftol(vx5.fpicv.x*fx + vx5.fpicv.y*fy + vx5.fpicv.z*fz,&v);
			break;
	}
	if ((unsigned long)(u-gxsizcache) >= (unsigned long)vx5.xsiz)
		if (u < 0) gxsizcache = u-(u+1)%vx5.xsiz-vx5.xsiz+1; else gxsizcache = u-(u%vx5.xsiz);
	if ((unsigned long)(v-gysizcache) >= (unsigned long)vx5.ysiz)
		if (v < 0) gysizcache = v-(v+1)%vx5.ysiz-vx5.ysiz+1; else gysizcache = v-(v%vx5.ysiz);
	return((vx5.pic[(v-gysizcache)*(vx5.bpl>>2)+(u-gxsizcache)]&0xffffff)|0x80000000);
}

	//Special case for SETSEC & SETCEI bumpmapping (vx5.picmode == 3)
	//no safety checks, returns alpha as signed char in range: (-128 to 127)
long hpngcolfunc (point3d *p)
{
	long u, v;
	float fx, fy, fz;

	fx = p->x-vx5.fpico.x;
	fy = p->y-vx5.fpico.y;
	fz = p->z-vx5.fpico.z;
	ftol(vx5.fpicu.x*fx + vx5.fpicu.y*fy + vx5.fpicu.z*fz,&u);
	ftol(vx5.fpicv.x*fx + vx5.fpicv.y*fy + vx5.fpicv.z*fz,&v);

	if ((unsigned long)(u-gxsizcache) >= (unsigned long)vx5.xsiz)
		if (u < 0) gxsizcache = u-(u+1)%vx5.xsiz-vx5.xsiz+1; else gxsizcache = u-(u%vx5.xsiz);
	if ((unsigned long)(v-gysizcache) >= (unsigned long)vx5.ysiz)
		if (v < 0) gysizcache = v-(v+1)%vx5.ysiz-vx5.ysiz+1; else gysizcache = v-(v%vx5.ysiz);
	return(vx5.pic[(v-gysizcache)*(vx5.bpl>>2)+(u-gxsizcache)]>>24);
}

static long slng (const char *s)
{
	const char *v;

	for(v=s;v[0];v+=v[0]*4);
	return((long)v-(long)s+(v[2]-v[1]+1)*4+4);
}

void voxdealloc (const char *v)
{
	long i, j;
	i = (((long)v-(long)vbuf)>>2); j = (slng(v)>>2)+i;
#if 0
	while (i < j) { vbit[i>>5] &= ~(1<<i); i++; }
#else
	if (!((j^i)&~31))
		vbit[i>>5] &= ~(p2m[j&31]^p2m[i&31]);
	else
	{
		vbit[i>>5] &=   p2m[i&31];  i >>= 5;
		vbit[j>>5] &= (~p2m[j&31]); j >>= 5;
		for(j--;j>i;j--) vbit[j] = 0;
	}
#endif
}

	//Note: danum MUST be a multiple of 4!
char *voxalloc (long danum)
{
	long i, badcnt, p0, p1, vend;

	badcnt = 0; danum >>= 2; vend = (VOXSIZ>>2)-danum;
	do
	{
		for(;vbiti<vend;vbiti+=danum)
		{
			if (vbit[vbiti>>5]&(1<<vbiti)) continue;
			for(p0=vbiti;(!(vbit[(p0-1)>>5]&(1<<(p0-1))));p0--);
			for(p1=p0+danum-1;p1>vbiti;p1--)
				if (vbit[p1>>5]&(1<<p1)) goto allocnothere;

			vbiti = p0+danum;
			for(i=p0;i<vbiti;i++) vbit[i>>5] |= (1<<i);
			return((char *)(&vbuf[p0]));
allocnothere:;
		}
		vbiti = 0; badcnt++;
	} while (badcnt < 2);
	evilquit("voxalloc: vbuf full"); return(0);
}

long isvoxelsolid (long x, long y, long z)
{
	char *v;

	if ((unsigned long)(x|y) >= VSID) return(0);
	v = sptr[y*VSID+x];
	while (1)
	{
		if (z < v[1]) return(0);
		if (!v[0]) return(1);
		v += v[0]*4;
		if (z < v[3]) return(1);
	}
}

	//Returns 1 if any voxels in range (x,y,z0) to (x,y,z1-1) are solid, else 0
long anyvoxelsolid (long x, long y, long z0, long z1)
{
	char *v;

		//         v1.....v3   v1.....v3    v1.......................>
		//                z0.........z1
	if ((unsigned long)(x|y) >= VSID) return(0);
	v = sptr[y*VSID+x];
	while (1)
	{
		if (z1 <= v[1]) return(0);
		if (!v[0]) return(1);
		v += v[0]*4;
		if (z0 < v[3]) return(1);
	}
}

	//Returns 1 if any voxels in range (x,y,z0) to (x,y,z1-1) are empty, else 0
long anyvoxelempty (long x, long y, long z0, long z1)
{
	char *v;

		//         v1.....v3   v1.....v3    v1.......................>
		//                z0.........z1
	if ((unsigned long)(x|y) >= VSID) return(1);
	v = sptr[y*VSID+x];
	while (1)
	{
		if (z0 < v[1]) return(1);
		if (!v[0]) return(0);
		v += v[0]*4;
		if (z1 <= v[3]) return(0);
	}
}

	//Returns z of first solid voxel under (x,y,z). Returns z if in solid.
long getfloorz (long x, long y, long z)
{
	char *v;

	if ((unsigned long)(x|y) >= VSID) return(z);
	v = sptr[y*VSID+x];
	while (1)
	{
		if (z <= v[1]) return(v[1]);
		if (!v[0]) break;
		v += v[0]*4;
		if (z < v[3]) break;
	}
	return(z);
}

	//Returns:
	//   0: air
	//   1: unexposed solid
	//else: address to color in vbuf (this can never be 0 or 1)
long getcube (long x, long y, long z)
{
	long ceilnum;
	char *v;

	if ((unsigned long)(x|y) >= VSID) return(0);
	v = sptr[y*VSID+x];
	while (1)
	{
		if (z <= v[2])
		{
			if (z < v[1]) return(0);
			return((long)&v[(z-v[1])*4+4]);
		}
		ceilnum = v[2]-v[1]-v[0]+2;

		if (!v[0]) return(1);
		v += v[0]*4;

		if (z < v[3])
		{
			if (z-v[3] < ceilnum) return(1);
			return((long)&v[(z-v[3])*4]);
		}
	}
}

	// Inputs: uind[MAXZDIM]: uncompressed 32-bit color buffer (-1: air)
	//         nind?[MAXZDIM]: neighbor buf:
	//            -2: unexposed solid
	//            -1: air
	//    0-16777215: exposed solid (color)
	//         px,py: parameters for setting unexposed voxel colors
	//Outputs: cbuf[MAXCSIZ]: compressed output buffer
	//Returns: n: length of compressed buffer (in bytes)
long compilestack (long *uind, long *n0, long *n1, long *n2, long *n3, char *cbuf, long px, long py)
{
	long oz, onext, n, cp2, cp1, cp0, rp1, rp0;
	lpoint3d p;

	p.x = px; p.y = py;

		//Do top slab (sky)
	oz = -1;
	p.z = -1; while (uind[p.z+1] == -1) p.z++;
	onext = 0;
	cbuf[1] = p.z+1;
	cbuf[2] = p.z+1;
	cbuf[3] = 0;  //Top z0 (filler, not used yet)
	n = 4;
	cp1 = 1; cp0 = 0;
	rp1 = -1; rp0 = -1;

	do
	{
			//cp2 = state at p.z-1 (0 = air, 1 = next2air, 2 = solid)
			//cp1 = state at p.z   (0 = air, 1 = next2air, 2 = solid)
			//cp0 = state at p.z+1 (0 = air, 1 = next2air, 2 = solid)
		cp2 = cp1; cp1 = cp0; cp0 = 2;
		if (p.z < MAXZDIM-2)  //Bottom must be solid!
		{
			if (uind[p.z+1] == -1)
				cp0 = 0;
			else if ((n0[p.z+1] == -1) || (n1[p.z+1] == -1) ||
						(n2[p.z+1] == -1) || (n3[p.z+1] == -1))
				cp0 = 1;
		}

			//Add slab
		if (cp1 != rp0)
		{
			if ((!cp1) && (rp0 > 0)) { oz = p.z; }
			else if ((rp0 < cp1) && (rp0 < rp1))
			{
				if (oz < 0) oz = p.z;
				cbuf[onext] = ((n-onext)>>2); onext = n;
				cbuf[n+1] = p.z;
				cbuf[n+2] = p.z-1;
				cbuf[n+3] = oz;
				n += 4; oz = -1;
			}
			rp1 = rp0; rp0 = cp1;
		}

			//Add color
		if ((cp1 == 1) || ((cp1 == 2) && ((!cp0) || (!cp2))))
		{
			if (cbuf[onext+2] == p.z-1) cbuf[onext+2] = p.z;
			if (uind[p.z] == -2) *(long *)&cbuf[n] = vx5.colfunc(&p);
								 else *(long *)&cbuf[n] = uind[p.z];
			n += 4;
		}

		p.z++;
	} while (p.z < MAXZDIM);
	cbuf[onext] = 0;
	return(n);
}

#ifdef _MSC_VER

static _inline void expandbit256 (void *s, void *d)
{
	_asm
	{
		push esi
		push edi
		mov esi, s
		mov edi, d
		mov ecx, 32   ;current bit index
		xor edx, edx  ;value of current 32-bit bits
		jmp short in2it
begit:lea esi, [esi+eax*4]
		movzx eax, byte ptr [esi+3]
		sub eax, ecx              ;xor mask [eax] for ceiling begins
		jl short xskpc
xdoc: mov [edi], edx
		add edi, 4
		mov edx, -1
		add ecx, 32
		sub eax, 32
		jge short xdoc
xskpc:and edx, xbsceil[eax*4+128] ;~(-1<<eax); xor mask [eax] for ceiling ends
in2it:movzx eax, byte ptr [esi+1]
		sub eax, ecx              ;xor mask [eax] for floor begins
		jl short xskpf
xdof: mov [edi], edx
		add edi, 4
		xor edx, edx
		add ecx, 32
		sub eax, 32
		jge short xdof
xskpf:or edx, xbsflor[eax*4+128] ;(-1<<eax); xor mask [eax] for floor ends
		movzx eax, byte ptr [esi]
		test eax, eax
		jnz short begit
		sub ecx, 256              ;finish writing buffer to [edi]
		jg short xskpe
xdoe: mov [edi], edx
		add edi, 4
		mov edx, -1
		add ecx, 32
		jle short xdoe
xskpe:pop edi
		pop esi
	}
}

#endif

void expandbitstack (long x, long y, __int64 *bind)
{
	if ((x|y)&(~(VSID-1))) { clearbuf((void *)bind,8,0L); return; }
	expandbit256(sptr[y*VSID+x],(void *)bind);
}

void expandstack (long x, long y, long *uind)
{
	long z, topz;
	char *v, *v2;

	if ((x|y)&(~(VSID-1))) { clearbuf((void *)uind,MAXZDIM,0); return; }

		//Expands compiled voxel info to 32-bit uind[?]
	v = sptr[y*VSID+x]; z = 0;
	while (1)
	{
		while (z < v[1]) { uind[z] = -1; z++; }
		while (z <= v[2]) { uind[z] = (*(long *)&v[(z-v[1])*4+4]); z++; }
		v2 = &v[(v[2]-v[1]+1)*4+4];

		if (!v[0]) break;
		v += v[0]*4;

		topz = v[3]+(((long)v2-(long)v)>>2);
		while (z < topz) { uind[z] = -2; z++; }
		while (z < v[3]) { uind[z] = *(long *)v2; z++; v2 += 4; }
	}
	while (z < MAXZDIM) { uind[z] = -2; z++; }
}

void gline (long leng, float x0, float y0, float x1, float y1)
{
	unsigned __int64 q;
	float f, f1, f2, vd0, vd1, vz0, vx1, vy1, vz1;
	long j;
	cftype *c;
#if (USEV5ASM == 0)
	long gx, ogx, gy, ixy, col, dax, day;
	cftype *c2, *ce;
	char *v;
#endif

	vd0 = x0*gistr.x + y0*gihei.x + gcorn[0].x;
	vd1 = x0*gistr.y + y0*gihei.y + gcorn[0].y;
	vz0 = x0*gistr.z + y0*gihei.z + gcorn[0].z;
	vx1 = x1*gistr.x + y1*gihei.x + gcorn[0].x;
	vy1 = x1*gistr.y + y1*gihei.y + gcorn[0].y;
	vz1 = x1*gistr.z + y1*gihei.z + gcorn[0].z;

	f = sqrt(vx1*vx1 + vy1*vy1);
	f1 = f / vx1;
	f2 = f / vy1;
	if (fabs(vx1) > fabs(vy1)) vd0 = vd0*f1; else vd0 = vd1*f2;
	if (*(long *)&vd0 < 0) vd0 = 0; //vd0 MUST NOT be negative: bad for asm
	vd1 = f;
	ftol(fabs(f1)*PREC,&gdz[0]);
	ftol(fabs(f2)*PREC,&gdz[1]);

	gixy[0] = (((*(signed long *)&vx1)>>31)<<3)+4; //=sgn(vx1)*4
	gixy[1] = gixyi[(*(unsigned long *)&vy1)>>31]; //=sgn(vy1)*4*VSID
	if (gdz[0] <= 0) { ftol(gposxfrac[(*(unsigned long *)&vx1)>>31]*fabs(f1)*PREC,&gpz[0]); if (gpz[0] <= 0) gpz[0] = 0x7fffffff; gdz[0] = 0x7fffffff-gpz[0]; } //Hack for divide overflow
	else ftol(gposxfrac[(*(unsigned long *)&vx1)>>31]*(float)gdz[0],&gpz[0]);
	if (gdz[1] <= 0) { ftol(gposyfrac[(*(unsigned long *)&vy1)>>31]*fabs(f2)*PREC,&gpz[1]); if (gpz[1] <= 0) gpz[1] = 0x7fffffff; gdz[1] = 0x7fffffff-gpz[1]; } //Hack for divide overflow
	else ftol(gposyfrac[(*(unsigned long *)&vy1)>>31]*(float)gdz[1],&gpz[1]);

	c = &cf[128];
	c->i0 = gscanptr; c->i1 = &gscanptr[leng];
	c->z0 = gstartz0; c->z1 = gstartz1;
	if (giforzsgn < 0)
	{
		ftol((vd1-vd0)*cmprecip[leng],&gi0); ftol(vd0*CMPPREC,&c->cx0);
		ftol((vz1-vz0)*cmprecip[leng],&gi1); ftol(vz0*CMPPREC,&c->cy0);
	}
	else
	{
		ftol((vd0-vd1)*cmprecip[leng],&gi0); ftol(vd1*CMPPREC,&c->cx0);
		ftol((vz0-vz1)*cmprecip[leng],&gi1); ftol(vz1*CMPPREC,&c->cy0);
	}
	c->cx1 = leng*gi0 + c->cx0;
	c->cy1 = leng*gi1 + c->cy0;

	gxmax = gmaxscandist;

		//Hack for early-out case when looking up towards sky
#if 0  //DOESN'T WORK WITH LOWER MIPS!
	if (c->cy1 < 0)
		if (gposz > 0)
		{
			if (dmulrethigh(-gposz,c->cx1,c->cy1,gxmax) >= 0)
			{
				j = scale(-gposz,c->cx1,c->cy1)+PREC; //+PREC for good luck
				if ((unsigned long)j < (unsigned long)gxmax) gxmax = j;
			}
		} else gxmax = 0;
#endif

		//Clip borders safely (MUST use integers!) - don't wrap around
	if (gixy[0] < 0) j = glipos.x; else j = VSID-1-glipos.x;
	q = mul64(gdz[0],j); q += (unsigned __int64)gpz[0];
	if (q < (unsigned __int64)gxmax)
	{
		gxmax = (long)q;
	}
	if (gixy[1] < 0) j = glipos.y; else j = VSID-1-glipos.y;
	q = mul64(gdz[1],j); q += (unsigned __int64)gpz[1];
	if (q < (unsigned __int64)gxmax)
	{
		gxmax = (long)q;
	}

	if (vx5.sideshademode)
	{
		gcsub[0] = gcsub[(((unsigned long)gixy[0])>>31)+4];
		gcsub[1] = gcsub[(((unsigned long)gixy[1])>>31)+6];
	}

#if USEV5ASM
	if (nskypic)
	{
		if (skycurlng < 0)
		{
			ftol((atan2(vy1,vx1)+PI)*skylngmul-.5,&skycurlng);
			if ((unsigned long)skycurlng >= skyysiz)
				skycurlng = ((skyysiz-1)&(j>>31));
		}
		else if (skycurdir < 0)
		{
			j = skycurlng+1; if (j >= skyysiz) j = 0;
			while (skylng[j].x*vy1 > skylng[j].y*vx1)
				{ skycurlng = j++; if (j >= skyysiz) j = 0; }
		}
		else
		{
			while (skylng[skycurlng].x*vy1 < skylng[skycurlng].y*vx1)
				{ skycurlng--; if (skycurlng < 0) skycurlng = skyysiz-1; }
		}
		skyoff = skycurlng*skybpl + nskypic;
	}

	//resp = 0;
	grouscanasm((long)gstartv);
	//if (resp)
	//{
	//   static char tempbuf[2048], tempbuf2[256];
	//   sprintf(tempbuf,"eax:%08x\tmm0:%08x%08x\nebx:%08x\tmm1:%08x%08x\necx:%08x\tmm2:%08x%08x\nedx:%08x\tmm3:%08x%08x\nesi:%08x\tmm4:%08x%08x\nedi:%08x\tmm5:%08x%08x\nebp:%08x\tmm6:%08x%08x\nesp:%08x\tmm7:%08x%08x\n",
	//      reax,remm[ 1],remm[ 0], rebx,remm[ 3],remm[ 2],
	//      recx,remm[ 5],remm[ 4], redx,remm[ 7],remm[ 6],
	//      resi,remm[ 9],remm[ 8], redi,remm[11],remm[10],
	//      rebp,remm[13],remm[12], resp,remm[15],remm[14]);
	//
	//   for(j=0;j<3;j++)
	//   {
	//      sprintf(tempbuf2,"%d i0:%d i1:%d z0:%ld z1:%ld cx0:%08x cy0:%08x cx1:%08x cy1:%08x\n",
	//         j,(long)cf[j].i0-(long)gscanptr,(long)cf[j].i1-(long)gscanptr,cf[j].z0,cf[j].z1,cf[j].cx0,cf[j].cy0,cf[j].cx1,cf[j].cy1);
	//      strcat(tempbuf,tempbuf2);
	//   }
	//   evilquit(tempbuf);
	//}
#else
//------------------------------------------------------------------------
	ce = c; v = gstartv;
	j = (((unsigned long)(gpz[1]-gpz[0]))>>31);
	gx = gpz[j];
	ixy = gpixy;
	if (v == (char *)*(long *)gpixy) goto drawflor; goto drawceil;

	while (1)
	{

drawfwall:;
		if (v[1] != c->z1)
		{
			if (v[1] > c->z1) c->z1 = v[1];
			else { do
			{
				c->z1--; col = *(long *)&v[(c->z1-v[1])*4+4];
				while (dmulrethigh(gylookup[c->z1],c->cx1,c->cy1,ogx) < 0)
				{
					c->i1->col = col; c->i1--; if (c->i0 > c->i1) goto deletez;
					c->cx1 -= gi0; c->cy1 -= gi1;
				}
			} while (v[1] != c->z1); }
		}

		if (v == (char *)*(long *)ixy) goto drawflor;

//drawcwall:;
		if (v[3] != c->z0)
		{
			if (v[3] < c->z0) c->z0 = v[3];
			else { do
			{
				c->z0++; col = *(long *)&v[(c->z0-v[3])*4-4];
				while (dmulrethigh(gylookup[c->z0],c->cx0,c->cy0,ogx) >= 0)
				{
					c->i0->col = col; c->i0++; if (c->i0 > c->i1) goto deletez;
					c->cx0 += gi0; c->cy0 += gi1;
				}
			} while (v[3] != c->z0); }
		}

drawceil:;
		while (dmulrethigh(gylookup[c->z0],c->cx0,c->cy0,gx) >= 0)
		{
			c->i0->col = (*(long *)&v[-4]); c->i0++; if (c->i0 > c->i1) goto deletez;
			c->cx0 += gi0; c->cy0 += gi1;
		}

drawflor:;
		while (dmulrethigh(gylookup[c->z1],c->cx1,c->cy1,gx) < 0)
		{
			c->i1->col = *(long *)&v[4]; c->i1--; if (c->i0 > c->i1) goto deletez;
			c->cx1 -= gi0; c->cy1 -= gi1;
		}

afterdelete:;
		c--;
		if (c < &cf[128])
		{
			ixy += gixy[j];
			gpz[j] += gdz[j];
			j = (((unsigned long)(gpz[1]-gpz[0]))>>31);
			ogx = gx; gx = gpz[j];

			if (gx > gxmax) break;
			v = (char *)*(long *)ixy; c = ce;
		}
			//Find highest intersecting vbuf slab
		while (1)
		{
			if (!v[0]) goto drawfwall;
			if (dmulrethigh(gylookup[v[2]+1],c->cx0,c->cy0,ogx) >= 0) break;
			v += v[0]*4;
		}
			//If next slab ALSO intersects, split cf!
		gy = gylookup[v[v[0]*4+3]];
		if (dmulrethigh(gy,c->cx1,c->cy1,ogx) < 0)
		{
			col = (long)c->i1; dax = c->cx1; day = c->cy1;
			while (dmulrethigh(gylookup[v[2]+1],dax,day,ogx) < 0)
				{ col -= sizeof(castdat); dax -= gi0; day -= gi1; }
			ce++; if (ce >= &cf[192]) return; //Give it max=64 entries like ASM
			for(c2=ce;c2>c;c2--) c2[0] = c2[-1];
			c[1].i1 = (castdat *)col; c->i0 = ((castdat *)col)+1;
			c[1].cx1 = dax; c->cx0 = dax+gi0;
			c[1].cy1 = day; c->cy0 = day+gi1;
			c[1].z1 = c->z0 = v[v[0]*4+3];
			c++;
		}
	}
//------------------------------------------------------------------------

	for(c=ce;c>=&cf[128];c--)
		while (c->i0 <= c->i1) { c->i0->col = 0; c->i0++; }
	return;

deletez:;
	ce--; if (ce < &cf[128]) return;
	for(c2=c;c2<=ce;c2++) c2[0] = c2[1];
	goto afterdelete;
#endif
}

#ifdef _MSC_VER

static _inline void mmxcoloradd (long *a)
{
	_asm
	{
		mov eax, a
		movd mm0, [eax]
		paddusb mm0, flashbrival
		movd [eax], mm0
	}
}

static _inline void mmxcolorsub (long *a)
{
	_asm
	{
		mov eax, a
		movd mm0, [eax]
		psubusb mm0, flashbrival
		movd [eax], mm0
	}
}

#endif

static _inline void addusb (char *a, long b)
{
	(*a) += b; if ((*a) < b) (*a) = 255;
}

#if (ESTNORMRAD == 2)
static signed char bitnum[32] =
{
	0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4,
	1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5
};
//static long bitsum[32] =
//{
//   0,-2,-1,-3, 0,-2,-1,-3, 1,-1, 0,-2, 1,-1, 0,-2,
//   2, 0, 1,-1, 2, 0, 1,-1, 3, 1, 2, 0, 3, 1, 2, 0
//};
static long bitsnum[32] =
{
	0        ,1-(2<<16),1-(1<<16),2-(3<<16),
	1        ,2-(2<<16),2-(1<<16),3-(3<<16),
	1+(1<<16),2-(1<<16),2        ,3-(2<<16),
	2+(1<<16),3-(1<<16),3        ,4-(2<<16),
	1+(2<<16),2        ,2+(1<<16),3-(1<<16),
	2+(2<<16),3        ,3+(1<<16),4-(1<<16),
	2+(3<<16),3+(1<<16),3+(2<<16),4,
	3+(3<<16),4+(1<<16),4+(2<<16),5
};
static float fsqrecip[5860]; //75*75 + 15*15 + 3*3 = 5859 is max value (5*5*5 box)
#endif

void estnorm (long x, long y, long z, point3d *fp)
{
	lpoint3d n;
	long *lptr, xx, yy, zz, b[5], i, j, k;
	float f;

	n.x = 0; n.y = 0; n.z = 0;

#if (ESTNORMRAD == 2)
	if (labs(x-xbsox) + labs(y-xbsoy) > 1)
	{
			//x,y not close enough to cache: calls expandbitstack 25 times :(
		xbsox = x; xbsoy = y; xbsof = 24*5;
		lptr = (long *)(&xbsbuf[24*5+1]);
		for(yy=-2;yy<=2;yy++)
			for(xx=-2;xx<=2;xx++,lptr-=10)
				expandbitstack(x+xx,y+yy,(__int64 *)lptr);
	}
	else if (x != xbsox)
	{
			//shift xbsbuf cache left/right: calls expandbitstack 5 times :)
		if (x < xbsox) { xx = -2; xbsof -= 24*5; lptr = (long *)(&xbsbuf[xbsof+1]); }
					 else { xx = 2; lptr = (long *)(&xbsbuf[xbsof-5*5+1]); xbsof -= 1*5; }
		xbsox = x; if (xbsof < 0) xbsof += 25*5;
		for(yy=-2;yy<=2;yy++)
		{
			if (lptr < (long *)&xbsbuf[1]) lptr += 25*10;
			expandbitstack(x+xx,y+yy,(__int64 *)lptr);
			lptr -= 5*10;
		}
	}
	else if (y != xbsoy)
	{
			//shift xbsbuf cache up/down: calls expandbitstack 5 times :)
		if (y < xbsoy) { yy = -2; xbsof -= 20*5; lptr = (long *)(&xbsbuf[xbsof+1]); }
					 else { yy = 2; lptr = (long *)(&xbsbuf[xbsof+1]); xbsof -= 5*5; }
		xbsoy = y; if (xbsof < 0) xbsof += 25*5;
		for(xx=-2;xx<=2;xx++)
		{
			if (lptr < (long *)&xbsbuf[1]) lptr += 25*10;
			expandbitstack(x+xx,y+yy,(__int64 *)lptr);
			lptr -= 1*10;
		}
	}

	z -= 2;
	if ((z&31) <= 27) //2 <= (z&31) <= 29
		{ lptr = (long *)((long)(&xbsbuf[xbsof+1]) + ((z&~31)>>3)); z &= 31; }
	else
		{ lptr = (long *)((long)(&xbsbuf[xbsof+1]) + (z>>3)); z &= 7; }

	for(yy=-2;yy<=2;yy++)
	{
		if (lptr >= (long *)&xbsbuf[1+10*5])
		{
			b[0] = ((lptr[  0]>>z)&31); b[1] = ((lptr[-10]>>z)&31);
			b[2] = ((lptr[-20]>>z)&31); b[3] = ((lptr[-30]>>z)&31);
			b[4] = ((lptr[-40]>>z)&31); lptr -= 50;
		}
		else
		{
			b[0] = ((lptr[0]>>z)&31); lptr -= 10; if (lptr < (long *)&xbsbuf[1]) lptr += 25*10;
			b[1] = ((lptr[0]>>z)&31); lptr -= 10; if (lptr < (long *)&xbsbuf[1]) lptr += 25*10;
			b[2] = ((lptr[0]>>z)&31); lptr -= 10; if (lptr < (long *)&xbsbuf[1]) lptr += 25*10;
			b[3] = ((lptr[0]>>z)&31); lptr -= 10; if (lptr < (long *)&xbsbuf[1]) lptr += 25*10;
			b[4] = ((lptr[0]>>z)&31); lptr -= 10; if (lptr < (long *)&xbsbuf[1]) lptr += 25*10;
		}

			//Make filter spherical
		//if (yy&1) { b[0] &= 0xe; b[4] &= 0xe; }
		//else if (yy) { b[0] &= 0x4; b[1] &= 0xe; b[3] &= 0xe; b[4] &= 0x4; }

		n.x += ((bitnum[b[4]]-bitnum[b[0]])<<1)+bitnum[b[3]]-bitnum[b[1]];
		j = bitsnum[b[0]]+bitsnum[b[1]]+bitsnum[b[2]]+bitsnum[b[3]]+bitsnum[b[4]];
		n.z += j; n.y += (*(signed short *)&j)*yy;
	}
	n.z >>= 16;
#else
	for(yy=-ESTNORMRAD;yy<=ESTNORMRAD;yy++)
		for(xx=-ESTNORMRAD;xx<=ESTNORMRAD;xx++)
			for(zz=-ESTNORMRAD;zz<=ESTNORMRAD;zz++)
				if (isvoxelsolid(x+xx,y+yy,z+zz))
					{ n.x += xx; n.y += yy; n.z += zz; }
#endif

#if 1
	f = fsqrecip[n.x*n.x + n.y*n.y + n.z*n.z];
	fp->x = ((float)n.x)*f; fp->y = ((float)n.y)*f; fp->z = ((float)n.z)*f;
#else

		//f = 1.0 / sqrt((double)(n.x*n.x + n.y*n.y + n.z*n.z));
		//fp->x = f*(float)n.x; fp->y = f*(float)n.y; fp->z = f*(float)n.z;
	zz = n.x*n.x + n.y*n.y + n.z*n.z;
	if (cputype&(1<<25))
	{
		_asm
		{
			cvtsi2ss xmm0, zz
			rsqrtss xmm0, xmm0
			;movss f, xmm0

				;fp->x = f*(float)n.x; fp->y = f*(float)n.y; fp->z = f*(float)n.z;
			cvtsi2ss xmm1, n.z
			shufps xmm0, xmm0, 0
			mov eax, fp
			movlhps xmm1, xmm1
			cvtpi2ps xmm1, n
			mulps xmm0, xmm1
			movlps [eax], xmm0
			movhlps xmm0, xmm0
			movss [eax+8], xmm0
		}
	}
	else
	{
		_asm
		{
			pi2fd mm0, zz       ;mm0:     0          zz
			pfrsqrt mm0, mm0    ;mm0: 1/sqrt(zz) 1/sqrt(zz)
			pi2fd mm1, n.x      ;mm1:     0         n.x
			pi2fd mm2, n.y      ;mm2:     0         n.y
			punpckldq mm1, mm2  ;mm1:    n.y        n.x
			pi2fd mm2, n.z      ;mm2:     0         n.z
			pfmul mm1, mm0      ;mm1:n.y/sqrt(zz) n.x/sqrt(zz)
			pfmul mm2, mm0      ;mm2:     0       n.z/sqrt(zz)
			mov eax, fp
			movq [eax], mm1
			movd [eax+8], mm2
			femms
		}
	}
#endif
}

static long vspan (long x, long y0, long y1)
{
	long y, yy, *bbufx;

	y = (y0>>5); bbufx = &bbuf[x][0];
	if ((y1>>5) == y)
	{
		yy = bbufx[y]; bbufx[y] &= ~(p2m[y1&31]^p2m[y0&31]);
		return(bbufx[y] ^ yy);
	}

	if (!(bbufx[y]&(~p2m[y0&31])))
		if (!(bbufx[y1>>5]&p2m[y1&31]))
		{
			for(yy=(y1>>5)-1;yy>y;yy--)
				if (bbufx[yy]) goto vspan_skip;
			return(0);
		}
vspan_skip:;
	bbufx[y] &= p2m[y0&31];
	bbufx[y1>>5] &= (~p2m[y1&31]);
	for(yy=(y1>>5)-1;yy>y;yy--) bbufx[yy] = 0;
	return(1);
}

static long docube (long x, long y, long z)
{
	long x0, y0, x1, y1, g;

	ffxptr = &ffx[(z+1)*z-1];
	x0 = (long)ffxptr[x].x; x1 = (long)ffxptr[x].y;
	y0 = (long)ffxptr[y].x; y1 = (long)ffxptr[y].y;
	for(g=0;x0<x1;x0++) g |= vspan(x0,y0,y1);
	return(g);
}

void hline (float x0, float y0, float x1, float y1, long *ix0, long *ix1)
{
	float dyx;

	dyx = (y1-y0) * grd; //grd = 1/(x1-x0)

		  if (y0 < wy0) ftol((wy0-y0)/dyx+x0,ix0);
	else if (y0 > wy1) ftol((wy1-y0)/dyx+x0,ix0);
	else ftol(x0,ix0);
		  if (y1 < wy0) ftol((wy0-y0)/dyx+x0,ix1);
	else if (y1 > wy1) ftol((wy1-y0)/dyx+x0,ix1);
	else ftol(x1,ix1);
	if ((*ix0) < iwx0) (*ix0) = iwx0;
	if ((*ix0) > iwx1) (*ix0) = iwx1; //(*ix1) = min(max(*ix1,wx0),wx1);
	gline(labs((*ix1)-(*ix0)),(float)(*ix0),((*ix0)-x1)*dyx + y1,
									  (float)(*ix1),((*ix1)-x1)*dyx + y1);
}

void vline (float x0, float y0, float x1, float y1, long *iy0, long *iy1)
{
	float dxy;

	dxy = (x1-x0) * grd; //grd = 1/(y1-y0)

		  if (x0 < wx0) ftol((wx0-x0)/dxy+y0,iy0);
	else if (x0 > wx1) ftol((wx1-x0)/dxy+y0,iy0);
	else ftol(y0,iy0);
		  if (x1 < wx0) ftol((wx0-x0)/dxy+y0,iy1);
	else if (x1 > wx1) ftol((wx1-x0)/dxy+y0,iy1);
	else ftol(y1,iy1);
	if ((*iy0) < iwy0) (*iy0) = iwy0;
	if ((*iy0) > iwy1) (*iy0) = iwy1;
	gline(labs((*iy1)-(*iy0)),((*iy0)-y1)*dxy + x1,(float)(*iy0),
									  ((*iy1)-y1)*dxy + x1,(float)(*iy1));
}

static float optistrx, optistry, optiheix, optiheiy, optiaddx, optiaddy;

static __int64 foglut[2048], fogcol;
static long ofogdist = -1;

#ifdef _MSC_VER

#ifdef __cplusplus
extern "C" {
#endif
extern void *opti4asm;
#define opti4 ((point4d *)&opti4asm)
#ifdef __cplusplus
}
#endif

void (*hrend)(long,long,long,long,long,long);
void (*vrend)(long,long,long,long,long);

#if (USEZBUFFER != 1)
void hrendnoz (long sx, long sy, long p1, long plc, long incr, long j)
{
	sy = ylookup[sy]+frameplace; p1 = sy+(p1<<2); sy += (sx<<2);
	do
	{
		*(long *)sy = angstart[plc>>16][j].col;
		plc += incr; sy += 4;
	} while (sy != p1);
}

void vrendnoz (long sx, long sy, long p1, long iplc, long iinc)
{
	sy = ylookup[sy]+(sx<<2)+frameplace;
	for(;sx<p1;sx++)
	{
		*(long *)sy = angstart[uurend[sx]>>16][iplc].col;
		uurend[sx] += uurend[sx+MAXXDIM]; sy += 4; iplc += iinc;
	}
}

#else

#if 0
	//Example C code
void hrendz (long sx, long sy, long p1, long plc, long incr, long j)
{
	long p0, i; float dirx, diry;
	p0 = ylookup[sy]+(sx<<2)+frameplace;
	p1 = ylookup[sy]+(p1<<2)+frameplace;
	dirx = optistrx*(float)sx + optiheix*(float)sy + optiaddx;
	diry = optistry*(float)sx + optiheiy*(float)sy + optiaddy;
	i = zbufoff;
	do
	{
		*(long *)p0 = angstart[plc>>16][j].col;
		*(float *)(p0+i) = (float)angstart[plc>>16][j].dist/sqrt(dirx*dirx+diry*diry);
		dirx += optistrx; diry += optistry; plc += incr; p0 += 4;
	} while (p0 != p1);
}

	//Example C code
void vrendz (long sx, long sy, long p1, long iplc, long iinc)
{
	float dirx, diry; long i, p0;
	p0 = ylookup[sy]+(sx<<2)+frameplace;
	p1 = ylookup[sy]+(p1<<2)+frameplace;
	dirx = optistrx*(float)sx + optiheix*(float)sy + optiaddx;
	diry = optistry*(float)sx + optiheiy*(float)sy + optiaddy;
	i = zbufoff;
	while (p0 < p1)
	{
		*(long *)p0 = angstart[uurend[sx]>>16][iplc].col;
		*(float *)(p0+i) = (float)angstart[uurend[sx]>>16][iplc].dist/sqrt(dirx*dirx+diry*diry);
		dirx += optistrx; diry += optistry; uurend[sx] += uurend[sx+MAXXDIM]; p0 += 4; iplc += iinc; sx++;
	}
}

	//Example C code
void hrendzfog (long sx, long sy, long p1, long plc, long incr, long j)
{
	long p0, i, k, l; float dirx, diry;
	p0 = ylookup[sy]+(sx<<2)+frameplace;
	p1 = ylookup[sy]+(p1<<2)+frameplace;
	dirx = optistrx*(float)sx + optiheix*(float)sy + optiaddx;
	diry = optistry*(float)sx + optiheiy*(float)sy + optiaddy;
	i = zbufoff;
	do
	{
		k = angstart[plc>>16][j].col;
		l = angstart[plc>>16][j].dist;
		l = (foglut[l>>20]&32767);
		*(long *)p0 = ((((( vx5.fogcol     &255)-( k     &255))*l)>>15)    ) +
						  ((((((vx5.fogcol>> 8)&255)-((k>> 8)&255))*l)>>15)<< 8) +
						  ((((((vx5.fogcol>>16)&255)-((k>>16)&255))*l)>>15)<<16)+k;
		*(float *)(p0+i) = (float)angstart[plc>>16][j].dist/sqrt(dirx*dirx+diry*diry);
		dirx += optistrx; diry += optistry; plc += incr; p0 += 4;
	} while (p0 != p1);
}

	//Example C code
void vrendzfog (long sx, long sy, long p1, long iplc, long iinc)
{
	float dirx, diry; long i, k, l, p0;
	p0 = ylookup[sy]+(sx<<2)+frameplace;
	p1 = ylookup[sy]+(p1<<2)+frameplace;
	dirx = optistrx*(float)sx + optiheix*(float)sy + optiaddx;
	diry = optistry*(float)sx + optiheiy*(float)sy + optiaddy;
	i = zbufoff;
	while (p0 < p1)
	{
		k = angstart[uurend[sx]>>16][iplc].col;
		l = angstart[uurend[sx]>>16][iplc].dist;
		l = (foglut[l>>20]&32767);
		*(long *)p0 = ((((( vx5.fogcol     &255)-( k     &255))*l)>>15)    ) +
						  ((((((vx5.fogcol>> 8)&255)-((k>> 8)&255))*l)>>15)<< 8) +
						  ((((((vx5.fogcol>>16)&255)-((k>>16)&255))*l)>>15)<<16)+k;
		*(float *)(p0+i) = (float)angstart[uurend[sx]>>16][iplc].dist/sqrt(dirx*dirx+diry*diry);
		dirx += optistrx; diry += optistry; uurend[sx] += uurend[sx+MAXXDIM]; p0 += 4; iplc += iinc; sx++;
	}
}

#endif

void hrendzsse (long sx, long sy, long p1, long plc, long incr, long j)
{
	_asm
	{
		push esi
		push edi
beghasm_p3:
		mov eax, sx
		mov ecx, sy
		mov esi, p1
		mov edx, ylookup[ecx*4]
		add edx, frameplace
		lea edi, [edx+eax*4]
		lea esi, [edx+esi*4]

		and eax, 0xfffffffc
		cvtsi2ss xmm0, eax
		cvtsi2ss xmm4, ecx
		movss xmm1, xmm0
		movss xmm5, xmm4
		mulss xmm0, optistrx
		mulss xmm1, optistry
		mulss xmm4, optiheix
		mulss xmm5, optiheiy
		addss xmm0, optiaddx
		addss xmm1, optiaddy
		addss xmm0, xmm4
		addss xmm1, xmm5

		mov ecx, zbufoff
		mov edx, j
		movd mm6, plc
		movd mm7, incr

		shufps xmm0, xmm0, 0
		shufps xmm1, xmm1, 0
		movaps xmm2, opti4asm[2*16]
		movaps xmm3, opti4asm[3*16]
		addps xmm0, opti4asm[0*16]
		addps xmm1, opti4asm[1*16]
			;xmm0 =  xmm0      ^2 +  xmm1      ^2        (p)
			;xmm2 = (xmm0+xmm2)^2 + (xmm1+xmm3)^2 - xmm0 (v)
			;xmm1 = ...                                  (a)
		addps xmm2, xmm0  ;This block converts inner loop...
		addps xmm3, xmm1  ;from: 1 / sqrt(x*x + y*y), x += xi, y += yi;
		mulps xmm0, xmm0  ;  to: 1 / sqrt(p), p += v, v += a;
		mulps xmm1, xmm1
		mulps xmm2, xmm2
		mulps xmm3, xmm3
		addps xmm0, xmm1
		movaps xmm1, opti4asm[4*16]
		addps xmm2, xmm3
		subps xmm2, xmm0

			;Do first 0-3 pixels to align unrolled loop of 4
		test edi, 15
		jz short skip1ha

		test edi, 8
		jz short skipshufa
		shufps xmm0, xmm0, 0x4e ;rotate right by 2
skipshufa:
		test edi, 4
		jz short skipshufb
		shufps xmm0, xmm0, 0x39 ;rotate right by 1
skipshufb:

beg1ha:
		pextrw eax, mm6, 1
		paddd mm6, mm7
		mov eax, angstart[eax*4]
		movd mm0, [eax+edx*8]
		movd [edi], mm0
		cvtsi2ss xmm7, [eax+edx*8+4]
		rsqrtss xmm3, xmm0
		mulss xmm7, xmm3
		shufps xmm0, xmm0, 0x39 ;rotate right by 1
		movss [edi+ecx], xmm7
		add edi, 4
		cmp edi, esi
		jz short endh
		test edi, 15
		jnz short beg1ha

		addps xmm0, xmm2
		addps xmm2, xmm1
skip1ha:
		lea eax, [edi+16]      ;these 3 lines re-ordered
		cmp eax, esi
		ja short skip4h

		movq mm0, mm6          ;mm0: 0,plc
		paddd mm0, mm7         ;mm0: 0,plc+inc
		punpckldq mm7, mm7     ;mm7: inc,inc
		punpckldq mm6, mm0     ;mm6: plc+inc,plc
		paddd mm7, mm7         ;mm7: inc+inc,inc+inc

		sub esi, 16

		 ;eax: temp   ³ mm0:  z0 argb0   argb1 argb0 ³ xmm0: plc3 plc2 plc1 plc0
		 ;ebx:  -     ³ mm1:  z1 argb1               ³ xmm1: acc3 acc2 acc1 acc0
		 ;ecx:zbufoff ³ mm2:  z2 argb2   argb3 argb2 ³ xmm2: inc3 inc2 inc1 inc0
		 ;edx:  j     ³ mm3:  z3 argb3               ³ xmm3:  r3   r2   r1   r0
		 ;esi:  -     ³ mm4:              z1    z0   ³ xmm4:            z3   z2
		 ;edi:scroff  ³ mm5:              z3    z2   ³ xmm5:
		 ;ebp:  -     ³ mm6: plc1 plc0               ³ xmm6:
beg4h: ;esp:  -     ³ mm7: inc1 inc0               ³ xmm7:  z3   z2   z1   z0
		pextrw eax, mm6, 1
		mov eax, angstart[eax*4]
		movq mm0, [eax+edx*8]
		pextrw eax, mm6, 3
		mov eax, angstart[eax*4]
		movq mm1, [eax+edx*8]
		paddd mm6, mm7
		pextrw eax, mm6, 1
		mov eax, angstart[eax*4]
		movq mm2, [eax+edx*8]
		pextrw eax, mm6, 3
		mov eax, angstart[eax*4]
		movq mm3, [eax+edx*8]
		paddd mm6, mm7

		movq mm4, mm0
		movq mm5, mm2
		punpckldq mm0, mm1
		punpckldq mm2, mm3
		movntq [edi], mm0
		movntq [edi+8], mm2

		punpckhdq mm4, mm1
		punpckhdq mm5, mm3
		cvtpi2ps xmm7, mm4
		cvtpi2ps xmm4, mm5
		rsqrtps xmm3, xmm0
		movlhps xmm7, xmm4
		mulps xmm7, xmm3
		movntps [edi+ecx], xmm7
		addps xmm0, xmm2
		addps xmm2, xmm1

		add edi, 16
		cmp edi, esi
		jbe short beg4h
		add esi, 16
		cmp edi, esi
		jae endh

		psrad mm7, 1    ;Restore mm7 from incr*2 to just incr for single loop
skip4h:
beg1h:
		pextrw eax, mm6, 1
		paddd mm6, mm7
		mov eax, angstart[eax*4]
		movd mm0, [eax+edx*8]
		movd [edi], mm0
		cvtsi2ss xmm7, [eax+edx*8+4]
		rsqrtss xmm3, xmm0
		mulss xmm7, xmm3
		shufps xmm0, xmm0, 0x39 ;rotate right by 1
		movss [edi+ecx], xmm7
		add edi, 4
		cmp edi, esi
		jb short beg1h
endh: pop edi
		pop esi
	}
}

void hrendzfogsse (long sx, long sy, long p1, long plc, long incr, long j)
{
	static __int64 mm7bak;
	_asm
	{
		push esi
		push edi
beghasm_p3:
		mov eax, sx
		mov ecx, sy
		mov esi, p1
		mov edx, ylookup[ecx*4]
		add edx, frameplace
		lea edi, [edx+eax*4]
		lea esi, [edx+esi*4]

		and eax, 0xfffffffc
		cvtsi2ss xmm0, eax
		cvtsi2ss xmm4, ecx
		movss xmm1, xmm0
		movss xmm5, xmm4
		mulss xmm0, optistrx
		mulss xmm1, optistry
		mulss xmm4, optiheix
		mulss xmm5, optiheiy
		addss xmm0, optiaddx
		addss xmm1, optiaddy
		addss xmm0, xmm4
		addss xmm1, xmm5

		mov ecx, zbufoff
		mov edx, j
		movd mm6, plc
		movd mm7, incr

		shufps xmm0, xmm0, 0
		shufps xmm1, xmm1, 0
		movaps xmm2, opti4asm[2*16]
		movaps xmm3, opti4asm[3*16]
		addps xmm0, opti4asm[0*16]
		addps xmm1, opti4asm[1*16]
			;xmm0 =  xmm0      ^2 +  xmm1      ^2        (p)
			;xmm2 = (xmm0+xmm2)^2 + (xmm1+xmm3)^2 - xmm0 (v)
			;xmm1 = ...                                  (a)
		addps xmm2, xmm0  ;This block converts inner loop...
		addps xmm3, xmm1  ;from: 1 / sqrt(x*x + y*y), x += xi, y += yi;
		mulps xmm0, xmm0  ;  to: 1 / sqrt(p), p += v, v += a;
		mulps xmm1, xmm1
		mulps xmm2, xmm2
		mulps xmm3, xmm3
		addps xmm0, xmm1
		movaps xmm1, opti4asm[4*16]
		addps xmm2, xmm3
		subps xmm2, xmm0

			;Do first 0-3 pixels to align unrolled loop of 4
		test edi, 15
		jz short skip1ha

		test edi, 8
		jz short skipshufa
		shufps xmm0, xmm0, 0x4e ;rotate right by 2
skipshufa:
		test edi, 4
		jz short skipshufb
		shufps xmm0, xmm0, 0x39 ;rotate right by 1
skipshufb:

beg1ha:
		pextrw eax, mm6, 1
		paddd mm6, mm7
		mov eax, angstart[eax*4]

			;Z
		cvtsi2ss xmm7, [eax+edx*8+4]
		rsqrtss xmm3, xmm0
		mulss xmm7, xmm3
		shufps xmm0, xmm0, 0x39 ;rotate right by 1
		movss [edi+ecx], xmm7

			;Col
		punpcklbw mm0, [eax+edx*8]
		psrlw mm0, 8
		movq mm1, fogcol
		psubw mm1, mm0
		paddw mm1, mm1
		mov eax, [eax+edx*8+4]
		shr eax, 16+4
		pmulhw mm1, foglut[eax*8]
		paddw mm0, mm1
		packuswb mm0, mm1
		movd [edi], mm0

		add edi, 4
		cmp edi, esi
		jz short endh
		test edi, 15
		jnz short beg1ha

		addps xmm0, xmm2
		addps xmm2, xmm1
skip1ha:
		lea eax, [edi+16]      ;these 3 lines re-ordered
		cmp eax, esi
		ja short skip4h

		movq mm0, mm6          ;mm0: 0,plc
		paddd mm0, mm7         ;mm0: 0,plc+inc
		punpckldq mm7, mm7     ;mm7: inc,inc
		punpckldq mm6, mm0     ;mm6: plc+inc,plc
		paddd mm7, mm7         ;mm7: inc+inc,inc+inc

		sub esi, 16

		 ;eax: temp   ³ mm0:  z0 argb0   argb1 argb0 ³ xmm0: plc3 plc2 plc1 plc0
		 ;ebx:  -     ³ mm1:  z1 argb1               ³ xmm1: acc3 acc2 acc1 acc0
		 ;ecx:zbufoff ³ mm2:  z2 argb2   argb3 argb2 ³ xmm2: inc3 inc2 inc1 inc0
		 ;edx:  j     ³ mm3:  z3 argb3               ³ xmm3:  r3   r2   r1   r0
		 ;esi:  -     ³ mm4:              z1    z0   ³ xmm4:            z3   z2
		 ;edi:scroff  ³ mm5:              z3    z2   ³ xmm5:
		 ;ebp:  -     ³ mm6: plc1 plc0               ³ xmm6:
		 ;esp:  -     ³ mm7: inc1 inc0               ³ xmm7:  z3   z2   z1   z0

		movq mm7bak, mm7
beg4h:pextrw eax, mm6, 1
		mov eax, angstart[eax*4]
		movq mm4, [eax+edx*8]
		pextrw eax, mm6, 3
		mov eax, angstart[eax*4]
		movq mm1, [eax+edx*8]
		paddd mm6, mm7bak
		pextrw eax, mm6, 1
		mov eax, angstart[eax*4]
		movq mm5, [eax+edx*8]
		pextrw eax, mm6, 3
		mov eax, angstart[eax*4]
		movq mm3, [eax+edx*8]
		paddd mm6, mm7bak

		movq mm0, mm4
		movq mm2, mm5

			;Do Z
		punpckhdq mm4, mm1
		punpckhdq mm5, mm3
		cvtpi2ps xmm7, mm4
		cvtpi2ps xmm4, mm5
		rsqrtps xmm3, xmm0
		movlhps xmm7, xmm4
		mulps xmm7, xmm3
		movntps [edi+ecx], xmm7
		addps xmm0, xmm2
		addps xmm2, xmm1

			;Do colors
			;mm4:dist1 dist0
			;mm5:dist3 dist2
		pxor mm7, mm7
		punpcklbw mm0, mm7
		punpcklbw mm1, mm7
		punpcklbw mm2, mm7
		punpcklbw mm3, mm7

		movq mm7, fogcol
		psubw mm7, mm0
		paddw mm7, mm7
		pextrw eax, mm4, 1
		shr eax, 4
		pmulhw mm7, foglut[eax*8]
		paddw mm0, mm7

		movq mm7, fogcol
		psubw mm7, mm1
		paddw mm7, mm7
		pextrw eax, mm4, 3
		shr eax, 4
		pmulhw mm7, foglut[eax*8]
		paddw mm1, mm7

		movq mm7, fogcol
		psubw mm7, mm2
		paddw mm7, mm7
		pextrw eax, mm5, 1
		shr eax, 4
		pmulhw mm7, foglut[eax*8]
		paddw mm2, mm7

		movq mm7, fogcol
		psubw mm7, mm3
		paddw mm7, mm7
		pextrw eax, mm5, 3
		shr eax, 4
		pmulhw mm7, foglut[eax*8]
		paddw mm3, mm7

		packuswb mm0, mm1
		packuswb mm2, mm3
		movntq [edi], mm0
		movntq [edi+8], mm2

		add edi, 16
		cmp edi, esi
		jbe short beg4h
		add esi, 16
		cmp edi, esi
		jae endh

		movq mm7, mm7bak
		psrad mm7, 1    ;Restore mm7 from incr*2 to just incr for single loop
skip4h:
beg1h:
		pextrw eax, mm6, 1
		paddd mm6, mm7
		mov eax, angstart[eax*4]

			;Z
		cvtsi2ss xmm7, [eax+edx*8+4]
		rsqrtss xmm3, xmm0
		mulss xmm7, xmm3
		shufps xmm0, xmm0, 0x39 ;rotate right by 1
		movss [edi+ecx], xmm7

			;Col
		punpcklbw mm0, [eax+edx*8]
		psrlw mm0, 8
		movq mm1, fogcol
		psubw mm1, mm0
		paddw mm1, mm1
		mov eax, [eax+edx*8+4]
		shr eax, 16+4
		pmulhw mm1, foglut[eax*8]
		paddw mm0, mm1
		packuswb mm0, mm1
		movd [edi], mm0

		add edi, 4
		cmp edi, esi
		jb short beg1h
endh: pop edi
		pop esi
	}
}

void hrendz3dn (long sx, long sy, long p1, long plc, long incr, long j)
{
	_asm
	{
		push esi
		push edi
		mov eax, sy
		mov eax, ylookup[eax*4]
		add eax, frameplace
		mov esi, p1
		lea esi, [eax+esi*4]    ;esi = p1
		mov edi, sx
		lea edi, [eax+edi*4]    ;edi = p0

		movd mm0, sx
		punpckldq mm0, sy
		pi2fd mm0, mm0          ;mm0: (float)sy (float)sx
		pshufw mm2, mm0, 0xee   ;mm2: (float)sy (float)sy
		punpckldq mm0, mm0      ;mm0: (float)sx (float)sx
		movd mm1, optistrx
		punpckldq mm1, optistry
		pfmul mm0, mm1          ;mm0: (float)sx*optistry (float)sx*optistrx
		movd mm3, optiheix
		punpckldq mm3, optiheiy
		pfmul mm2, mm3          ;mm2: (float)sy*optiheiy (float)sy*optiheix
		pfadd mm0, mm2
		movd mm3, optiaddx
		punpckldq mm3, optiaddy ;mm3: optiaddy optiaddx
		pfadd mm0, mm3          ;mm0: diry diry

		movd mm6, plc
		movd mm7, incr
		mov ecx, zbufoff
		mov edx, j

beg:  pextrw eax, mm6, 1
		mov eax, angstart[eax*4]
		movq mm2, [eax+edx*8]   ;mm2:      dist       col
		pshufw mm3, mm2, 0xee   ;mm3:         ?      dist
		pi2fd mm3, mm3          ;mm3:         ?   (f)dist
		movq mm4, mm0           ;mm4:      diry      dirx
		pfmul mm4, mm4          ;mm4:    diry^2    dirx^2
		pfadd mm0, mm1          ;mm0: dirx+optx diry+opty (unrelated)
		pfacc mm4, mm4          ;mm4: (x^2+y^2)   x^2+y^2
		pfrsqrt mm4, mm4        ;mm4: 1/sqrt(*) 1/sqrt(*)
		pfmul mm3, mm4          ;mm3:         0    zvalue
		paddd mm6, mm7          ;mm6:            plc+incr (unrelated)
		movd [edi], mm2
		movd [edi+ecx], mm3
		add edi, 4
		cmp edi, esi
		jb short beg
		pop edi
		pop esi
	}
}

void hrendzfog3dn (long sx, long sy, long p1, long plc, long incr, long j)
{
	_asm
	{
		push esi
		push edi
		mov eax, sy
		mov eax, ylookup[eax*4]
		add eax, frameplace
		mov esi, p1
		lea esi, [eax+esi*4]    ;esi = p1
		mov edi, sx
		lea edi, [eax+edi*4]    ;edi = p0

		movd mm0, sx
		punpckldq mm0, sy
		pi2fd mm0, mm0          ;mm0: (float)sy (float)sx
		pshufw mm2, mm0, 0xee   ;mm2: (float)sy (float)sy
		punpckldq mm0, mm0      ;mm0: (float)sx (float)sx
		movd mm1, optistrx
		punpckldq mm1, optistry
		pfmul mm0, mm1          ;mm0: (float)sx*optistry (float)sx*optistrx
		movd mm3, optiheix
		punpckldq mm3, optiheiy
		pfmul mm2, mm3          ;mm2: (float)sy*optiheiy (float)sy*optiheix
		pfadd mm0, mm2
		movd mm3, optiaddx
		punpckldq mm3, optiaddy ;mm3: optiaddy optiaddx
		pfadd mm0, mm3          ;mm0: diry diry

		pxor mm5, mm5

		movd mm6, plc
		movd mm7, incr
		mov ecx, zbufoff
		mov edx, j

beg:  pextrw eax, mm6, 1
		mov eax, angstart[eax*4]
		movq mm2, [eax+edx*8]   ;mm2:      dist       col
		pshufw mm3, mm2, 0xee   ;mm3:         ?      dist
		pi2fd mm3, mm3          ;mm3:         ?   (f)dist
		movq mm4, mm0           ;mm4:      diry      dirx
		pfmul mm4, mm4          ;mm4:    diry^2    dirx^2
		pfadd mm0, mm1          ;mm0: dirx+optx diry+opty (unrelated)
		pfacc mm4, mm4          ;mm4: (x^2+y^2)   x^2+y^2
		pfrsqrt mm4, mm4        ;mm4: 1/sqrt(*) 1/sqrt(*)
		pfmul mm3, mm4          ;mm3:         0    zvalue
		paddd mm6, mm7          ;mm6:            plc+incr (unrelated)

			;Extra calculations for fog
		pextrw eax, mm2, 3
		punpcklbw mm2, mm5
		movq mm4, fogcol
		psubw mm4, mm2
		paddw mm4, mm4
		shr eax, 4
		pmulhw mm4, foglut[eax*8]
		paddw mm2, mm4
		packuswb mm2, mm4

		movd [edi], mm2
		movd [edi+ecx], mm3
		add edi, 4
		cmp edi, esi
		jb short beg
		pop edi
		pop esi
	}
}

void vrendzsse (long sx, long sy, long p1, long iplc, long iinc)
{
	_asm
	{
		push ebx
		push esi
		push edi
begvasm_p3:
		mov esi, sx
		mov eax, sy
		mov edx, p1
		mov ecx, ylookup[eax*4]
		add ecx, frameplace
		lea edx, [ecx+edx*4]
		lea edi, [ecx+esi*4]

		mov ecx, esi
		and ecx, 0xfffffffc
		cvtsi2ss xmm0, ecx
		cvtsi2ss xmm4, eax
		movss xmm1, xmm0
		movss xmm5, xmm4
		mulss xmm0, optistrx
		mulss xmm1, optistry
		mulss xmm4, optiheix
		mulss xmm5, optiheiy
		addss xmm0, optiaddx
		addss xmm1, optiaddy
		addss xmm0, xmm4
		addss xmm1, xmm5

		shufps xmm0, xmm0, 0
		shufps xmm1, xmm1, 0
		movaps xmm2, opti4asm[2*16]
		movaps xmm3, opti4asm[3*16]
		addps xmm0, opti4asm[0*16]
		addps xmm1, opti4asm[1*16]
			;xmm0 =  xmm0      ^2 +  xmm1      ^2        (p)
			;xmm2 = (xmm0+xmm2)^2 + (xmm1+xmm3)^2 - xmm0 (v)
			;xmm1 = ...                                  (a)
		addps xmm2, xmm0  ;This block converts inner loop...
		addps xmm3, xmm1  ;from: 1 / sqrt(x*x + y*y), x += xi, y += yi;
		mulps xmm0, xmm0  ;  to: 1 / sqrt(p), p += v, v += a;
		mulps xmm1, xmm1
		mulps xmm2, xmm2
		mulps xmm3, xmm3
		addps xmm0, xmm1
		movaps xmm1, opti4asm[4*16]
		addps xmm2, xmm3
		subps xmm2, xmm0

		mov p1, edx
		mov ecx, zbufoff
		shl esi, 2
		add esi, uurend
		mov ebx, iplc

		cmp edi, edx
		jae short endv

			;Do first 0-3 pixels to align unrolled loop of 4
		test edi, 15
		jz short skip1va

		test edi, 8
		jz short skipshufc
		shufps xmm0, xmm0, 0x4e ;rotate right by 2
skipshufc:
		test edi, 4
		jz short skipshufd
		shufps xmm0, xmm0, 0x39 ;rotate right by 1
skipshufd:

beg1va:
		mov edx, [esi]
		mov eax, [esi+MAXXDIM*4]
		add eax, edx
		sar edx, 16
		mov edx, angstart[edx*4]
		mov [esi], eax
		mov eax, [edx+ebx*8]
		mov [edi], eax
		cvtsi2ss xmm7, [edx+ebx*8+4]
		rsqrtss xmm3, xmm0
		mulss xmm7, xmm3
		shufps xmm0, xmm0, 0x39 ;rotate right by 1
		movss [edi+ecx], xmm7
		add ebx, iinc
		add esi, 4
		add edi, 4
		cmp edi, p1
		jz short endv
		test edi, 15
		jnz short beg1va

		addps xmm0, xmm2
		addps xmm2, xmm1
skip1va:
		lea edx, [edi+16]
		cmp edx, p1
		ja short prebeg1v

		cmp iinc, 0
		jl short beg4vn

beg4vp:
		movq mm6, [esi]
		movq mm7, [esi+8]
		pextrw eax, mm6, 1
		pextrw edx, mm6, 3
		paddd mm6, [esi+MAXXDIM*4]
		mov eax, angstart[eax*4]
		mov edx, angstart[edx*4]
		movq mm0, [eax+ebx*8]
		movq mm1, [edx+ebx*8+8]
		pextrw eax, mm7, 1
		pextrw edx, mm7, 3
		paddd mm7, [esi+8+MAXXDIM*4]
		mov eax, angstart[eax*4]
		mov edx, angstart[edx*4]
		movq mm2, [eax+ebx*8+16]
		movq mm3, [edx+ebx*8+24]
		add ebx, 4

		movq mm4, mm0
		movq mm5, mm2
		punpckldq mm0, mm1
		punpckldq mm2, mm3
		movntq [edi], mm0
		movntq [edi+8], mm2

		punpckhdq mm4, mm1
		punpckhdq mm5, mm3
		cvtpi2ps xmm7, mm4
		cvtpi2ps xmm4, mm5
		rsqrtps xmm3, xmm0
		movlhps xmm7, xmm4
		mulps xmm7, xmm3
		movntps [edi+ecx], xmm7
		addps xmm0, xmm2
		addps xmm2, xmm1

		movq [esi], mm6
		movq [esi+8], mm7

		add esi, 16
		add edi, 16
		lea edx, [edi+16]
		cmp edx, p1
		jbe short beg4vp
		cmp edi, p1
		jae short endv
		jmp short prebeg1v

beg4vn:
		movq mm6, [esi]
		movq mm7, [esi+8]
		pextrw eax, mm6, 1
		pextrw edx, mm6, 3
		paddd mm6, [esi+MAXXDIM*4]
		mov eax, angstart[eax*4]
		mov edx, angstart[edx*4]
		movq mm0, [eax+ebx*8]
		movq mm1, [edx+ebx*8-8]
		pextrw eax, mm7, 1
		pextrw edx, mm7, 3
		paddd mm7, [esi+8+MAXXDIM*4]
		mov eax, angstart[eax*4]
		mov edx, angstart[edx*4]
		movq mm2, [eax+ebx*8-16]
		movq mm3, [edx+ebx*8-24]
		sub ebx, 4

		movq mm4, mm0
		movq mm5, mm2
		punpckldq mm0, mm1
		punpckldq mm2, mm3
		movntq [edi], mm0
		movntq [edi+8], mm2

		punpckhdq mm4, mm1
		punpckhdq mm5, mm3
		cvtpi2ps xmm7, mm4
		cvtpi2ps xmm4, mm5
		rsqrtps xmm3, xmm0
		movlhps xmm7, xmm4
		mulps xmm7, xmm3
		movntps [edi+ecx], xmm7
		addps xmm0, xmm2
		addps xmm2, xmm1

		movq [esi], mm6
		movq [esi+8], mm7

		add esi, 16
		add edi, 16
		lea edx, [edi+16]
		cmp edx, p1
		jbe short beg4vn
		cmp edi, p1
		jae short endv

prebeg1v:
beg1v:
		mov edx, [esi]
		mov eax, [esi+MAXXDIM*4]
		add eax, edx
		sar edx, 16
		mov edx, angstart[edx*4]
		mov [esi], eax
		mov eax, [edx+ebx*8]
		mov [edi], eax
		cvtsi2ss xmm7, [edx+ebx*8+4]
		rsqrtss xmm3, xmm0
		mulss xmm7, xmm3
		shufps xmm0, xmm0, 0x39 ;rotate right by 1
		movss [edi+ecx], xmm7
		add ebx, iinc
		add esi, 4
		add edi, 4
		cmp edi, p1
		jne short beg1v
endv: pop edi
		pop esi
		pop ebx
	}
}

void vrendzfogsse (long sx, long sy, long p1, long iplc, long iinc)
{
	_asm
	{
		push ebx
		push esi
		push edi
begvasm_p3:
		mov esi, sx
		mov eax, sy
		mov edx, p1
		mov ecx, ylookup[eax*4]
		add ecx, frameplace
		lea edx, [ecx+edx*4]
		lea edi, [ecx+esi*4]

		mov ecx, esi
		and ecx, 0xfffffffc
		cvtsi2ss xmm0, ecx
		cvtsi2ss xmm4, eax
		movss xmm1, xmm0
		movss xmm5, xmm4
		mulss xmm0, optistrx
		mulss xmm1, optistry
		mulss xmm4, optiheix
		mulss xmm5, optiheiy
		addss xmm0, optiaddx
		addss xmm1, optiaddy
		addss xmm0, xmm4
		addss xmm1, xmm5

		shufps xmm0, xmm0, 0
		shufps xmm1, xmm1, 0
		movaps xmm2, opti4asm[2*16]
		movaps xmm3, opti4asm[3*16]
		addps xmm0, opti4asm[0*16]
		addps xmm1, opti4asm[1*16]
			;xmm0 =  xmm0      ^2 +  xmm1      ^2        (p)
			;xmm2 = (xmm0+xmm2)^2 + (xmm1+xmm3)^2 - xmm0 (v)
			;xmm1 = ...                                  (a)
		addps xmm2, xmm0  ;This block converts inner loop...
		addps xmm3, xmm1  ;from: 1 / sqrt(x*x + y*y), x += xi, y += yi;
		mulps xmm0, xmm0  ;  to: 1 / sqrt(p), p += v, v += a;
		mulps xmm1, xmm1
		mulps xmm2, xmm2
		mulps xmm3, xmm3
		addps xmm0, xmm1
		movaps xmm1, opti4asm[4*16]
		addps xmm2, xmm3
		subps xmm2, xmm0

		mov p1, edx
		mov ecx, zbufoff
		shl esi, 2
		add esi, uurend
		mov ebx, iplc

		cmp edi, edx
		jae short endv

			;Do first 0-3 pixels to align unrolled loop of 4
		test edi, 15
		jz short skip1va

		test edi, 8
		jz short skipshufc
		shufps xmm0, xmm0, 0x4e ;rotate right by 2
skipshufc:
		test edi, 4
		jz short skipshufd
		shufps xmm0, xmm0, 0x39 ;rotate right by 1
skipshufd:

beg1va:
		mov edx, [esi]
		mov eax, [esi+MAXXDIM*4]
		add eax, edx
		sar edx, 16
		mov edx, angstart[edx*4]
		mov [esi], eax

			;Z
		cvtsi2ss xmm7, [edx+ebx*8+4]
		rsqrtss xmm3, xmm0
		mulss xmm7, xmm3
		shufps xmm0, xmm0, 0x39 ;rotate right by 1
		movss [edi+ecx], xmm7

			;Col
		punpcklbw mm0, [edx+ebx*8]
		psrlw mm0, 8
		movq mm1, fogcol
		psubw mm1, mm0
		paddw mm1, mm1
		mov eax, [edx+ebx*8+4]
		shr eax, 16+4
		pmulhw mm1, foglut[eax*8]
		paddw mm0, mm1
		packuswb mm0, mm1
		movd [edi], mm0

		add ebx, iinc
		add esi, 4
		add edi, 4
		cmp edi, p1
		jz short endv
		test edi, 15
		jnz short beg1va

		addps xmm0, xmm2
		addps xmm2, xmm1
skip1va:
		lea edx, [edi+16]
		cmp edx, p1
		ja short prebeg1v

		cmp iinc, 0
		jl short beg4vn

beg4vp:
		movq mm6, [esi]
		movq mm7, [esi+8]
		pextrw eax, mm6, 1
		pextrw edx, mm6, 3
		paddd mm6, [esi+MAXXDIM*4]
		mov eax, angstart[eax*4]
		mov edx, angstart[edx*4]
		movq mm4, [eax+ebx*8]
		movq mm1, [edx+ebx*8+8]
		pextrw eax, mm7, 1
		pextrw edx, mm7, 3
		paddd mm7, [esi+8+MAXXDIM*4]
		mov eax, angstart[eax*4]
		mov edx, angstart[edx*4]
		movq mm5, [eax+ebx*8+16]
		movq mm3, [edx+ebx*8+24]
		add ebx, 4

			;Do Z
		movq mm0, mm4
		movq mm2, mm5
		punpckhdq mm4, mm1
		punpckhdq mm5, mm3
		cvtpi2ps xmm7, mm4
		cvtpi2ps xmm4, mm5
		rsqrtps xmm3, xmm0
		movlhps xmm7, xmm4
		mulps xmm7, xmm3
		movntps [edi+ecx], xmm7
		addps xmm0, xmm2
		addps xmm2, xmm1

		movq [esi], mm6
		movq [esi+8], mm7

			;Do color
		pxor mm7, mm7
		punpcklbw mm0, mm7
		punpcklbw mm1, mm7
		punpcklbw mm2, mm7
		punpcklbw mm3, mm7

		movq mm7, fogcol
		psubw mm7, mm0
		paddw mm7, mm7
		pextrw eax, mm4, 1
		shr eax, 4
		pmulhw mm7, foglut[eax*8]
		paddw mm0, mm7

		movq mm7, fogcol
		psubw mm7, mm1
		paddw mm7, mm7
		pextrw eax, mm4, 3
		shr eax, 4
		pmulhw mm7, foglut[eax*8]
		paddw mm1, mm7

		movq mm7, fogcol
		psubw mm7, mm2
		paddw mm7, mm7
		pextrw eax, mm5, 1
		shr eax, 4
		pmulhw mm7, foglut[eax*8]
		paddw mm2, mm7

		movq mm7, fogcol
		psubw mm7, mm3
		paddw mm7, mm7
		pextrw eax, mm5, 3
		shr eax, 4
		pmulhw mm7, foglut[eax*8]
		paddw mm3, mm7

		packuswb mm0, mm1
		packuswb mm2, mm3
		movntq [edi], mm0
		movntq [edi+8], mm2

		add esi, 16
		add edi, 16
		lea edx, [edi+16]
		cmp edx, p1
		jbe short beg4vp
		cmp edi, p1
		jae short endv
		jmp short prebeg1v

beg4vn:
		movq mm6, [esi]
		movq mm7, [esi+8]
		pextrw eax, mm6, 1
		pextrw edx, mm6, 3
		paddd mm6, [esi+MAXXDIM*4]
		mov eax, angstart[eax*4]
		mov edx, angstart[edx*4]
		movq mm4, [eax+ebx*8]
		movq mm1, [edx+ebx*8-8]
		pextrw eax, mm7, 1
		pextrw edx, mm7, 3
		paddd mm7, [esi+8+MAXXDIM*4]
		mov eax, angstart[eax*4]
		mov edx, angstart[edx*4]
		movq mm5, [eax+ebx*8-16]
		movq mm3, [edx+ebx*8-24]
		sub ebx, 4

			;Do Z
		movq mm0, mm4
		movq mm2, mm5
		punpckhdq mm4, mm1
		punpckhdq mm5, mm3
		cvtpi2ps xmm7, mm4
		cvtpi2ps xmm4, mm5
		rsqrtps xmm3, xmm0
		movlhps xmm7, xmm4
		mulps xmm7, xmm3
		movntps [edi+ecx], xmm7
		addps xmm0, xmm2
		addps xmm2, xmm1

		movq [esi], mm6
		movq [esi+8], mm7

			;Do color
		pxor mm7, mm7
		punpcklbw mm0, mm7
		punpcklbw mm1, mm7
		punpcklbw mm2, mm7
		punpcklbw mm3, mm7

		movq mm7, fogcol
		psubw mm7, mm0
		paddw mm7, mm7
		pextrw eax, mm4, 1
		shr eax, 4
		pmulhw mm7, foglut[eax*8]
		paddw mm0, mm7

		movq mm7, fogcol
		psubw mm7, mm1
		paddw mm7, mm7
		pextrw eax, mm4, 3
		shr eax, 4
		pmulhw mm7, foglut[eax*8]
		paddw mm1, mm7

		movq mm7, fogcol
		psubw mm7, mm2
		paddw mm7, mm7
		pextrw eax, mm5, 1
		shr eax, 4
		pmulhw mm7, foglut[eax*8]
		paddw mm2, mm7

		movq mm7, fogcol
		psubw mm7, mm3
		paddw mm7, mm7
		pextrw eax, mm5, 3
		shr eax, 4
		pmulhw mm7, foglut[eax*8]
		paddw mm3, mm7

		packuswb mm0, mm1
		packuswb mm2, mm3
		movntq [edi], mm0
		movntq [edi+8], mm2

		add esi, 16
		add edi, 16
		lea edx, [edi+16]
		cmp edx, p1
		jbe short beg4vn
		cmp edi, p1
		jae short endv

prebeg1v:
beg1v:
		mov edx, [esi]
		mov eax, [esi+MAXXDIM*4]
		add eax, edx
		sar edx, 16
		mov edx, angstart[edx*4]
		mov [esi], eax

			;Z
		cvtsi2ss xmm7, [edx+ebx*8+4]
		rsqrtss xmm3, xmm0
		mulss xmm7, xmm3
		shufps xmm0, xmm0, 0x39 ;rotate right by 1
		movss [edi+ecx], xmm7

			;Col
		punpcklbw mm0, [edx+ebx*8]
		psrlw mm0, 8
		movq mm1, fogcol
		psubw mm1, mm0
		paddw mm1, mm1
		mov eax, [edx+ebx*8+4]
		shr eax, 16+4
		pmulhw mm1, foglut[eax*8]
		paddw mm0, mm1
		packuswb mm0, mm1
		movd [edi], mm0

		add ebx, iinc
		add esi, 4
		add edi, 4
		cmp edi, p1
		jne short beg1v
endv: pop edi
		pop esi
		pop ebx
	}
}

void vrendz3dn (long sx, long sy, long p1, long iplc, long iinc)
{
	_asm
	{
		push ebx
		push esi
		push edi
		mov esi, p1
		mov edi, sx
		cmp edi, esi
		jae short endv
		mov eax, sy
		mov eax, ylookup[eax*4]
		add eax, frameplace
		lea esi, [eax+esi*4]    ;esi = p1
		lea edi, [eax+edi*4]    ;edi = p0

		movd mm0, sx
		punpckldq mm0, sy
		pi2fd mm0, mm0          ;mm0: (float)sy (float)sx
		pshufw mm2, mm0, 0xee   ;mm2: (float)sy (float)sy
		punpckldq mm0, mm0      ;mm0: (float)sx (float)sx
		movd mm1, optistrx
		punpckldq mm1, optistry
		pfmul mm0, mm1          ;mm0: (float)sx*optistry (float)sx*optistrx
		movd mm3, optiheix
		punpckldq mm3, optiheiy
		pfmul mm2, mm3          ;mm2: (float)sy*optiheiy (float)sy*optiheix
		pfadd mm0, mm2
		movd mm3, optiaddx
		punpckldq mm3, optiaddy ;mm3: optiaddy optiaddx
		pfadd mm0, mm3          ;mm0: diry diry

		mov ecx, zbufoff
		mov edx, iplc
		mov ebx, sx
		mov eax, uurend
		lea ebx, [eax+ebx*4]

begv_3dn:
		movd mm5, [ebx]
		pextrw eax, mm5, 1
		paddd mm5, [ebx+MAXXDIM*4]
		movd [ebx], mm5
		mov eax, angstart[eax*4]
		movq mm2, [eax+edx*8]   ;mm2:      dist       col
		pshufw mm3, mm2, 0xee   ;mm3:         ?      dist
		pi2fd mm3, mm3          ;mm3:         ?   (f)dist
		movq mm4, mm0           ;mm4:      diry      dirx
		pfmul mm4, mm4          ;mm4:    diry^2    dirx^2
		pfadd mm0, mm1          ;mm0: dirx+optx diry+opty (unrelated)
		pfacc mm4, mm4          ;mm4: (x^2+y^2)   x^2+y^2
		pfrsqrt mm4, mm4        ;mm4: 1/sqrt(*) 1/sqrt(*)
		pfmul mm3, mm4          ;mm3:         0    zvalue
		movd [edi], mm2
		movd [edi+ecx], mm3
		add edx, iinc
		add ebx, 4
		add edi, 4
		cmp edi, esi
		jb short begv_3dn
endv: pop edi
		pop esi
		pop ebx
	}
}

void vrendzfog3dn (long sx, long sy, long p1, long iplc, long iinc)
{
	_asm
	{
		push ebx
		push esi
		push edi
		mov esi, p1
		mov edi, sx
		cmp edi, esi
		jae short endv
		mov eax, sy
		mov eax, ylookup[eax*4]
		add eax, frameplace
		lea esi, [eax+esi*4]    ;esi = p1
		lea edi, [eax+edi*4]    ;edi = p0

		movd mm0, sx
		punpckldq mm0, sy
		pi2fd mm0, mm0          ;mm0: (float)sy (float)sx
		pshufw mm2, mm0, 0xee   ;mm2: (float)sy (float)sy
		punpckldq mm0, mm0      ;mm0: (float)sx (float)sx
		movd mm1, optistrx
		punpckldq mm1, optistry
		pfmul mm0, mm1          ;mm0: (float)sx*optistry (float)sx*optistrx
		movd mm3, optiheix
		punpckldq mm3, optiheiy
		pfmul mm2, mm3          ;mm2: (float)sy*optiheiy (float)sy*optiheix
		pfadd mm0, mm2
		movd mm3, optiaddx
		punpckldq mm3, optiaddy ;mm3: optiaddy optiaddx
		pfadd mm0, mm3          ;mm0: diry diry

		pxor mm6, mm6

		mov ecx, zbufoff
		mov edx, iplc
		mov ebx, sx
		mov eax, uurend
		lea ebx, [eax+ebx*4]

begv_3dn:
		movd mm5, [ebx]
		pextrw eax, mm5, 1
		paddd mm5, [ebx+MAXXDIM*4]
		movd [ebx], mm5
		mov eax, angstart[eax*4]
		movq mm2, [eax+edx*8]   ;mm2:      dist       col
		pshufw mm3, mm2, 0xee   ;mm3:         ?      dist
		pi2fd mm3, mm3          ;mm3:         ?   (f)dist
		movq mm4, mm0           ;mm4:      diry      dirx
		pfmul mm4, mm4          ;mm4:    diry^2    dirx^2
		pfadd mm0, mm1          ;mm0: dirx+optx diry+opty (unrelated)
		pfacc mm4, mm4          ;mm4: (x^2+y^2)   x^2+y^2
		pfrsqrt mm4, mm4        ;mm4: 1/sqrt(*) 1/sqrt(*)
		pfmul mm3, mm4          ;mm3:         0    zvalue

			;Extra calculations for fog
		pextrw eax, mm2, 3
		punpcklbw mm2, mm6
		movq mm4, fogcol
		psubw mm4, mm2
		paddw mm4, mm4
		shr eax, 4
		pmulhw mm4, foglut[eax*8]
		paddw mm2, mm4
		packuswb mm2, mm4

		movd [edi], mm2
		movd [edi+ecx], mm3
		add edx, iinc
		add ebx, 4
		add edi, 4
		cmp edi, esi
		jb short begv_3dn
endv: pop edi
		pop esi
		pop ebx
	}
}

#endif

void setcamera (dpoint3d *ipo, dpoint3d *ist, dpoint3d *ihe, dpoint3d *ifo,
					 float dahx, float dahy, float dahz)
{
	long i, j;

	gipos.x = ipo->x; gipos.y = ipo->y; gipos.z = ipo->z;
	gistr.x = ist->x; gistr.y = ist->y; gistr.z = ist->z;
	gihei.x = ihe->x; gihei.y = ihe->y; gihei.z = ihe->z;
	gifor.x = ifo->x; gifor.y = ifo->y; gifor.z = ifo->z;
	gihx = dahx; gihy = dahy; gihz = dahz;

	gixs.x = gistr.x; gixs.y = gihei.x; gixs.z = gifor.x;
	giys.x = gistr.y; giys.y = gihei.y; giys.z = gifor.y;
	gizs.x = gistr.z; gizs.y = gihei.z; gizs.z = gifor.z;
	giadd.x = -(gipos.x*gistr.x + gipos.y*gistr.y + gipos.z*gistr.z);
	giadd.y = -(gipos.x*gihei.x + gipos.y*gihei.y + gipos.z*gihei.z);
	giadd.z = -(gipos.x*gifor.x + gipos.y*gifor.y + gipos.z*gifor.z);

	gcorn[0].x = -gihx*gistr.x - gihy*gihei.x + gihz*gifor.x;
	gcorn[0].y = -gihx*gistr.y - gihy*gihei.y + gihz*gifor.y;
	gcorn[0].z = -gihx*gistr.z - gihy*gihei.z + gihz*gifor.z;
	gcorn[1].x = xres*gistr.x+gcorn[0].x;
	gcorn[1].y = xres*gistr.y+gcorn[0].y;
	gcorn[1].z = xres*gistr.z+gcorn[0].z;
	gcorn[2].x = yres*gihei.x+gcorn[1].x;
	gcorn[2].y = yres*gihei.y+gcorn[1].y;
	gcorn[2].z = yres*gihei.z+gcorn[1].z;
	gcorn[3].x = yres*gihei.x+gcorn[0].x;
	gcorn[3].y = yres*gihei.y+gcorn[0].y;
	gcorn[3].z = yres*gihei.z+gcorn[0].z;
	for(j=0,i=3;j<4;i=j++)
	{
		ginor[i].x = gcorn[i].y*gcorn[j].z - gcorn[i].z*gcorn[j].y;
		ginor[i].y = gcorn[i].z*gcorn[j].x - gcorn[i].x*gcorn[j].z;
		ginor[i].z = gcorn[i].x*gcorn[j].y - gcorn[i].y*gcorn[j].x;
	}
}

void opticast ()
{
	float f, ff, cx, cy, fx, fy, gx, gy, x0, y0, x1, y1, x2, y2, x3, y3;
	long i, j, sx, sy, p0, p1, cx16, cy16, kadd, kmul, u, u1, ui;

	if (gifor.z < 0) giforzsgn = -1; else giforzsgn = 1; //giforzsgn = (gifor.z < 0);

	gixyi[0] = (VSID<<2); gixyi[1] = -gixyi[0];
	glipos.x = ((long)gipos.x);
	glipos.y = ((long)gipos.y);
	glipos.z = ((long)gipos.z);
	gpixy = (long)&sptr[glipos.y*VSID + glipos.x];
	ftol(gipos.z*PREC-.5f,&gposz);
	gposxfrac[1] = gipos.x - (float)glipos.x; gposxfrac[0] = 1-gposxfrac[1];
	gposyfrac[1] = gipos.y - (float)glipos.y; gposyfrac[0] = 1-gposyfrac[1];
#if USEV5ASM
	for(j=u=0;j<gmipnum;j++,u+=i)
		for(i=0;i<(256>>j)+4;i++)
			gylookup[i+u] = ((((gposz>>j)-i*PREC)>>(16-j))&0x0000ffff);
	gxmip = max(vx5.mipscandist,4)*PREC;
#else
	for(i=0;i<256+4;i++) gylookup[i] = (i*PREC-gposz);
#endif
	gmaxscandist = min(max(vx5.maxscandist,1),2047)*PREC;

#if (USEZBUFFER != 1)
	hrend = hrendnoz; vrend = vrendnoz;
#else
	if (ofogdist < 0)
	{
		if (cputype&(1<<25)) { hrend = hrendzsse; vrend = vrendzsse; }
							 else { hrend = hrendz3dn; vrend = vrendz3dn; }
	}
	else
	{
		if (cputype&(1<<25)) { hrend = hrendzfogsse; vrend = vrendzfogsse; }
							 else { hrend = hrendzfog3dn; vrend = vrendzfog3dn; }

	}
#endif
	if (ofogdist < 0) nskypic = skypic;
				  else { nskypic = skyoff = 0; } //Optimization hack: draw sky as pure black when using fog

	gstartv = (char *)*(long *)gpixy;
	if (glipos.z >= gstartv[1])
	{
		do
		{
			if (!gstartv[0]) return;
			gstartv += gstartv[0]*4;
		} while (glipos.z >= gstartv[1]);
		if (glipos.z < gstartv[3]) return;
		gstartz0 = gstartv[3];
	} else gstartz0 = 0;
	gstartz1 = gstartv[1];

	if (gifor.z == 0) f = 32000; else f = gihz/gifor.z;
	f = min(max(f,-32000),32000);
	cx = gistr.z*f + gihx;
	cy = gihei.z*f + gihy;

	wx0 = (float)(-(vx5.anginc)); wx1 = (float)(xres-1+(vx5.anginc));
	wy0 = (float)(-(vx5.anginc)); wy1 = (float)(yres-1+(vx5.anginc));
	ftol(wx0,&iwx0); ftol(wx1,&iwx1);
	ftol(wy0,&iwy0); ftol(wy1,&iwy1);

	fx = wx0-cx; fy = wy0-cy; gx = wx1-cx; gy = wy1-cy;
	x0 = x3 = wx0; y0 = y1 = wy0; x1 = x2 = wx1; y2 = y3 = wy1;
	if (fy < 0)
	{
		if (fx < 0) { f = sqrt(fx*fy); x0 = cx-f; y0 = cy-f; }
		if (gx > 0) { f = sqrt(-gx*fy); x1 = cx+f; y1 = cy-f; }
	}
	if (gy > 0)
	{
		if (gx > 0) { f = sqrt(gx*gy); x2 = cx+f; y2 = cy+f; }
		if (fx < 0) { f = sqrt(-fx*gy); x3 = cx-f; y3 = cy+f; }
	}
	if (x0 > x1) { if (fx < 0) y0 = fx/gx*fy + cy; else y1 = gx/fx*fy + cy; }
	if (y1 > y2) { if (fy < 0) x1 = fy/gy*gx + cx; else x2 = gy/fy*gx + cx; }
	if (x2 < x3) { if (fx < 0) y3 = fx/gx*gy + cy; else y2 = gx/fx*gy + cy; }
	if (y3 < y0) { if (fy < 0) x0 = fy/gy*fx + cx; else x3 = gy/fy*fx + cx; }
		//This makes precision errors cause pixels to overwrite rather than omit
	x0 -= .01; x1 += .01;
	y1 -= .01; y2 += .01;
	x3 -= .01; x2 += .01;
	y0 -= .01; y3 += .01;

	f = (float)PREC / gihz;
	optistrx = gistr.x*f; optiheix = gihei.x*f; optiaddx = gcorn[0].x*f;
	optistry = gistr.y*f; optiheiy = gihei.y*f; optiaddy = gcorn[0].y*f;
#ifdef _MSC_VER
	opti4[0].y = optistrx; opti4[0].z = optistrx*2; opti4[0].z2 = optistrx*3;
	opti4[1].y = optistry; opti4[1].z = optistry*2; opti4[1].z2 = optistry*3;
	opti4[2].x = opti4[2].y = opti4[2].z = opti4[2].z2 = optistrx*4.0f;
	opti4[3].x = opti4[3].y = opti4[3].z = opti4[3].z2 = optistry*4.0f;
	opti4[4].x = opti4[4].y = opti4[4].z = opti4[4].z2 = (optistrx*optistrx + optistry*optistry)*32.0f; //NEW ALGO!
#endif

	ftol(cx*65536,&cx16);
	ftol(cy*65536,&cy16);

	ftol((x1-x0)/vx5.anginc,&j);
	if ((fy < 0) && (j > 0)) //(cx,cy),(x0,wy0),(x1,wy0)
	{
		ff = (x1-x0) / (float)j; grd = 1.0f / (wy0-cy);
		gscanptr = (castdat *)radar; skycurlng = -1; skycurdir = -giforzsgn;
		for(i=0,f=x0+ff*.5f;i<j;f+=ff,i++)
		{
			vline(cx,cy,f,wy0,&p0,&p1);
			if (giforzsgn < 0) angstart[i] = gscanptr+p0; else angstart[i] = gscanptr-p1;
			gscanptr += labs(p1-p0)+1;
		}

		j <<= 16; f = (float)j / ((x1-x0)*grd); ftol((cx-x0)*grd*f,&kadd);
		ftol(cx-.5f,&p1); p0 = lbound0(p1+1,xres); p1 = lbound0(p1,xres);
		ftol(cy-0.50005f,&sy); if (sy >= yres) sy = yres-1;
		ff = (fabs((float)p1-cx)+1)*f/2147483647.0 + cy; //Anti-crash hack
		while ((ff < sy) && (sy >= 0)) sy--;
		if (sy >= 0)
		{
			ftol(f,&kmul);
			for(;sy>=0;sy--) if (isshldiv16safe(kmul,(sy<<16)-cy16)) break; //Anti-crash hack
			if (giforzsgn < 0) i = -sy; else i = sy;
			for(;sy>=0;sy--,i-=giforzsgn)
			{
				ui = shldiv16(kmul,(sy<<16)-cy16);
				u = mulshr16((p0<<16)-cx16,ui)+kadd;
				while ((p0 > 0) && (u >= ui)) { u -= ui; p0--; }
				u1 = (p1-p0)*ui + u;
				while ((p1 < xres) && (u1 < j)) { u1 += ui; p1++; }
				if (p0 < p1) hrend(p0,sy,p1,u,ui,i);
			}
			_asm emms
		}
	}

	ftol((y2-y1)/vx5.anginc,&j);
	if ((gx > 0) && (j > 0)) //(cx,cy),(wx1,y1),(wx1,y2)
	{
		ff = (y2-y1) / (float)j; grd = 1.0f / (wx1-cx);
		gscanptr = (castdat *)radar; skycurlng = -1; skycurdir = -giforzsgn;
		for(i=0,f=y1+ff*.5f;i<j;f+=ff,i++)
		{
			hline(cx,cy,wx1,f,&p0,&p1);
			if (giforzsgn < 0) angstart[i] = gscanptr-p0; else angstart[i] = gscanptr+p1;
			gscanptr += labs(p1-p0)+1;
		}

		j <<= 16; f = (float)j / ((y2-y1)*grd); ftol((cy-y1)*grd*f,&kadd);
		ftol(cy-.5f,&p1); p0 = lbound0(p1+1,yres); p1 = lbound0(p1,yres);
		ftol(cx+0.50005f,&sx); if (sx < 0) sx = 0;
		ff = (fabs((float)p1-cy)+1)*f/2147483647.0 + cx; //Anti-crash hack
		while ((ff > sx) && (sx < xres)) sx++;
		if (sx < xres)
		{
			ftol(f,&kmul);
			for(;sx<xres;sx++) if (isshldiv16safe(kmul,(sx<<16)-cx16)) break; //Anti-crash hack
			for(;sx<xres;sx++)
			{
				ui = shldiv16(kmul,(sx<<16)-cx16);
				u = mulshr16((p0<<16)-cy16,ui)+kadd;
				while ((p0 > 0) && (u >= ui)) { u -= ui; lastx[--p0] = sx; }
				uurend[sx] = u; uurend[sx+MAXXDIM] = ui; u += (p1-p0)*ui;
				while ((p1 < yres) && (u < j)) { u += ui; lastx[p1++] = sx; }
			}
			if (giforzsgn < 0)
				  { for(sy=p0;sy<p1;sy++) vrend(lastx[sy],sy,xres,lastx[sy],1); }
			else { for(sy=p0;sy<p1;sy++) vrend(lastx[sy],sy,xres,-lastx[sy],-1); }
			_asm emms
		}
	}

	ftol((x2-x3)/vx5.anginc,&j);
	if ((gy > 0) && (j > 0)) //(cx,cy),(x2,wy1),(x3,wy1)
	{
		ff = (x2-x3) / (float)j; grd = 1.0f / (wy1-cy);
		gscanptr = (castdat *)radar; skycurlng = -1; skycurdir = giforzsgn;
		for(i=0,f=x3+ff*.5f;i<j;f+=ff,i++)
		{
			vline(cx,cy,f,wy1,&p0,&p1);
			if (giforzsgn < 0) angstart[i] = gscanptr-p0; else angstart[i] = gscanptr+p1;
			gscanptr += labs(p1-p0)+1;
		}

		j <<= 16; f = (float)j / ((x2-x3)*grd); ftol((cx-x3)*grd*f,&kadd);
		ftol(cx-.5f,&p1); p0 = lbound0(p1+1,xres); p1 = lbound0(p1,xres);
		ftol(cy+0.50005f,&sy); if (sy < 0) sy = 0;
		ff = (fabs((float)p1-cx)+1)*f/2147483647.0 + cy; //Anti-crash hack
		while ((ff > sy) && (sy < yres)) sy++;
		if (sy < yres)
		{
			ftol(f,&kmul);
			for(;sy<yres;sy++) if (isshldiv16safe(kmul,(sy<<16)-cy16)) break; //Anti-crash hack
			if (giforzsgn < 0) i = sy; else i = -sy;
			for(;sy<yres;sy++,i-=giforzsgn)
			{
				ui = shldiv16(kmul,(sy<<16)-cy16);
				u = mulshr16((p0<<16)-cx16,ui)+kadd;
				while ((p0 > 0) && (u >= ui)) { u -= ui; p0--; }
				u1 = (p1-p0)*ui + u;
				while ((p1 < xres) && (u1 < j)) { u1 += ui; p1++; }
				if (p0 < p1) hrend(p0,sy,p1,u,ui,i);
			}
			_asm emms
		}
	}

	ftol((y3-y0)/vx5.anginc,&j);
	if ((fx < 0) && (j > 0)) //(cx,cy),(wx0,y3),(wx0,y0)
	{
		ff = (y3-y0) / (float)j; grd = 1.0f / (wx0-cx);
		gscanptr = (castdat *)radar; skycurlng = -1; skycurdir = giforzsgn;
		for(i=0,f=y0+ff*.5f;i<j;f+=ff,i++)
		{
			hline(cx,cy,wx0,f,&p0,&p1);
			if (giforzsgn < 0) angstart[i] = gscanptr+p0; else angstart[i] = gscanptr-p1;
			gscanptr += labs(p1-p0)+1;
		}

		j <<= 16; f = (float)j / ((y3-y0)*grd); ftol((cy-y0)*grd*f,&kadd);
		ftol(cy-.5f,&p1); p0 = lbound0(p1+1,yres); p1 = lbound0(p1,yres);
		ftol(cx-0.50005f,&sx); if (sx >= xres) sx = xres-1;
		ff = (fabs((float)p1-cy)+1)*f/2147483647.0 + cx; //Anti-crash hack
		while ((ff < sx) && (sx >= 0)) sx--;
		if (sx >= 0)
		{
			ftol(f,&kmul);
			for(;sx>=0;sx--) if (isshldiv16safe(kmul,(sx<<16)-cx16)) break; //Anti-crash hack
			for(;sx>=0;sx--)
			{
				ui = shldiv16(kmul,(sx<<16)-cx16);
				u = mulshr16((p0<<16)-cy16,ui)+kadd;
				while ((p0 > 0) && (u >= ui)) { u -= ui; lastx[--p0] = sx; }
				uurend[sx] = u; uurend[sx+MAXXDIM] = ui; u += (p1-p0)*ui;
				while ((p1 < yres) && (u < j)) { u += ui; lastx[p1++] = sx; }
			}
			for(sy=p0;sy<p1;sy++) vrend(0,sy,lastx[sy]+1,0,giforzsgn);
			_asm emms
		}
	}
}

	//0: asm temp for current x
	//1: asm temp for current y
	//2: bottom (28)
	//3: top    ( 0)
	//4: left   ( 8)
	//5: right  (24)
	//6: up     (12)
	//7: down   (12)
	//setsideshades(0,0,0,0,0,0);
	//setsideshades(0,28,8,24,12,12);
void setsideshades (char sto, char sbo, char sle, char sri, char sup, char sdo)
{
	((char *)&gcsub[2])[7] = sbo; ((char *)&gcsub[3])[7] = sto;
	((char *)&gcsub[4])[7] = sle; ((char *)&gcsub[5])[7] = sri;
	((char *)&gcsub[6])[7] = sup; ((char *)&gcsub[7])[7] = sdo;
	if (!(sto|sbo|sle|sri|sup|sdo))
	{
		vx5.sideshademode = 0;
		((char *)&gcsub[0])[7] = ((char *)&gcsub[1])[7] = 0x00;
	}
	else vx5.sideshademode = 1;
}

	//MUST have more than: CEILING(max possible CLIPRADIUS) * 4 entries!
#define MAXCLIPIT (VSID*4) //VSID*2+4 is not a power of 2!
static lpoint2d clipit[MAXCLIPIT];

double findmaxcr (double px, double py, double pz, double cr)
{
	double f, g, maxcr, thresh2;
	long x, y, z, i0, i1, ix, y0, y1, z0, z1;
	char *v;

	thresh2 = cr+1.7321+1; thresh2 *= thresh2;
	maxcr = cr*cr;

		//Find closest point of all nearby cubes to (px,py,pz)
	x = (long)px; y = (long)py; z = (long)pz; i0 = i1 = 0; ix = x; y0 = y1 = y;
	while (1)
	{
		f = max(fabs((double)x+.5-px)-.5,0);
		g = max(fabs((double)y+.5-py)-.5,0);
		f = f*f + g*g;
		if (f < maxcr)
		{
			if (((unsigned long)x >= VSID) || ((unsigned long)y >= VSID))
				{ z0 = z1 = 0; }
			else
			{
				v = sptr[y*VSID+x];
				if (z >= v[1])
				{
					while (1)
					{
						if (!v[0]) { z0 = z1 = 0; break; }
						v += v[0]*4;
						if (z < v[1]) { z0 = v[3]; z1 = v[1]; break; }
					}
				}
				else { z0 = MAXZDIM-2048; z1 = v[1]; }
			}

			if ((pz <= z0) || (pz >= z1))
				maxcr = f;
			else
			{
				g = min(pz-(double)z0,(double)z1-pz);
				f += g*g; if (f < maxcr) maxcr = f;
			}
		}

		if ((x-px)*(x-px)+(y-py)*(y-py) < thresh2)
		{
			if ((x <= ix) && (x > 0))
				{ clipit[i1].x = x-1; clipit[i1].y = y; i1 = ((i1+1)&(MAXCLIPIT-1)); }
			if ((x >= ix) && (x < VSID-1))
				{ clipit[i1].x = x+1; clipit[i1].y = y; i1 = ((i1+1)&(MAXCLIPIT-1)); }
			if ((y <= y0) && (y > 0))
				{ clipit[i1].x = x; clipit[i1].y = y-1; i1 = ((i1+1)&(MAXCLIPIT-1)); y0--; }
			if ((y >= y1) && (y < VSID-1))
				{ clipit[i1].x = x; clipit[i1].y = y+1; i1 = ((i1+1)&(MAXCLIPIT-1)); y1++; }
		}
		if (i0 == i1) break;
		x = clipit[i0].x; y = clipit[i0].y; i0 = ((i0+1)&(MAXCLIPIT-1));
	}
	return(sqrt(maxcr));
}

static double gx0, gy0, gcrf2, grdst, gendt, gux, guy;
static long gdist2square (double x, double y)
{
	double t;
	x -= gx0; y -= gy0; t = x*gux + y*guy; if (t <= 0) t = gcrf2;
	else if (t*grdst >= gendt) { x -= gux*gendt; y -= guy*gendt; t = gcrf2; }
	else t = t*t*grdst + gcrf2;
	return(x*x + y*y <= t);
}

long sphtrace (double x0, double y0, double z0,          //start pt
					double vx, double vy, double vz,          //move vector
					double *hitx, double *hity, double *hitz, //new pt after collision
					double *clpx, double *clpy, double *clpz, //pt causing collision
					double cr, double acr)
{
	double f, t, dax, day, daz, vyx, vxy, vxz, vyz, rvz, cr2, fz, fc;
	double dx, dy, dx1, dy1;
	double nx, ny, intx, inty, intz, dxy, dxz, dyz, dxyz, rxy, rxz, ryz, rxyz;
	long i, j, x, y, ix, iy0, iy1, i0, i1, iz[2], cz0, cz1;
	char *v;

		 //Precalculate global constants for ins & getval functions
	if ((vx == 0) && (vy == 0) && (vz == 0))
		{ (*hitx) = x0; (*hity) = y0; (*hitz) = z0; return(1); }
	gux = vx; guy = vy; gx0 = x0; gy0 = y0; dxy = vx*vx + vy*vy;
	if (dxy != 0) rxy = 1.0 / dxy; else rxy = 0;
	grdst = rxy; gendt = 1; cr2 = cr*cr; t = cr + 0.7072; gcrf2 = t*t;

	if (((long *)&vz)[1] >= 0) { dtol(   z0-cr-.5,&cz0); dtol(vz+z0+cr-.5,&cz1); }
								 else { dtol(vz+z0-cr-.5,&cz0); dtol(   z0+cr-.5,&cz1); }

		//Precalculate stuff for closest point on cube finder
	dax = 0; day = 0; vyx = 0; vxy = 0; rvz = 0; vxz = 0; vyz = 0;
	if (vx != 0) { vyx = vy/vx; if (((long *)&vx)[1] >= 0) dax = x0+cr; else dax = x0-cr-1; }
	if (vy != 0) { vxy = vx/vy; if (((long *)&vy)[1] >= 0) day = y0+cr; else day = y0-cr-1; }
	if (vz != 0)
	{
		rvz = 1.0/vz; vxz = vx*rvz; vyz = vy*rvz;
		if (((long *)&vz)[1] >= 0) daz = z0+cr; else daz = z0-cr;
	}

	dxyz = vz*vz;
	dxz = vx*vx+dxyz; if (dxz != 0) rxz = 1.0 / dxz;
	dyz = vy*vy+dxyz; if (dyz != 0) ryz = 1.0 / dyz;
	dxyz += dxy; rxyz = 1.0 / dxyz;

	dtol(x0-.5,&x); dtol(y0-.5,&y);
	ix = x; iy0 = iy1 = y;
	i0 = 0; clipit[0].x = x; clipit[0].y = y; i1 = 1;
	do
	{
		x = clipit[i0].x; y = clipit[i0].y; i0 = ((i0+1)&(MAXCLIPIT-1));

		dx = (double)x; dx1 = (double)(x+1);
		dy = (double)y; dy1 = (double)(y+1);

			//closest point on cube finder
			//Plane intersection (both vertical planes)
#if 0
		intx = dbound((dy-day)*vxy + x0,dx,dx1);
		inty = dbound((dx-dax)*vyx + y0,dy,dy1);
#else
		intx = (dy-day)*vxy + x0;
		inty = (dx-dax)*vyx + y0;
		if (((long *)&intx)[1] < ((long *)&dx)[1]) intx = dx;
		if (((long *)&inty)[1] < ((long *)&dy)[1]) inty = dy;
		if (((long *)&intx)[1] >= ((long *)&dx1)[1]) intx = dx1;
		if (((long *)&inty)[1] >= ((long *)&dy1)[1]) inty = dy1;
		//if (intx < (double)x) intx = (double)x;
		//if (inty < (double)y) inty = (double)y;
		//if (intx > (double)(x+1)) intx = (double)(x+1);
		//if (inty > (double)(y+1)) inty = (double)(y+1);
#endif

		do
		{
			if (((long *)&dxy)[1] == 0) { t = -1.0; continue; }
			nx = intx-x0; ny = inty-y0; t = vx*nx + vy*ny; if (((long *)&t)[1] < 0) continue;
			f = cr2 - nx*nx - ny*ny; if (((long *)&f)[1] >= 0) { t = -1.0; continue; }
			f = f*dxy + t*t; if (((long *)&f)[1] < 0) { t = -1.0; continue; }
			t = (t-sqrt(f))*rxy;
		} while (0);
		if (t >= gendt) goto sphtracecont;
		if (((long *)&t)[1] < 0) intz = z0; else intz = vz*t + z0;

			//Find closest ceil(iz[0]) & flor(iz[1]) in (x,y) column
		dtol(intz-.5,&i);
		if ((unsigned long)(x|y) < VSID)
		{
			v = sptr[y*VSID+x]; iz[0] = MAXZDIM-2048; iz[1] = v[1];
			while (i >= iz[1])
			{
				if (!v[0]) { iz[1] = -1; break; }
				v += v[0]*4;
				iz[0] = v[3]; if (i < iz[0]) { iz[1] = -1; break; }
				iz[1] = v[1];
			}
		}
		else iz[1] = -1;

			//hit xz plane, yz plane or z-axis edge?
		if (iz[1] < 0) //Treat whole column as solid
		{
			if (((long *)&t)[1] >= 0) { gendt = t; (*clpx) = intx; (*clpy) = inty; (*clpz) = intz; goto sphtracecont; }
		}

			//Must check tops & bottoms of slab
		for(i=1;i>=0;i--)
		{
				//Ceil/flor outside of quick&dirty bounding box
			if ((iz[i] < cz0) || (iz[i] > cz1)) continue;

				//Plane intersection (parallel to ground)
			intz = (double)iz[i]; t = intz-daz;
			intx = t*vxz + x0;
			inty = t*vyz + y0;

			j = 0;                         // A ³ 8 ³ 9
			//     if (intx < dx)  j |= 2; //ÄÄÄÅÄÄÄÅÄÄÄ
			//else if (intx > dx1) j |= 1; // 2 ³ 0 ³ 1
			//     if (inty < dy)  j |= 8; //ÄÄÄÅÄÄÄÅÄÄÄ
			//else if (inty > dy1) j |= 4; // 6 ³ 4 ³ 5
				  if (((long *)&intx)[1] <  ((long *)&dx)[1])  j |= 2;
			else if (((long *)&intx)[1] >= ((long *)&dx1)[1]) j |= 1;
				  if (((long *)&inty)[1] <  ((long *)&dy)[1])  j |= 8;
			else if (((long *)&inty)[1] >= ((long *)&dy1)[1]) j |= 4;

				//NOTE: only need to check once per "for"!
			if ((!j) && (vz != 0)) //hit xy plane?
			{
				t *= rvz;
				if ((((long *)&t)[1] >= 0) && (t < gendt)) { gendt = t; (*clpx) = intx; (*clpy) = inty; (*clpz) = intz; }
				continue;
			}

				//common calculations used for rest of checks...
			fz = intz-z0; fc = cr2-fz*fz; fz *= vz;

			if (j&3)
			{
				nx = (double)((j&1)+x);
				if (((long *)&dxz)[1] != 0) //hit y-axis edge?
				{
					f = nx-x0; t = vx*f + fz; f = (fc - f*f)*dxz + t*t;
					if (((long *)&f)[1] >= 0) t = (t-sqrt(f))*rxz; else t = -1.0;
				} else t = -1.0;
				ny = vy*t + y0;
					  if (((long *)&ny)[1] > ((long *)&dy1)[1]) j |= 0x10;
				else if (((long *)&ny)[1] >= ((long *)&dy)[1])
				{
					if ((((long *)&t)[1] >= 0) && (t < gendt)) { gendt = t; (*clpx) = nx; (*clpy) = ny; (*clpz) = intz; }
					continue;
				}
				inty = (double)(((j>>4)&1)+y);
			}
			else inty = (double)(((j>>2)&1)+y);

			if (j&12)
			{
				ny = (double)(((j>>2)&1)+y);
				if (((long *)&dyz)[1] != 0) //hit x-axis edge?
				{
					f = ny-y0; t = vy*f + fz; f = (fc - f*f)*dyz + t*t;
					if (((long *)&f)[1] >= 0) t = (t-sqrt(f))*ryz; else t = -1.0;
				} else t = -1.0;
				nx = vx*t + x0;
					  if (((long *)&nx)[1] > ((long *)&dx1)[1]) j |= 0x20;
				else if (((long *)&nx)[1] >= ((long *)&dx)[1])
				{
					if ((((long *)&t)[1] >= 0) && (t < gendt)) { gendt = t; (*clpx) = nx; (*clpy) = ny; (*clpz) = intz; }
					continue;
				}
				intx = (double)(((j>>5)&1)+x);
			}
			else intx = (double)((j&1)+x);

				//hit corner?
			nx = intx-x0; ny = inty-y0;
			t = vx*nx + vy*ny + fz; if (((long *)&t)[1] < 0) continue;
			f = fc - nx*nx - ny*ny; if (((long *)&f)[1] >= 0) continue;
			f = f*dxyz + t*t; if (((long *)&f)[1] < 0) continue;
			t = (t-sqrt(f))*rxyz;
			if (t < gendt) { gendt = t; (*clpx) = intx; (*clpy) = inty; (*clpz) = intz; }
		}
sphtracecont:;
		if ((x <= ix)  && (x >      0) && (gdist2square(dx- .5,dy+ .5))) { clipit[i1].x = x-1; clipit[i1].y = y; i1 = ((i1+1)&(MAXCLIPIT-1)); }
		if ((x >= ix)  && (x < VSID-1) && (gdist2square(dx+1.5,dy+ .5))) { clipit[i1].x = x+1; clipit[i1].y = y; i1 = ((i1+1)&(MAXCLIPIT-1)); }
		if ((y <= iy0) && (y >      0) && (gdist2square(dx+ .5,dy- .5))) { clipit[i1].x = x; clipit[i1].y = y-1; i1 = ((i1+1)&(MAXCLIPIT-1)); iy0 = y-1; }
		if ((y >= iy1) && (y < VSID-1) && (gdist2square(dx+ .5,dy+1.5))) { clipit[i1].x = x; clipit[i1].y = y+1; i1 = ((i1+1)&(MAXCLIPIT-1)); iy1 = y+1; }
	} while (i0 != i1);
#if 1
	(*hitx) = dbound(vx*gendt + x0,acr,VSID-acr);
	(*hity) = dbound(vy*gendt + y0,acr,VSID-acr);
	(*hitz) = dbound(vz*gendt + z0,MAXZDIM-2048+acr,MAXZDIM-1-acr);
#else
	(*hitx) = min(max(vx*gendt + x0,acr),VSID-acr);
	(*hity) = min(max(vy*gendt + y0,acr),VSID-acr);
	(*hitz) = min(max(vz*gendt + z0,MAXZDIM-2048+acr),MAXZDIM-1-acr);
#endif
	return(gendt == 1);
}

void clipmove (dpoint3d *p, dpoint3d *v, double acr)
{
	double f, gx, gy, gz, nx, ny, nz, ex, ey, ez, hitx, hity, hitz, cr;
	//double nx2, ny2, nz2, ex2, ey2, ez2; //double ox, oy, oz;
	long i, j, k;

	//ox = p->x; oy = p->y; oz = p->z;
	gx = p->x+v->x; gy = p->y+v->y; gz = p->z+v->z;

	cr = findmaxcr(p->x,p->y,p->z,acr);
	vx5.clipmaxcr = cr;

	vx5.cliphitnum = 0;
	for(i=0;i<3;i++)
	{
		if ((v->x == 0) && (v->y == 0) && (v->z == 0)) break;

		cr -= 1e-7;  //Shrinking radius error control hack

		//j = sphtraceo(p->x,p->y,p->z,v->x,v->y,v->z,&nx,&ny,&nz,&ex,&ey,&ez,cr,acr);
		//k = sphtraceo(p->x,p->y,p->z,v->x,v->y,v->z,&nx2,&ny2,&nz2,&ex2,&ey2,&ez2,cr,acr);

		j = sphtrace(p->x,p->y,p->z,v->x,v->y,v->z,&nx,&ny,&nz,&ex,&ey,&ez,cr,acr);

		//if ((j != k) || (fabs(nx-nx2) > .000001) || (fabs(ny-ny2) > .000001) || (fabs(nz-nz2) > .000001) ||
		//   ((j == 0) && ((fabs(ex-ex2) > .000001) || (fabs(ey-ey2) > .000001) || (fabs(ez-ez2) > .000001))))
		//{
		//   printf("%d %f %f %f %f %f %f\n",i,p->x,p->y,p->z,v->x,v->y,v->z);
		//   printf("%f %f %f ",nx,ny,nz); if (!j) printf("%f %f %f\n",ex,ey,ez); else printf("\n");
		//   printf("%f %f %f ",nx2,ny2,nz2); if (!k) printf("%f %f %f\n",ex2,ey2,ez2); else printf("\n");
		//   printf("\n");
		//}
		if (j) { p->x = nx; p->y = ny; p->z = nz; break; }

		vx5.cliphit[i].x = ex; vx5.cliphit[i].y = ey; vx5.cliphit[i].z = ez;
		vx5.cliphitnum = i+1;
		p->x = nx; p->y = ny; p->z = nz;

			//Calculate slide vector
		v->x = gx-nx; v->y = gy-ny; v->z = gz-nz;
		switch(i)
		{
			case 0:
				hitx = ex-nx; hity = ey-ny; hitz = ez-nz;
				f = (v->x*hitx + v->y*hity + v->z*hitz) / (cr * cr);
				v->x -= hitx*f; v->y -= hity*f; v->z -= hitz*f;
				break;
			case 1:
				nx -= ex; ny -= ey; nz -= ez;
				ex = hitz*ny - hity*nz;
				ey = hitx*nz - hitz*nx;
				ez = hity*nx - hitx*ny;
				f = ex*ex + ey*ey + ez*ez; if (f <= 0) break;
				f = (v->x*ex + v->y*ey + v->z*ez) / f;
				v->x = ex*f; v->y = ey*f; v->z = ez*f;
				break;
			default: break;
		}
	}

		//If you didn't move much, then don't move at all. This helps prevents
		//cliprad from shrinking, but you get stuck too much :(
	//if ((p->x-ox)*(p->x-ox) + (p->y-oy)*(p->y-oy) + (p->z-oz)*(p->z-oz) < 1e-12)
	//   { p->x = ox; p->y = oy; p->z = oz; }
}

unsigned long calcglobalmass ()
{
	unsigned long i, j;
	char *v;

	j = VSID*VSID*256;
	for(i=0;i<VSID*VSID;i++)
	{
		v = sptr[i]; j -= v[1];
		while (v[0]) { v += v[0]*4; j += v[3]-v[1]; }
	}
	return(j);
}

long loadvxl (const char *lodfilnam, dpoint3d *ipo, dpoint3d *ist, dpoint3d *ihe, dpoint3d *ifo)
{
	FILE *fil;
	long i, j, fsiz;
	char *v, *v2;

	if (!vbuf) { vbuf = (long *)malloc((VOXSIZ>>2)<<2); if (!vbuf) evilquit("vbuf malloc failed"); }
	if (!vbit) { vbit = (long *)malloc((VOXSIZ>>7)<<2); if (!vbit) evilquit("vbuf malloc failed"); }

	if (!kzopen(lodfilnam)) return(0);
	fsiz = kzfilelength();

	kzread(&i,4); if (i != 0x09072000) return(0);
	kzread(&i,4); if (i != VSID) return(0);
	kzread(&i,4); if (i != VSID) return(0);
	kzread(ipo,24);
	kzread(ist,24);
	kzread(ihe,24);
	kzread(ifo,24);

	v = (char *)(&vbuf[1]); //1st dword for voxalloc compare logic optimization
	kzread((void *)v,fsiz-kztell());

	for(i=0;i<VSID*VSID;i++)
	{
		sptr[i] = v;
		while (v[0]) v += (((long)v[0])<<2);
		v += ((((long)v[2])-((long)v[1])+2)<<2);
	}
	kzclose();

	memset(&sptr[VSID*VSID],0,sizeof(sptr)-VSID*VSID*4);
	vbiti = (((long)v-(long)vbuf)>>2); //# vbuf longs/vbit bits allocated
	clearbuf((void *)vbit,vbiti>>5,-1);
	clearbuf((void *)&vbit[vbiti>>5],(VOXSIZ>>7)-(vbiti>>5),0);
	vbit[vbiti>>5] = (1<<vbiti)-1;

	vx5.globalmass = calcglobalmass();
	backedup = -1;

	gmipnum = 1; vx5.flstnum = 0;
	updatebbox(0,0,0,VSID,VSID,MAXZDIM,0);
	return(1);
}

void dorthorotate (double ox, double oy, double oz, dpoint3d *ist, dpoint3d *ihe, dpoint3d *ifo)
{
	double f, t, dx, dy, dz, rr[9];

	dcossin(ox,&ox,&dx);
	dcossin(oy,&oy,&dy);
	dcossin(oz,&oz,&dz);
	f = ox*oz; t = dx*dz; rr[0] =  t*dy + f; rr[7] = -f*dy - t;
	f = ox*dz; t = dx*oz; rr[1] = -f*dy + t; rr[6] =  t*dy - f;
	rr[2] = dz*oy; rr[3] = -dx*oy; rr[4] = ox*oy; rr[8] = oz*oy; rr[5] = dy;
	ox = ist->x; oy = ihe->x; oz = ifo->x;
	ist->x = ox*rr[0] + oy*rr[3] + oz*rr[6];
	ihe->x = ox*rr[1] + oy*rr[4] + oz*rr[7];
	ifo->x = ox*rr[2] + oy*rr[5] + oz*rr[8];
	ox = ist->y; oy = ihe->y; oz = ifo->y;
	ist->y = ox*rr[0] + oy*rr[3] + oz*rr[6];
	ihe->y = ox*rr[1] + oy*rr[4] + oz*rr[7];
	ifo->y = ox*rr[2] + oy*rr[5] + oz*rr[8];
	ox = ist->z; oy = ihe->z; oz = ifo->z;
	ist->z = ox*rr[0] + oy*rr[3] + oz*rr[6];
	ihe->z = ox*rr[1] + oy*rr[4] + oz*rr[7];
	ifo->z = ox*rr[2] + oy*rr[5] + oz*rr[8];
}

static long mixc[MAXZDIM>>1][8]; //4K
static long mixn[MAXZDIM>>1];    //0.5K
void genmipvxl (long x0, long y0, long x1, long y1)
{
	long i, n, oldn, x, y, z, xsiz, ysiz, zsiz, oxsiz, oysiz;
	long cz, oz, nz, zz, besti, cstat, curz[4], curzn[4][4], mipnum, mipmax;
	char *v[4], *tv, **sr, **sw, **ssr, **ssw;

	if ((!(x0|y0)) && (x1 == VSID) && (y1 == VSID)) mipmax = vx5.vxlmipuse;
															 else mipmax = gmipnum;
	if (mipmax <= 0) return;
	mipnum = 1;

	vx5.colfunc = curcolfunc;
	xsiz = VSID; ysiz = VSID; zsiz = MAXZDIM;
	ssr = sptr; ssw = sptr+xsiz*ysiz;
	while ((xsiz > 1) && (ysiz > 1) && (zsiz > 1) && (mipnum < mipmax))
	{
		oxsiz = xsiz; xsiz >>= 1;
		oysiz = ysiz; ysiz >>= 1;
						  zsiz >>= 1;

		x0--; if (x0 < 0) x0 = 0;
		y0--; if (y0 < 0) y0 = 0;
		x1++; if (x1 > VSID) x1 = VSID;
		y1++; if (y1 > VSID) y1 = VSID;

		x0 >>= 1; x1 = ((x1+1)>>1);
		y0 >>= 1; y1 = ((y1+1)>>1);
		for(y=y0;y<y1;y++)
		{
			sr = ssr+oxsiz*(y<<1)+(x0<<1);
			sw = ssw+xsiz*y+x0;
			for(x=x0;x<x1;x++)
			{
					//ÚÄÄÄÂÄÄÄÂÄÄÄÂÄÄÄ¿
					//³npt³z1 ³z1c³dum³
					//³ b ³ g ³ r ³ i ³
					//³ b ³ g ³ r ³ i ³
					//³npt³z1 ³z1c³z0 ³
					//³ b ³ g ³ r ³ i ³
					//ÀÄÄÄÁÄÄÄÁÄÄÄÁÄÄÄÙ
				v[0] = sr[      0];
				v[1] = sr[      1];
				v[2] = sr[oysiz  ];
				v[3] = sr[oysiz+1];
				for(i=3;i>=0;i--)
				{
					curz[i] = curzn[i][0] = (long)v[i][1];
					curzn[i][1] = ((long)v[i][2])+1;

					tv = v[i];
					while (1)
					{
						oz = (long)tv[1];
						for(z=oz;z<=((long)tv[2]);z++)
						{
							nz = (z>>1);
							mixc[nz][mixn[nz]++] = *(long *)(&tv[((z-oz)<<2)+4]);
						}
						z = (z-oz) - (((long)tv[0])-1);
						if (!tv[0]) break;
						tv += (((long)tv[0])<<2);
						oz = (long)tv[3];
						for(;z<0;z++)
						{
							nz = ((z+oz)>>1);
							mixc[nz][mixn[nz]++] = *(long *)(&tv[z<<2]);
						}
					}
				}
				cstat = 0; oldn = 0; n = 4; tbuf[3] = 0; z = 0x80000000;
				while (1)
				{
					oz = z;

						//z,besti = min,argmin(curz[0],curz[1],curz[2],curz[3])
					besti = (((unsigned long)(curz[1]-curz[    0]))>>31);
						 i = (((unsigned long)(curz[3]-curz[    2]))>>31)+2;
					besti +=(((( signed long)(curz[i]-curz[besti]))>>31)&(i-besti));
					z = curz[besti]; if (z >= MAXZDIM) break;

					if ((!cstat) && ((z>>1) >= ((oz+1)>>1)))
					{
						if (oz >= 0)
						{
							tbuf[oldn] = ((n-oldn)>>2);
							tbuf[oldn+2]--;
							tbuf[n+3] = ((oz+1)>>1);
							oldn = n; n += 4;
						}
						tbuf[oldn] = 0;
						tbuf[oldn+1] = tbuf[oldn+2] = (z>>1); cz = -1;
					}
					if (cstat&0x1111)
					{
						if (((((long)tbuf[oldn+2])<<1)+1 >= oz) && (cz < 0))
						{
							while ((((long)tbuf[oldn+2])<<1) < z)
							{
								zz = (long)tbuf[oldn+2];

								*(long *)&tbuf[n] = mixc[zz][rand()%mixn[zz]];
								mixn[zz] = 0;

								tbuf[oldn+2]++; n += 4;
							}
						}
						else
						{
							if (cz < 0) cz = (oz>>1);
							else if ((cz<<1)+1 < oz)
							{
									//Insert fake slab
								tbuf[oldn] = ((n-oldn)>>2);
								tbuf[oldn+2]--;
								tbuf[n] = 0;
								tbuf[n+1] = tbuf[n+2] = tbuf[n+3] = cz;
								oldn = n; n += 4;
								cz = (oz>>1);
							}
							while ((cz<<1) < z)
							{
								*(long *)&tbuf[n] = mixc[cz][rand()%mixn[cz]];
								mixn[cz] = 0;

								cz++; n += 4;
							}
						}
					}

					i = (besti<<2);
					cstat = (((1<<i)+cstat)&0x3333); //--33--22--11--00
					switch ((cstat>>i)&3)
					{
						case 0: curz[besti] = curzn[besti][0]; break;
						case 1: curz[besti] = curzn[besti][1]; break;
						case 2:
							if (!(v[besti][0])) { curz[besti] = MAXZDIM; }
							else
							{
								tv = v[besti]; i = (((long)tv[2])-((long)tv[1])+1)-(((long)tv[0])-1);
								tv += (((long)tv[0])<<2);
								curz[besti] = ((long)(tv[3])) + i;
								curzn[besti][3] = (long)(tv[3]);
								curzn[besti][0] = (long)(tv[1]);
								curzn[besti][1] = ((long)tv[2])+1;
								v[besti] = tv;
							}
							break;
						case 3: curz[besti] = curzn[besti][3]; break;
						//default: __assume(0); //tells MSVC default can't be reached
					}
				}
				tbuf[oldn+2]--;
				if (cz >= 0)
				{
					tbuf[oldn] = ((n-oldn)>>2);
					tbuf[n] = 0;
					tbuf[n+1] = tbuf[n+3] = cz;
					tbuf[n+2] = cz-1;
					n += 4;
				}

					//De-allocate column (x,y) if it exists
				if (sw[0]) voxdealloc(sw[0]);

					//Allocate & copy to new column (x,y)
				sw[0] = voxalloc(n);
				copybuf((void *)tbuf,(void *)sw[0],n>>2);
				sw++; sr += 2;
			}
			sr += ysiz*2;
		}
		ssr = ssw; ssw += xsiz*ysiz;
		mipnum++; if (mipnum > gmipnum) gmipnum = mipnum;
	}

		//Remove extra mips (bbox must be 0,0,VSID,VSID to get inside this)
	while ((xsiz > 1) && (ysiz > 1) && (zsiz > 1) && (mipnum < gmipnum))
	{
		xsiz >>= 1; ysiz >>= 1; zsiz >>= 1;
		for(i=xsiz*ysiz;i>0;i--)
		{
			if (ssw[0]) voxdealloc(ssw[0]); //De-allocate column if it exists
			ssw++;
		}
		gmipnum--;
	}

	_asm emms

#if 0 //TEMP HACK!!!
	{
	FILE *fil;
	dpoint3d dp;
	if (!(fil = fopen("temp512.vxl","wb"))) return;
	i = 0x09072000; fwrite(&i,4,1,fil);  //Version
	i = (VSID>>1); fwrite(&i,4,1,fil);
	i = (VSID>>1); fwrite(&i,4,1,fil);
	dp.x = (double)i*.5; dp.y = (double)i*.5; dp.z = (double)i*.5;
	fwrite(&dp,24,1,fil);
	dp.x = 1.0; dp.y = 0.0; dp.z = 0.0; fwrite(&dp,24,1,fil);
	dp.x = 0.0; dp.y = 0.0; dp.z = 1.0; fwrite(&dp,24,1,fil);
	dp.x = 0.0; dp.y =-1.0; dp.z = 0.0; fwrite(&dp,24,1,fil);
	for(i=0;i<(VSID>>1)*(VSID>>1);i++)
		fwrite((void *)sptr[i+VSID*VSID],slng(sptr[i+VSID*VSID]),1,fil);
	fclose(fil);
	}
	gmipnum = 1;
#endif

}

//------------------------- SXL parsing code begins --------------------------

static char *sxlbuf = 0;
static long sxlparspos, sxlparslen;

long loadsxl (const char *sxlnam, char **vxlnam, char **skynam, char **globst)
{
	long j, k, m, n;

	printf("loadsxl %s\n", sxlnam);

		//NOTE: MUST buffer file because insertsprite uses kz file code :/
	if (!kzopen(sxlnam)) return(0);
	sxlparslen = kzfilelength();
	if (sxlbuf) { free(sxlbuf); sxlbuf = 0; }
	if (!(sxlbuf = (char *)malloc(sxlparslen))) return(0);
	kzread(sxlbuf,sxlparslen);
	kzclose();

	j = n = 0;

		//parse vxlnam
	(*vxlnam) = &sxlbuf[j];
	while ((sxlbuf[j]!=13)&&(sxlbuf[j]!=10) && (j < sxlparslen)) j++; sxlbuf[j++] = 0;
	while (((sxlbuf[j]==13)||(sxlbuf[j]==10)) && (j < sxlparslen)) j++;

		//parse skynam
	(*skynam) = &sxlbuf[j];
	while ((sxlbuf[j]!=13)&&(sxlbuf[j]!=10) && (j < sxlparslen)) j++; sxlbuf[j++] = 0;
	while (((sxlbuf[j]==13)||(sxlbuf[j]==10)) && (j < sxlparslen)) j++;

		//parse globst
	m = n = j; (*globst) = &sxlbuf[n];
	while (((sxlbuf[j] == ' ') || (sxlbuf[j] == 9)) && (j < sxlparslen))
	{
		j++;
		while ((sxlbuf[j]!=13) && (sxlbuf[j]!=10) && (j < sxlparslen)) sxlbuf[n++] = sxlbuf[j++];
		sxlbuf[n++] = 13; j++;
		while (((sxlbuf[j]==13) || (sxlbuf[j]==10)) && (j < sxlparslen)) j++;
	}
	if (n > m) sxlbuf[n-1] = 0; else (*globst) = &nullst;

		//Count number of sprites in .SXL file (helpful for GAME)
	sxlparspos = j;
	return(1);
}

#ifdef __cplusplus
extern "C" {
#endif

extern void *caddasm;
#define cadd4 ((point4d *)&caddasm)
extern void *ztabasm;
#define ztab4 ((point4d *)&ztabasm)
extern short qsum0[4], qsum1[4], qbplbpp[4];
extern float scisdist;
extern __int64 kv6colmul[256], kv6coladd[256];

char ptfaces16[43][8] =
{
	0, 0, 0,  0,  0, 0, 0,0,  4, 0,32,96, 64, 0,32,0,  4,16,80,112,48, 16,80,0,  0,0,0,0,0,0,0,0,
	4,64,96,112, 80,64,96,0,  6, 0,32,96,112,80,64,0,  6,16,80, 64,96,112,48,0,  0,0,0,0,0,0,0,0,
	4, 0,16, 48, 32, 0,16,0,  6, 0,16,48, 32,96,64,0,  6, 0,16, 80,112,48,32,0,  0,0,0,0,0,0,0,0,
	0, 0, 0,  0,  0, 0, 0,0,  0, 0, 0, 0,  0, 0, 0,0,  0, 0, 0,  0,  0, 0, 0,0,  0,0,0,0,0,0,0,0,
	4, 0,64, 80, 16, 0,64,0,  6, 0,32,96, 64,80,16,0,  6, 0,64, 80,112,48,16,0,  0,0,0,0,0,0,0,0,
	6, 0,64, 96,112,80,16,0,  6, 0,32,96,112,80,16,0,  6, 0,64, 96,112,48,16,0,  0,0,0,0,0,0,0,0,
	6, 0,64, 80, 16,48,32,0,  6,16,48,32, 96,64,80,0,  6, 0,64, 80,112,48,32,0,  0,0,0,0,0,0,0,0,
	0, 0, 0,  0,  0, 0, 0,0,  0, 0, 0, 0,  0, 0, 0,0,  0, 0, 0,  0,  0, 0, 0,0,  0,0,0,0,0,0,0,0,
	4,32,48,112, 96,32,48,0,  6, 0,32,48,112,96,64,0,  6,16,80,112, 96,32,48,0,  0,0,0,0,0,0,0,0,
	6,32,48,112, 80,64,96,0,  6, 0,32,48,112,80,64,0,  6,16,80, 64, 96,32,48,0,  0,0,0,0,0,0,0,0,
	6, 0,16, 48,112,96,32,0,  6, 0,16,48,112,96,64,0,  6, 0,16, 80,112,96,32,0,
};

#ifdef __cplusplus
}
#endif

static __int64 all32767 = 0x7fff7fff7fff7fff;

#endif

	//Updates mip-mapping
typedef struct { long x0, y0, z0, x1, y1, z1, csgdel; } bboxtyp;
#define BBOXSIZ 256
static bboxtyp bbox[BBOXSIZ];
static long bboxnum = 0;
void updatevxl ()
{
	long i;

	for(i=bboxnum-1;i>=0;i--)
	{
		if (vx5.vxlmipuse > 1)
			genmipvxl(bbox[i].x0,bbox[i].y0,bbox[i].x1,bbox[i].y1);
	}
	bboxnum = 0;
}

void updatebbox (long x0, long y0, long z0, long x1, long y1, long z1, long csgdel)
{
	long i;

	if ((x0 >= x1) || (y0 >= y1) || (z0 >= z1)) return;
	for(i=bboxnum-1;i>=0;i--)
	{
		if ((x0 >= bbox[i].x1) || (bbox[i].x0 >= x1)) continue;
		if ((y0 >= bbox[i].y1) || (bbox[i].y0 >= y1)) continue;
		if ((z0 >= bbox[i].z1) || (bbox[i].z0 >= z1)) continue;
		if (bbox[i].x0 < x0) x0 = bbox[i].x0;
		if (bbox[i].y0 < y0) y0 = bbox[i].y0;
		if (bbox[i].z0 < z0) z0 = bbox[i].z0;
		if (bbox[i].x1 > x1) x1 = bbox[i].x1;
		if (bbox[i].y1 > y1) y1 = bbox[i].y1;
		if (bbox[i].z1 > z1) z1 = bbox[i].z1;
		csgdel |= bbox[i].csgdel;
		bboxnum--; bbox[i] = bbox[bboxnum];
	}
	bbox[bboxnum].x0 = x0; bbox[bboxnum].x1 = x1;
	bbox[bboxnum].y0 = y0; bbox[bboxnum].y1 = y1;
	bbox[bboxnum].z0 = z0; bbox[bboxnum].z1 = z1;
	bbox[bboxnum].csgdel = csgdel; bboxnum++;
	if (bboxnum >= BBOXSIZ) updatevxl();
}

//----------------------------------------------------------------------------

void voxsetframebuffer (long p, long b, long x, long y)
{
	long i;

	frameplace = p;
	if (x > MAXXDIM) x = MAXXDIM; //This sucks, but it crashes without it
	if (y > MAXYDIM) y = MAXYDIM;

		//Set global variables used by kv6draw's PIII asm (drawboundcube)
	qsum1[3] = qsum1[1] = 0x7fff-y; qsum1[2] = qsum1[0] = 0x7fff-x;

	if ((b != ylookup[1]) || (x != xres) || (y != yres))
	{
		bytesperline = b; xres = x; yres = y; xres4 = (xres<<2);
		ylookup[0] = 0; for(i=0;i<yres;i++) ylookup[i+1] = ylookup[i]+bytesperline;
		//gihx = gihz = (float)xres*.5f; gihy = (float)yres*.5f; //BAD!!!
#if (USEZBUFFER == 1)
		if ((ylookup[yres]+256 > zbuffersiz) || (!zbuffermem))  //Increase Z buffer size if too small
		{
			if (zbuffermem) { free(zbuffermem); zbuffermem = 0; }
			zbuffersiz = ylookup[yres]+256;
			if (!(zbuffermem = (long *)malloc(zbuffersiz))) evilquit("voxsetframebuffer: allocation too big");
		}
#endif
	}
#if (USEZBUFFER == 1)
		//zbuffer aligns its memory to the same pixel boundaries as the screen!
		//WARNING: Pentium 4's L2 cache has severe slowdowns when 65536-64 <= (zbufoff&65535) < 64
	zbufoff = (((((long)zbuffermem)-frameplace-128)+255)&~255)+128;
#endif
	uurend = &uurendmem[((frameplace&4)^(((long)uurendmem)&4))>>2];

	if (vx5.fogcol >= 0)
	{
		fogcol = (((__int64)(vx5.fogcol&0xff0000))<<16) +
					(((__int64)(vx5.fogcol&0x00ff00))<< 8) +
					(((__int64)(vx5.fogcol&0x0000ff))    );

		if (vx5.maxscandist > 2047) vx5.maxscandist = 2047;
		if ((vx5.maxscandist != ofogdist) && (vx5.maxscandist > 0))
		{
			ofogdist = vx5.maxscandist;

			//foglut[?>>20] = min(?*32767/vx5.maxscandist,32767)
#if 0
			long j, k, l;
			j = 0; l = 0x7fffffff/vx5.maxscandist;
			for(i=0;i<2048;i++)
			{
				k = (j>>16); j += l;
				if (k < 0) break;
				foglut[i] = (((__int64)k)<<32)+(((__int64)k)<<16)+((__int64)k);
			}
			while (i < 2048) foglut[i++] = all32767;
#else
			i = 0x7fffffff/vx5.maxscandist;
			_asm
			{
				xor eax, eax
				mov ecx, -2048*8
				mov edx, i
fogbeg:     movd mm0, eax
				add eax, edx
				jo short fogend
				pshufw mm0, mm0, 0x55
				movq foglut[ecx+2048*8], mm0
				add ecx, 8
				js short fogbeg
				jmp short fogend2
fogend:     movq mm0, all32767
fogbeg2:    movq foglut[ecx+2048*8], mm0
				add ecx, 8
				js short fogbeg2
fogend2:    emms
			}
#endif
		}
	} else ofogdist = -1;
}

void uninitvoxlap ()
{
	if (sxlbuf) { free(sxlbuf); sxlbuf = 0; }

	if (vbuf) { free(vbuf); vbuf = 0; }
	if (vbit) { free(vbit); vbit = 0; }

	if (skylng) { free((void *)skylng); skylng = 0; }
	if (skylat) { free((void *)skylat); skylat = 0; }
	if (skypic) { free((void *)skypic); skypic = skyoff = 0; }

	if (vx5.pic) { free(vx5.pic); vx5.pic = 0; }
#if (USEZBUFFER == 1)
	if (zbuffermem) { free(zbuffermem); zbuffermem = 0; }
#endif
	if (radarmem) { free(radarmem); radarmem = 0; radar = 0; }
}

long initvoxlap ()
{
	__int64 q;
	long i, j, k, z, zz;
	float f, ff;

	v5_asm_dep_unlock();

		//CPU Must have: FPU,RDTSC,CMOV,MMX,MMX+
	if ((cputype&((1<<0)|(1<<4)|(1<<15)|(1<<22)|(1<<23))) !=
					 ((1<<0)|(1<<4)|(1<<15)|(1<<22)|(1<<23))) return(-1);
		//CPU UNSUPPORTED!
	if ((!(cputype&(1<<25))) && //SSE
		(!((cputype&((1<<30)|(1<<31))) == ((1<<30)|(1<<31))))) //3DNow!+
		return(-1);
	//if (cputype&(1<<25)) fixsse(); //SSE

	  //WARNING: xres&yres are local to VOXLAP5.C so don't rely on them here!
	if (!(radarmem = (long *)malloc(max((((MAXXDIM*MAXYDIM*27)>>1)+7)&~7,(VSID+4)*3*SCPITCH*4+8))))
		return(-1);
	radar = (long *)((((long)radarmem)+7)&~7);

	for(i=0;i<32;i++) { xbsflor[i] = (-1<<i); xbsceil[i] = ~xbsflor[i]; }

#if (ESTNORMRAD == 2)
		//LUT for ESTNORM
	fsqrecip[0] = 0.f; fsqrecip[1] = 1.f;
	fsqrecip[2] = (float)(1.f/sqrt(2.f)); fsqrecip[3] = (float)1.f/sqrt(3.f);
	for(z=4,i=3;z<sizeof(fsqrecip)/sizeof(fsqrecip[0]);z+=6) //fsqrecip[z] = 1/sqrt(z);
	{
		fsqrecip[z+0] = fsqrecip[(z+0)>>1]*fsqrecip[2];
		fsqrecip[z+2] = fsqrecip[(z+2)>>1]*fsqrecip[2];
		fsqrecip[z+4] = fsqrecip[(z+4)>>1]*fsqrecip[2];
		fsqrecip[z+5] = fsqrecip[i]*fsqrecip[3]; i += 2;

		f = (fsqrecip[z+0]+fsqrecip[z+2])*.5f;
		if (z <= 22) f = (1.5f-(.5f*((float)(z+1))) * f*f)*f;
		fsqrecip[z+1] = (1.5f-(.5f*((float)(z+1))) * f*f)*f;

		f = (fsqrecip[z+2]+fsqrecip[z+4])*.5f;
		if (z <= 22) f = (1.5f-(.5f*((float)(z+3))) * f*f)*f;
		fsqrecip[z+3] = (1.5f-(.5f*((float)(z+3))) * f*f)*f;
	}
#endif

		//Lookup table to save 1 divide for gline()
	for(i=1;i<CMPRECIPSIZ;i++) cmprecip[i] = CMPPREC/(float)i;

		//Flashscan equal-angle compare table
	for(i=0;i<(1<<LOGFLASHVANG)*8;i++)
	{
		if (!(i&((1<<LOGFLASHVANG)-1)))
			j = (gfclookup[i>>LOGFLASHVANG]<<4)+8 - (1<<LOGFLASHVANG)*64;
		gfc[i].y = j; j += 64*2;
		ftol(sqrt((1<<(LOGFLASHVANG<<1))*64.f*64.f-gfc[i].y*gfc[i].y),&gfc[i].x);
	}

		//Init norm flash variables:
	ff = (float)GSIZ*.5f; // /(1);
	for(z=1;z<(GSIZ>>1);z++)
	{
		ffxptr = &ffx[(z+1)*z-1];
		f = ff; ff = (float)GSIZ*.5f/((float)z+1);
		for(zz=-z;zz<=z;zz++)
		{
			if (zz <= 0) i = (long)(((float)zz-.5f)*f); else i = (long)(((float)zz-.5f)*ff);
			if (zz >= 0) j = (long)(((float)zz+.5f)*f); else j = (long)(((float)zz+.5f)*ff);
			ffxptr[zz].x = (unsigned short)max(i+(GSIZ>>1),0);
			ffxptr[zz].y = (unsigned short)min(j+(GSIZ>>1),GSIZ);
		}
	}
	for(i=0;i<=25*5;i+=5) xbsbuf[i] = 0x00000000ffffffff;
	for(z=0;z<32;z++) { p2c[z] = (1<<z); p2m[z] = p2c[z]-1; }

	ucossininit();

	memset(mixn,0,sizeof(mixn));

	vx5.anginc = 1; //Higher=faster (1:full,2:half)
	vx5.sideshademode = 0; setsideshades(0,0,0,0,0,0);
	vx5.mipscandist = 128;
	vx5.maxscandist = 256; //must be <= 2047
	vx5.colfunc = curcolfunc; //This prevents omission bugs from crashing voxlap5
	vx5.curcol = 0x80804c33;
	vx5.currad = 8;
	vx5.curhei = 0;
	vx5.curpow = 2.0;
	vx5.amount = 0x70707;
	vx5.pic = 0;
	vx5.cliphitnum = 0;
	vx5.flstnum = 0;
	vx5.kv6col = 0x808080;
	vx5.vxlmipuse = 1;
	vx5.fogcol = -1;
	vx5.fallcheck = 0;

	gmipnum = 0;

	return(0);
}

// ----- GAME.C code begins

	//Player position variables:
#define CLIPRAD 5
dpoint3d ipos, istr, ihei, ifor, ivel;

long lockanginc = 0;

	//Mouse button state global variables:
long obstatus = 0, bstatus = 0;

	//Timer global variables:
double odtotclk, dtotclk;
float fsynctics;
long totclk;

long initmap ()
{
	char cursxlnam[MAX_PATH+1] = "vxl/default.sxl";
	char *vxlnam, *skynam, *userst;

	if (!loadsxl(cursxlnam,&vxlnam,&skynam,&userst)) {
		printf("failed to load '%s'\n", cursxlnam);
		exit(1);
	}

	if (!loadvxl(vxlnam,&ipos,&istr,&ihei,&ifor)) {
		printf("failed to load '%s'\n", vxlnam);
		exit(1);
	}

	vx5.vxlmipuse = 9;
	vx5.mipscandist = 192;
	vx5.fallcheck = 1;

	updatevxl();

	vx5.maxscandist = (long)(VSID*1.42);

	ivel.x = ivel.y = ivel.z = 0;

	return 0;
}

long initapp (long argc, char **argv)
{
	prognam = "\"Ken-VOXLAP\" test game";
	xres = 640; yres = 480; colbits = 32; fullscreen = 0;

	if (initvoxlap() < 0) return(-1);

	setsideshades(0,4,1,3,2,2);

		//AthlonXP 2000+: SSE:26.76ms, 3DN:27.34ms, SSE2:28.93ms
	extern long cputype;
	switch(2) // SSE2
	{
		case 0: cputype &= ~((1<<25)|(1<<26)); cputype |= ((1<<30)|(1<<31)); break;
		case 1: cputype |= (1<<25); cputype &= ~(1<<26); cputype &= ~((1<<30)|(1<<31)); break;
		case 2: cputype |= ((1<<25)|(1<<26)); cputype &= ~((1<<30)|(1<<31)); break;
		default:;
	}

	if (initmap() < 0) return(-1);

		//Init klock
	readklock(&dtotclk);
	totclk = (long)(dtotclk*1000.0);

	fsynctics = 1.f;

	return 0;
}

void uninitapp ()
{
	long i;
	uninitvoxlap();
}

void doframe ()
{
	dpoint3d dp, dp2, dp3, dpos;
	point3d fp, fp2, fp3, fpos;
	lpoint3d lp, lp2;
	double d;
	float f, fmousx, fmousy;
	long i, j, k, l, m, *hind, hdir;
	char tempbuf[260];

	if (!startdirectdraw(&i,&j,&k,&l)) goto skipalldraw;

	voxsetframebuffer(i,j,k,l);
	setcamera(&ipos,&istr,&ihei,&ifor,xres*.5,yres*.5,xres*.5);
	setears3d(ipos.x,ipos.y,ipos.z,ifor.x,ifor.y,ifor.z,ihei.x,ihei.y,ihei.z);
	opticast();

	stopdirectdraw();
	nextpage();
skipalldraw:;

		//Read keyboard, mouse, and timer
	readkeyboard();
	obstatus = bstatus; readmouse(&fmousx,&fmousy,0,&bstatus);
	odtotclk = dtotclk; readklock(&dtotclk);
	totclk = (long)(dtotclk*1000.0); fsynctics = (float)(dtotclk-odtotclk);

		//Rotate player's view
	dp.x = istr.z*.1; dp.y = fmousy*.008; dp.z = fmousx*.008;
	dorthorotate(dp.x,dp.y,dp.z,&istr,&ihei,&ifor);

		//Draw less when turning fast to increase the frame rate
	if (!lockanginc)
	{
		i = 1; if (xres*yres >= 640*350) i++;
		f = dp.x*dp.x + dp.y*dp.y + dp.z*dp.z;
		if (f >= .01*.01)
		{
			i++;
			if (f >= .08*.08)
			{
				i++; if (f >= .64*.64) i++;
			}
		}
		vx5.anginc = (float)i;
	}

		//Move player and perform simple physics (gravity,momentum,friction)
	f = fsynctics*12.0;
	if (keystatus[0x1e]) { ivel.x -= istr.x*f; ivel.y -= istr.y*f; ivel.z -= istr.z*f; } // A
	if (keystatus[0x20]) { ivel.x += istr.x*f; ivel.y += istr.y*f; ivel.z += istr.z*f; } // D
	if (keystatus[0x11]) { ivel.x += ifor.x*f; ivel.y += ifor.y*f; ivel.z += ifor.z*f; } // W
	if (keystatus[0x1f]) { ivel.x -= ifor.x*f; ivel.y -= ifor.y*f; ivel.z -= ifor.z*f; } // S
	if (keystatus[0x12]) { ivel.x -= ihei.x*f; ivel.y -= ihei.y*f; ivel.z -= ihei.z*f; } // E
	if (keystatus[0x10]) { ivel.x += ihei.x*f; ivel.y += ihei.y*f; ivel.z += ihei.z*f; } // Q
	//ivel.z += fsynctics*2.0; //Gravity (used to be *4.0)
	f = fsynctics*64.0;
	dp.x = ivel.x*f;
	dp.y = ivel.y*f;
	dp.z = ivel.z*f;
	dp2 = ipos; clipmove(&ipos,&dp,CLIPRAD);
	if (f != 0)
	{
		f = .9/f; //Friction
		ivel.x = (ipos.x-dp2.x)*f;
		ivel.y = (ipos.y-dp2.y)*f;
		ivel.z = (ipos.z-dp2.z)*f;
	}

	if (vx5.clipmaxcr < CLIPRAD*.9)
	{
		if (vx5.clipmaxcr >= CLIPRAD*.1) //Try to push player to safety
		{
			if (vx5.cliphitnum > 0)
			{
				switch (vx5.cliphitnum)
				{
					case 1: dpos.x = dp2.x-vx5.cliphit[0].x;
							  dpos.y = dp2.y-vx5.cliphit[0].y;
							  dpos.z = dp2.z-vx5.cliphit[0].z;
							  break;
					case 2: dpos.x = dp2.x-(vx5.cliphit[0].x+vx5.cliphit[1].x)*.5;
							  dpos.y = dp2.y-(vx5.cliphit[0].y+vx5.cliphit[1].y)*.5;
							  dpos.z = dp2.z-(vx5.cliphit[0].z+vx5.cliphit[1].z)*.5;
							  break;
					case 3:
						 dp.x = (vx5.cliphit[1].x-vx5.cliphit[0].x);
						 dp.y = (vx5.cliphit[1].y-vx5.cliphit[0].y);
						 dp.z = (vx5.cliphit[1].z-vx5.cliphit[0].z);
						dp3.x = (vx5.cliphit[2].x-vx5.cliphit[0].x);
						dp3.y = (vx5.cliphit[2].y-vx5.cliphit[0].y);
						dp3.z = (vx5.cliphit[2].z-vx5.cliphit[0].z);
						dpos.x = dp.y*dp3.z - dp.z*dp3.y;
						dpos.y = dp.z*dp3.x - dp.x*dp3.z;
						dpos.z = dp.x*dp3.y - dp.y*dp3.x;

						  //Ugly hack making sure cross product points in right direction
						if ((dp2.x*3-(vx5.cliphit[0].x+vx5.cliphit[1].x+vx5.cliphit[2].x))*dpos.x +
							 (dp2.y*3-(vx5.cliphit[0].y+vx5.cliphit[1].y+vx5.cliphit[2].y))*dpos.y +
							 (dp2.z*3-(vx5.cliphit[0].z+vx5.cliphit[1].z+vx5.cliphit[2].z))*dpos.z < 0)
							{ dpos.x = -dpos.x; dpos.y = -dpos.y; dpos.z = -dpos.z; }

						break;
				}
				d = dpos.x*dpos.x + dpos.y*dpos.y + dpos.z*dpos.z;
				if (d > 0)
				{
					d = 1.0/sqrt(d);
					ivel.x += dpos.x*d;
					ivel.y += dpos.y*d;
					ivel.z += dpos.z*d;
				}
			}
		}
		else //Out of room... squish player!
		{
			puts("You got squished!");
			exit(0);
		}
	}
	f = ivel.x*ivel.x + ivel.y*ivel.y + ivel.z*ivel.z; //Limit maximum velocity
	if (f > 8.0*8.0) { f = 8.0/sqrt(f); ivel.x *= f; ivel.y *= f; ivel.z *= f; }

	updatevxl();
}

/// ------- KPLIB code begins

typedef struct
{
	FILE *fil;   //0:no file open, !=0:open file (either stand-alone or zip)
	int comptyp; //0:raw data (can be ZIP or stand-alone), 8:PKZIP LZ77 *flate
	int seek0;   //0:stand-alone file, !=0: start of zip compressed stream data
	int compleng;//Global variable for compression FIFO
	int comptell;//Global variable for compression FIFO
	int leng;    //Uncompressed file size (bytes)
	int pos;     //Current uncompressed relative file position (0<=pos<=leng)
	int endpos;  //Temp global variable for kzread
	int jmpplc;  //Store place where decompression paused
	int i;       //For stand-alone/ZIP comptyp#0, this is like "uncomptell"
					  //For ZIP comptyp#8&btype==0 "<64K store", this saves i state
	int bfinal;  //LZ77 decompression state (for later calls)
} kzfilestate;
static kzfilestate kzfs;

INT_PTR kzopen (const char *filnam)
{
	kzfs.fil = fopen(filnam,"rb");

	if (!kzfs.fil)
	{
		printf("file '%s' not found\n", filnam);
		exit(1);
	}

	kzfs.comptyp = 0;
	kzfs.seek0 = 0;
	kzfs.leng = filelength(_fileno(kzfs.fil));
	kzfs.pos = 0;
	kzfs.i = 0;

	return (INT_PTR)kzfs.fil;
}

	//returns number of bytes copied
int kzread (void *buffer, int leng)
{
	int i, j, k, bfinal, btype, hlit, hdist;

	if ((!kzfs.fil) || (leng <= 0)) return(0);

	if (kzfs.pos != kzfs.i) //Seek only when position changes
		fseek(kzfs.fil,kzfs.seek0+kzfs.pos,SEEK_SET);
	i = min(kzfs.leng-kzfs.pos,leng);
	fread(buffer,i,1,kzfs.fil);
	kzfs.i += i; //kzfs.i is a local copy of ftell(kzfs.fil);

	i = kzfs.pos;
	kzfs.pos += leng; if (kzfs.pos > kzfs.leng) kzfs.pos = kzfs.leng;
	return(kzfs.pos-i);
}

int kzfilelength ()
{
	if (!kzfs.fil) return(0);
	return(kzfs.leng);
}

int kztell ()
{
	if (!kzfs.fil) return(-1);
	return(kzfs.pos);
}

void kzclose ()
{
	if (kzfs.fil) { fclose(kzfs.fil); kzfs.fil = 0; }
}
