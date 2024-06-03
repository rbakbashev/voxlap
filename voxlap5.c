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

typedef struct
{
	long parent;      //index to parent sprite (-1=none)
	point3d p[2];     //"velcro" point of each object
	point3d v[2];     //axis of rotation for each object
	short vmin, vmax; //min value / max value
	char htype, filler[7];
} hingetype;

typedef struct { long tim, frm; } seqtyp;

typedef struct
{
	long numspr, numhin, numfrm, seqnum;
	long namoff;
	kv6data *basekv6;      //Points to original unconnected KV6 (maybe helpful?)
	struct vx5sprite *spr; //[numspr]
	hingetype *hinge;      //[numhin]
	long *hingesort;       //[numhin]
	short *frmval;         //[numfrm][numhin]
	seqtyp *seq;           //[seqnum]
} kfatype;

	//Notice that I aligned each point3d on a 16-byte boundary. This will be
	//   helpful when I get around to implementing SSE instructions someday...
typedef struct vx5sprite
{
	point3d p; //position in VXL coordinates
	long flags; //flags bit 0:0=use normal shading, 1=disable normal shading
					//flags bit 1:0=points to kv6data, 1=points to kfatype
					//flags bit 2:0=normal, 1=invisible sprite
	static union { point3d s, x; }; //kv6data.xsiz direction in VXL coordinates
	static union
	{
		kv6data *voxnum; //pointer to KV6 voxel data (bit 1 of flags = 0)
		kfatype *kfaptr; //pointer to KFA animation  (bit 1 of flags = 1)
	};
	static union { point3d h, y; }; //kv6data.ysiz direction in VXL coordinates
	long kfatim;        //time (in milliseconds) of KFA animation
	static union { point3d f, z; }; //kv6data.zsiz direction in VXL coordinates
	long okfatim;       //make vx5sprite exactly 64 bytes :)
} vx5sprite;

	//Falling voxels shared data: (flst = float list)
#define FLPIECES 256 //Max # of separate falling pieces
typedef struct //(68 bytes)
{
	lpoint3d chk; //a solid point on piece (x,y,pointer) (don't touch!)
	long i0, i1; //indices to start&end of slab list (don't touch!)
	long x0, y0, z0, x1, y1, z1; //bounding box, written by startfalls
	long mass; //mass of piece, written by startfalls (1 unit per voxel)
	point3d centroid; //centroid of piece, written by startfalls

		//userval is set to -1 when a new piece is spawned. Voxlap does not
		//read or write these values after that point. You should use these to
		//play an initial sound and track velocity
	long userval, userval2;
} flstboxtype;

	//Lighting variables: (used by updatelighting)
#define MAXLIGHTS 256
typedef struct { point3d p; float r2, sc; } lightsrctype;

	//Used by setspans/meltspans. Ordered this way to allow sorting as longs!
typedef struct { char z1, z0, x, y; } vspans;

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
	flstboxtype flstcnt[FLPIECES];

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

		//Lighting variables: (used by updatelighting)
	long lightmode; //0 (default), 1:simple lighting, 2:lightsrc lighting
	lightsrctype lightsrc[MAXLIGHTS]; //(?,?,?),128*128,262144
	long numlights;

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
extern void drawpoint2d (long, long, long);
extern void drawpoint3d (float, float, float, long);
extern void drawline2d (float, float, float, float, long);
extern void drawline3d (float, float, float, float, float, float, long);
extern long project2d (float, float, float, float *, float *, float *);
extern void drawpicinquad (long, long, long, long, long, long, long, long, float, float, float, float, float, float, float, float);
extern void drawpolyquad (long, long, long, long, float, float, float, float, float, float, float, float, float, float, float, float, float, float, float, float, float, float);

	//Sprite related functions:
extern void freekv6 (kv6data *kv6);
extern char *getkfilname (long);
extern long meltsphere (vx5sprite *, lpoint3d *, long);
extern long meltspans (vx5sprite *, vspans *, long, lpoint3d *);

	//Physics helper functions:
extern void dorthorotate (double, double, double, dpoint3d *, dpoint3d *, dpoint3d *);
extern long cansee (point3d *, point3d *, lpoint3d *);
extern double findmaxcr (double, double, double, double);
extern void clipmove (dpoint3d *, dpoint3d *, double);
extern void estnorm (long, long, long, point3d *);

	//VXL reading functions (fast!):
extern long isvoxelsolid (long, long, long);
extern long anyvoxelsolid (long, long, long, long);
extern long anyvoxelempty (long, long, long, long);
extern long getfloorz (long, long, long);
extern long getcube (long, long, long);

	//VXL writing functions (optimized & bug-free):
extern void setcube (long, long, long, long);
extern void setsphere (lpoint3d *, long, long);
extern void setellipsoid (lpoint3d *, lpoint3d *, long, long, long);
extern void setcylinder (lpoint3d *, lpoint3d *, long, long, long);
extern void setrect (lpoint3d *, lpoint3d *, long);
extern void settri (point3d *, point3d *, point3d *, long);
extern void setsector (point3d *, long *, long, float, long, long);
extern void setspans (vspans *, long, lpoint3d *, long);
extern void setheightmap (const unsigned char *, long, long, long, long, long, long, long);
extern void setkv6 (vx5sprite *, long);

	//VXL MISC functions:
extern void updatebbox (long, long, long, long, long, long, long);
extern void updatevxl ();
extern void genmipvxl (long, long, long, long);
extern void updatelighting (long, long, long, long, long, long);

	//Falling voxels functions:
extern void checkfloatinbox (long, long, long, long, long, long);
extern void startfalls ();
extern void dofall (long);
extern long meltfall (vx5sprite *, long, long);
extern void finishfalls ();

	//Procedural texture functions:
extern long curcolfunc (lpoint3d *);
extern long floorcolfunc (lpoint3d *);
extern long jitcolfunc (lpoint3d *);
extern long manycolfunc (lpoint3d *);
extern long sphcolfunc (lpoint3d *);
extern long woodcolfunc (lpoint3d *);
extern long pngcolfunc (lpoint3d *);
extern long kv6colfunc (lpoint3d *);

	//Editing backup/restore functions
extern void voxbackup (long, long, long, long, long);
extern void voxdontrestore ();
extern void voxrestore ();
extern void voxredraw ();

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

#define SETSPHMAXRAD 256
static double logint[SETSPHMAXRAD];
static float tempfloatbuf[SETSPHMAXRAD];
static long factr[SETSPHMAXRAD][2];

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

long lightvox (long i)
{
	long r, g, b;

	b = ((unsigned long)i>>24);
	r = min((((i>>16)&255)*b)>>7,255);
	g = min((((i>>8 )&255)*b)>>7,255);
	b = min((((i    )&255)*b)>>7,255);
	return((r<<16)+(g<<8)+b);
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

#if 0

	//Point: (x,y), line segment: (px,py)-(px+vx,py+vy)
	//Returns 1 if point is closer than sqrt(cr2) to line
long dist2linept2d (double x, double y, double px, double py, double vx, double vy, double cr2)
{
	double f, g;
	x -= px; y -= py; f = x*vx + y*vy; if (f <= 0) return(x*x + y*y <= cr2);
	g = vx*vx + vy*vy; if (f >= g) { x -= vx; y -= vy; return(x*x + y*y <= cr2); }
	x = x*g-vx*f; y = y*g-vy*f; return(x*x + y*y <= cr2*g*g);
}

static char clipbuf[MAXZDIM+16]; //(8 extra on each side)
long sphtraceo (double px, double py, double pz,    //start pt
					double vx, double vy, double vz,    //move vector
					double *nx, double *ny, double *nz, //new pt after collision
					double *fx, double *fy, double *fz, //pt that caused collision
					double cr, double acr)
{
	double t, u, ex, ey, ez, Za, Zb, Zc, thresh2;
	double vxyz, vyz, vxz, vxy, rvxyz, rvyz, rvxz, rvxy, rvx, rvy, rvz, cr2;
	long i, i0, i1, x, y, z, xx, yy, zz, v, vv, ix, y0, y1, z0, z1;
	char *vp;

	t = 1;
	(*nx) = px + vx;
	(*ny) = py + vy;
	(*nz) = pz + vz;

	z0 = max((long)(min(pz,*nz)-cr)-2,-1);
	z1 = min((long)(max(pz,*nz)+cr)+2,MAXZDIM);

	thresh2 = cr+1.7321+1; thresh2 *= thresh2;

	vyz = vz*vz; vxz = vx*vx; vxy = vy*vy; cr2 = cr*cr;
	vyz += vxy; vxy += vxz; vxyz = vyz + vxz; vxz += vz*vz;
	rvx = 1.0 / vx; rvy = 1.0 / vy; rvz = 1.0 / vz;
	rvyz = 1.0 / vyz; rvxz = 1.0 / vxz; rvxy = 1.0 / vxy;
	rvxyz = 1.0 / vxyz;

		//Algorithm fails (stops short) if cr < 2 :(
	i0 = i1 = 0; ix = x = (long)px; y = y0 = y1 = (long)py;
	while (1)
	{
		for(z=z0;z<=z1;z++) clipbuf[z+8] = 0;
		i = 16;
		for(yy=y;yy<y+2;yy++)
			for(xx=x;xx<x+2;xx++,i<<=1)
			{
				z = z0;
				if ((unsigned long)(xx|yy) < VSID)
				{
					vp = sptr[yy*VSID+xx];
					while (1)
					{
						if (vp[1] > z) z = vp[1];
						if (!vp[0]) break;
						vp += vp[0]*4;
						zz = vp[3]; if (zz > z1) zz = z1;
						while (z < zz) clipbuf[(z++)+8] |= i;
					}
				}
				while (z <= z1) clipbuf[(z++)+8] |= i;
			}

		xx = x+1; yy = y+1; v = clipbuf[z0+8];
		for(z=z0;z<z1;z++)
		{
			zz = z+1; v = (v>>4)|clipbuf[zz+8];
			if ((!v) || (v == 255)) continue;

//---------------Check 1(8) corners of cube (sphere intersection)-------------

			//if (((v-1)^v) >= v)  //True if v is: {1,2,4,8,16,32,64,128}
			if (!(v&(v-1)))      //Same as above, but {0,1,2,4,...} (v's never 0)
			{
				ex = xx-px; ey = yy-py; ez = zz-pz;
				Zb = ex*vx + ey*vy + ez*vz;
				Zc = ex*ex + ey*ey + ez*ez - cr2;
				u = Zb*Zb - vxyz*Zc;
				if ((((long *)&u)[1] | ((long *)&Zb)[1]) >= 0)
				//if ((u >= 0) && (Zb >= 0))
				{
						//   //Proposed compare optimization:
						//f = Zb*Zb-u; g = vxyz*t; h = (Zb*2-g)*g;
						//if ((unsigned __int64 *)&f < (unsigned __int64 *)&h)
					u = (Zb - sqrt(u)) * rvxyz;
					if ((u >= 0) && (u < t))
					{
						*fx = xx; *fy = yy; *fz = zz; t = u;
						*nx = vx*u + px; *ny = vy*u + py; *nz = vz*u + pz;
					}
				}
			}

//---------------Check 3(12) edges of cube (cylinder intersection)-----------

			vv = v&0x55; if (((vv-1)^vv) >= vv)  //True if (v&0x55)={1,4,16,64}
			{
				ey = yy-py; ez = zz-pz;
				Zb = ey*vy + ez*vz;
				Zc = ey*ey + ez*ez - cr2;
				u = Zb*Zb - vyz*Zc;
				if ((((long *)&u)[1] | ((long *)&Zb)[1]) >= 0)
				//if ((u >= 0) && (Zb >= 0))
				{
					u = (Zb - sqrt(u)) * rvyz;
					if ((u >= 0) && (u < t))
					{
						ex = vx*u + px;
						if ((ex >= x) && (ex <= xx))
						{
							*fx = ex; *fy = yy; *fz = zz; t = u;
							*nx = ex; *ny = vy*u + py; *nz = vz*u + pz;
						}
					}
				}
			}
			vv = v&0x33; if (((vv-1)^vv) >= vv) //True if (v&0x33)={1,2,16,32}
			{
				ex = xx-px; ez = zz-pz;
				Zb = ex*vx + ez*vz;
				Zc = ex*ex + ez*ez - cr2;
				u = Zb*Zb - vxz*Zc;
				if ((((long *)&u)[1] | ((long *)&Zb)[1]) >= 0)
				//if ((u >= 0) && (Zb >= 0))
				{
					u = (Zb - sqrt(u)) * rvxz;
					if ((u >= 0) && (u < t))
					{
						ey = vy*u + py;
						if ((ey >= y) && (ey <= yy))
						{
							*fx = xx; *fy = ey; *fz = zz; t = u;
							*nx = vx*u + px; *ny = ey; *nz = vz*u + pz;
						}
					}
				}
			}
			vv = v&0x0f; if (((vv-1)^vv) >= vv) //True if (v&0x0f)={1,2,4,8}
			{
				ex = xx-px; ey = yy-py;
				Zb = ex*vx + ey*vy;
				Zc = ex*ex + ey*ey - cr2;
				u = Zb*Zb - vxy*Zc;
				if ((((long *)&u)[1] | ((long *)&Zb)[1]) >= 0)
				//if ((u >= 0) && (Zb >= 0))
				{
					u = (Zb - sqrt(u)) * rvxy;
					if ((u >= 0) && (u < t))
					{
						ez = vz*u + pz;
						if ((ez >= z) && (ez <= zz))
						{
							*fx = xx; *fy = yy; *fz = ez; t = u;
							*nx = vx*u + px; *ny = vy*u + py; *nz = ez;
						}
					}
				}
			}

//---------------Check 3(6) faces of cube (plane intersection)---------------

			if (vx)
			{
				switch(v&0x03)
				{
					case 0x01: ex = xx+cr; if ((vx > 0) || (px < ex)) goto skipfacex; break;
					case 0x02: ex = xx-cr; if ((vx < 0) || (px > ex)) goto skipfacex; break;
					default: goto skipfacex;
				}
				u = (ex - px) * rvx;
				if ((u >= 0) && (u < t))
				{
					ey = vy*u + py;
					ez = vz*u + pz;
					if ((ey >= y) && (ey <= yy) && (ez >= z) && (ez <= zz))
					{
						*fx = xx; *fy = ey; *fz = ez; t = u;
						*nx = ex; *ny = ey; *nz = ez;
					}
				}
			}
skipfacex:;
			if (vy)
			{
				switch(v&0x05)
				{
					case 0x01: ey = yy+cr; if ((vy > 0) || (py < ey)) goto skipfacey; break;
					case 0x04: ey = yy-cr; if ((vy < 0) || (py > ey)) goto skipfacey; break;
					default: goto skipfacey;
				}
				u = (ey - py) * rvy;
				if ((u >= 0) && (u < t))
				{
					ex = vx*u + px;
					ez = vz*u + pz;
					if ((ex >= x) && (ex <= xx) && (ez >= z) && (ez <= zz))
					{
						*fx = ex; *fy = yy; *fz = ez; t = u;
						*nx = ex; *ny = ey; *nz = ez;
					}
				}
			}
skipfacey:;
			if (vz)
			{
				switch(v&0x11)
				{
					case 0x01: ez = zz+cr; if ((vz > 0) || (pz < ez)) goto skipfacez; break;
					case 0x10: ez = zz-cr; if ((vz < 0) || (pz > ez)) goto skipfacez; break;
					default: goto skipfacez;
				}
				u = (ez - pz) * rvz;
				if ((u >= 0) && (u < t))
				{
					ex = vx*u + px;
					ey = vy*u + py;
					if ((ex >= x) && (ex <= xx) && (ey >= y) && (ey <= yy))
					{
						*fx = ex; *fy = ey; *fz = zz; t = u;
						*nx = ex; *ny = ey; *nz = ez;
					}
				}
			}
skipfacez:;
		}

		if ((x <= ix) && (x > 0) && (dist2linept2d(x-1,y,px,py,vx,vy,thresh2)))
			{ clipit[i1].x = x-1; clipit[i1].y = y; i1 = ((i1+1)&(MAXCLIPIT-1)); }
		if ((x >= ix) && (x < VSID-1) && (dist2linept2d(x+1,y,px,py,vx,vy,thresh2)))
			{ clipit[i1].x = x+1; clipit[i1].y = y; i1 = ((i1+1)&(MAXCLIPIT-1)); }
		if ((y <= y0) && (y > 0) && (dist2linept2d(x,y-1,px,py,vx,vy,thresh2)))
			{ clipit[i1].x = x; clipit[i1].y = y-1; i1 = ((i1+1)&(MAXCLIPIT-1)); y0--; }
		if ((y >= y1) && (y < VSID-1) && (dist2linept2d(x,y+1,px,py,vx,vy,thresh2)))
			{ clipit[i1].x = x; clipit[i1].y = y+1; i1 = ((i1+1)&(MAXCLIPIT-1)); y1++; }
		if (i0 == i1) break;
		x = clipit[i0].x; y = clipit[i0].y; i0 = ((i0+1)&(MAXCLIPIT-1));
	}

	if ((*nx) < acr) (*nx) = acr;
	if ((*ny) < acr) (*ny) = acr;
	if ((*nx) > VSID-acr) (*nx) = VSID-acr;
	if ((*ny) > VSID-acr) (*ny) = VSID-acr;
	if ((*nz) > MAXZDIM-1-acr) (*nz) = MAXZDIM-1-acr;
	if ((*nz) < MAXZDIM-2048) (*nz) = MAXZDIM-2048;

	return (t == 1);
}

#endif

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

long cansee (point3d *p0, point3d *p1, lpoint3d *hit)
{
	lpoint3d a, c, d, p, i;
	point3d f, g;
	long cnt;

	ftol(p0->x-.5,&a.x); ftol(p0->y-.5,&a.y); ftol(p0->z-.5,&a.z);
	if (isvoxelsolid(a.x,a.y,a.z)) { hit->x = a.x; hit->y = a.y; hit->z = a.z; return(0); }
	ftol(p1->x-.5,&c.x); ftol(p1->y-.5,&c.y); ftol(p1->z-.5,&c.z);
	cnt = 0;

		  if (c.x <  a.x) { d.x = -1; f.x = p0->x-a.x;   g.x = (p0->x-p1->x)*1024; cnt += a.x-c.x; }
	else if (c.x != a.x) { d.x =  1; f.x = a.x+1-p0->x; g.x = (p1->x-p0->x)*1024; cnt += c.x-a.x; }
	else f.x = g.x = 0;
		  if (c.y <  a.y) { d.y = -1; f.y = p0->y-a.y;   g.y = (p0->y-p1->y)*1024; cnt += a.y-c.y; }
	else if (c.y != a.y) { d.y =  1; f.y = a.y+1-p0->y; g.y = (p1->y-p0->y)*1024; cnt += c.y-a.y; }
	else f.y = g.y = 0;
		  if (c.z <  a.z) { d.z = -1; f.z = p0->z-a.z;   g.z = (p0->z-p1->z)*1024; cnt += a.z-c.z; }
	else if (c.z != a.z) { d.z =  1; f.z = a.z+1-p0->z; g.z = (p1->z-p0->z)*1024; cnt += c.z-a.z; }
	else f.z = g.z = 0;

	ftol(f.x*g.z - f.z*g.x,&p.x); ftol(g.x,&i.x);
	ftol(f.y*g.z - f.z*g.y,&p.y); ftol(g.y,&i.y);
	ftol(f.y*g.x - f.x*g.y,&p.z); ftol(g.z,&i.z);

		//NOTE: GIGO! This can happen if p0,p1 (cansee input) is NaN, Inf, etc...
	if ((unsigned long)cnt > (VSID+VSID+2048)*2) cnt = (VSID+VSID+2048)*2;
	while (cnt > 0)
	{
		if (((p.x|p.y) >= 0) && (a.z != c.z)) { a.z += d.z; p.x -= i.x; p.y -= i.y; }
		else if ((p.z >= 0) && (a.x != c.x))  { a.x += d.x; p.x += i.z; p.z -= i.y; }
		else                                  { a.y += d.y; p.y += i.z; p.z += i.x; }
		if (isvoxelsolid(a.x,a.y,a.z)) break;
		cnt--;
	}
	hit->x = a.x; hit->y = a.y; hit->z = a.z; return(!cnt);
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

	printf("loadvxl %s\n", lodfilnam);

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

void expandrle (long x, long y, long *uind)
{
	long i;
	char *v;

	if ((x|y)&(~(VSID-1))) { uind[0] = 0; uind[1] = MAXZDIM; return; }

	v = sptr[y*VSID+x]; uind[0] = v[1]; i = 2;
	while (v[0])
	{
		v += v[0]*4; if (v[3] >= v[1]) continue;
		uind[i-1] = v[3]; uind[i] = v[1]; i += 2;
	}
	uind[i-1] = MAXZDIM;
}

	//Inputs:  n0[<=MAXZDIM]: rle buffer of column to compress
	//         n1-4[<=MAXZDIM]: neighboring rle buffers
	//         top,bot,top,bot,... (ends when bot == MAXZDIM)
	//         px,py: takes color from original column (before modification)
	//            If originally unexposed, calls vx5.colfunc(.)
	//Outputs: cbuf[MAXCSIZ]: compressed output buffer
	//Returns: n: length of compressed buffer (in bytes)
long compilerle (long *n0, long *n1, long *n2, long *n3, long *n4, char *cbuf, long px, long py)
{
	long i, ia, ze, zend, onext, dacnt, n, *ic;
	lpoint3d p;
	char *v;

	p.x = px; p.y = py;

		//Generate pointers to color slabs in this format:
		//   0:z0,  1:z1,  2:(pointer to z0's color)-z0
	v = sptr[py*VSID+px]; ic = tbuf2;
	while (1)
	{
		ia = v[1]; p.z = v[2];
		ic[0] = ia; ic[1] = p.z+1; ic[2] = ((long)v)-(ia<<2)+4; ic += 3;
		i = v[0]; if (!i) break;
		v += i*4; ze = v[3];
		ic[0] = ze+p.z-ia-i+2; ic[1] = ze; ic[2] = ((long)v)-(ze<<2); ic += 3;
	}
	ic[0] = MAXZDIM; ic[1] = MAXZDIM;

	p.z = n0[0]; cbuf[1] = n0[0];
	ze = n0[1]; cbuf[2] = ze-1;
	cbuf[3] = 0;
	i = onext = 0; ic = tbuf2; ia = 15; n = 4;
	if (ze != MAXZDIM) zend = ze-1; else zend = -1;
	while (1)
	{
		dacnt = 0;
		while (1)
		{
			do
			{
				while (p.z >= ic[1]) ic += 3;
				if (p.z >= ic[0]) *(long *)&cbuf[n] = *(long *)(ic[2]+(p.z<<2));
								 else *(long *)&cbuf[n] = vx5.colfunc(&p);
				n += 4; p.z++; if (p.z >= ze) goto rlendit2;
				while (p.z >= n1[0]) { n1++; ia ^= 1; }
				while (p.z >= n2[0]) { n2++; ia ^= 2; }
				while (p.z >= n3[0]) { n3++; ia ^= 4; }
				while (p.z >= n4[0]) { n4++; ia ^= 8; }
			} while ((ia) || (p.z == zend));

			if (!dacnt) { cbuf[onext+2] = p.z-1; dacnt = 1; }
			else
			{
				cbuf[onext] = ((n-onext)>>2); onext = n;
				cbuf[n+1] = p.z; cbuf[n+2] = p.z-1; cbuf[n+3] = p.z; n += 4;
			}

			if ((n1[0] < n2[0]) && (n1[0] < n3[0]) && (n1[0] < n4[0]))
				{ if (n1[0] >= ze) { p.z = ze-1; } else { p.z = *n1++; ia ^= 1; } }
			else if ((n2[0] < n3[0]) && (n2[0] < n4[0]))
				{ if (n2[0] >= ze) { p.z = ze-1; } else { p.z = *n2++; ia ^= 2; } }
			else if (n3[0] < n4[0])
				{ if (n3[0] >= ze) { p.z = ze-1; } else { p.z = *n3++; ia ^= 4; } }
			else
				{ if (n4[0] >= ze) { p.z = ze-1; } else { p.z = *n4++; ia ^= 8; } }

			if (p.z == MAXZDIM-1) goto rlenditall;
		}
rlendit2:;
		if (ze >= MAXZDIM) break;

		i += 2;
		cbuf[onext] = ((n-onext)>>2); onext = n;
		p.z = n0[i]; cbuf[n+1] = n0[i]; cbuf[n+3] = ze;
		ze = n0[i+1]; cbuf[n+2] = ze-1;
		n += 4;
	}
rlenditall:;
	cbuf[onext] = 0;
	return(n);
}

	//Delete everything on b2() in y0<=y<y1
void delslab (long *b2, long y0, long y1)
{
	long i, j, z;

	if (y1 >= MAXZDIM) y1 = MAXZDIM-1;
	if ((y0 >= y1) || (!b2)) return;
	for(z=0;y0>=b2[z+1];z+=2);
	if (y0 > b2[z])
	{
		if (y1 < b2[z+1])
		{
			for(i=z;b2[i+1]<MAXZDIM;i+=2);
			while (i > z) { b2[i+3] = b2[i+1]; b2[i+2] = b2[i]; i -= 2; }
			b2[z+3] = b2[z+1]; b2[z+1] = y0; b2[z+2] = y1; return;
		}
		b2[z+1] = y0; z += 2;
	}
	if (y1 >= b2[z+1])
	{
		for(i=z+2;y1>=b2[i+1];i+=2);
		j = z-i; b2[z] = b2[i]; b2[z+1] = b2[i+1];
		while (b2[i+1] < MAXZDIM)
			{ i += 2; b2[i+j] = b2[i]; b2[i+j+1] = b2[i+1]; }
	}
	if (y1 > b2[z]) b2[z] = y1;
}

	//Insert everything on b2() in y0<=y<y1
void insslab (long *b2, long y0, long y1)
{
	long i, j, z;

	if ((y0 >= y1) || (!b2)) return;
	for(z=0;y0>b2[z+1];z+=2);
	if (y1 < b2[z])
	{
		for(i=z;b2[i+1]<MAXZDIM;i+=2);
		do { b2[i+3] = b2[i+1]; b2[i+2] = b2[i]; i -= 2; } while (i >= z);
		b2[z+1] = y1; b2[z] = y0; return;
	}
	if (y0 < b2[z]) b2[z] = y0;
	if ((y1 >= b2[z+2]) && (b2[z+1] < MAXZDIM))
	{
		for(i=z+2;(y1 >= b2[i+2]) && (b2[i+1] < MAXZDIM);i+=2);
		j = z-i; b2[z+1] = b2[i+1];
		while (b2[i+1] < MAXZDIM)
			{ i += 2; b2[i+j] = b2[i]; b2[i+j+1] = b2[i+1]; }
	}
	if (y1 > b2[z+1]) b2[z+1] = y1;
}

//------------------------ SETCOLUMN CODE BEGINS ----------------------------

static long scx0, scx1, scox0, scox1, scoox0, scoox1;
static long scex0, scex1, sceox0, sceox1, scoy = 0x80000000, *scoym3;

void scumline ()
{
	long i, j, k, x, y, x0, x1, *mptr, *uptr;
	char *v;

	x0 = min(scox0-1,min(scx0,scoox0)); scoox0 = scox0; scox0 = scx0;
	x1 = max(scox1+1,max(scx1,scoox1)); scoox1 = scox1; scox1 = scx1;

	uptr = &scoym3[SCPITCH]; if (uptr == &radar[SCPITCH*9]) uptr = &radar[SCPITCH*6];
	mptr = &uptr[SCPITCH];   if (mptr == &radar[SCPITCH*9]) mptr = &radar[SCPITCH*6];

	if ((x1 < sceox0) || (x0 > sceox1))
	{
		for(x=x0;x<=x1;x++) expandstack(x,scoy-2,&uptr[x*SCPITCH*3]);
	}
	else
	{
		for(x=x0;x<sceox0;x++) expandstack(x,scoy-2,&uptr[x*SCPITCH*3]);
		for(x=x1;x>sceox1;x--) expandstack(x,scoy-2,&uptr[x*SCPITCH*3]);
	}

	if ((scex1|x1) >= 0)
	{
		for(x=x1+2;x<scex0;x++) expandstack(x,scoy-1,&mptr[x*SCPITCH*3]);
		for(x=x0-2;x>scex1;x--) expandstack(x,scoy-1,&mptr[x*SCPITCH*3]);
	}
	if ((x1+1 < scex0) || (x0-1 > scex1))
	{
		for(x=x0-1;x<=x1+1;x++) expandstack(x,scoy-1,&mptr[x*SCPITCH*3]);
	}
	else
	{
		for(x=x0-1;x<scex0;x++) expandstack(x,scoy-1,&mptr[x*SCPITCH*3]);
		for(x=x1+1;x>scex1;x--) expandstack(x,scoy-1,&mptr[x*SCPITCH*3]);
	}
	sceox0 = min(x0-1,scex0);
	sceox1 = max(x1+1,scex1);

	if ((x1 < scx0) || (x0 > scx1))
	{
		for(x=x0;x<=x1;x++) expandstack(x,scoy,&scoym3[x*SCPITCH*3]);
	}
	else
	{
		for(x=x0;x<scx0;x++) expandstack(x,scoy,&scoym3[x*SCPITCH*3]);
		for(x=x1;x>scx1;x--) expandstack(x,scoy,&scoym3[x*SCPITCH*3]);
	}
	scex0 = x0;
	scex1 = x1;

	y = scoy-1; if (y&(~(VSID-1))) return;
	if (x0 < 0) x0 = 0;
	if (x1 >= VSID) x1 = VSID-1;
	i = y*VSID+x0; k = x0*SCPITCH*3;
	for(x=x0;x<=x1;x++,i++,k+=SCPITCH*3)
	{
		v = sptr[i]; vx5.globalmass += v[1];
		while (v[0]) { v += v[0]*4; vx5.globalmass += v[1]-v[3]; }

			//De-allocate column (x,y)
		voxdealloc(sptr[i]);

		j = compilestack(&mptr[k],&mptr[k-SCPITCH*3],&mptr[k+SCPITCH*3],&uptr[k],&scoym3[k],
							  tbuf,x,y);

			//Allocate & copy to new column (x,y)
		sptr[i] = v = voxalloc(j); copybuf((void *)tbuf,(void *)v,j>>2);

		vx5.globalmass -= v[1];
		while (v[0]) { v += v[0]*4; vx5.globalmass += v[3]-v[1]; }
	}
}

	//x: x on voxel map
	//y: y on voxel map
	//z0: highest z on column
	//z1: lowest z(+1) on column
	//nbuf: buffer of color data from nbuf[z0] to nbuf[z1-1];
	//           -3: don't modify voxel
	//           -2: solid voxel (unexposed): to be calculated in compilestack
	//           -1: write air voxel
	//   0-16777215: write solid voxel (exposed)
void scum (long x, long y, long z0, long z1, long *nbuf)
{
	long z, *mptr;

	if ((x|y)&(~(VSID-1))) return;

	if (y != scoy)
	{
		if (scoy >= 0)
		{
			scumline();
			while (scoy < y-1)
			{
				scx0 = 0x7fffffff; scx1 = 0x80000000;
				scoy++; scoym3 += SCPITCH; if (scoym3 == &radar[SCPITCH*9]) scoym3 = &radar[SCPITCH*6];
				scumline();
			}
			scoy++; scoym3 += SCPITCH; if (scoym3 == &radar[SCPITCH*9]) scoym3 = &radar[SCPITCH*6];
		}
		else
		{
			scoox0 = scox0 = 0x7fffffff;
			sceox0 = scex0 = x+1;
			sceox1 = scex1 = x;
			scoy = y; scoym3 = &radar[SCPITCH*6];
		}
		scx0 = x;
	}
	else
	{
		while (scx1 < x-1) { scx1++; expandstack(scx1,y,&scoym3[scx1*SCPITCH*3]); }
	}

	mptr = &scoym3[x*SCPITCH*3]; scx1 = x; expandstack(x,y,mptr);

		//Modify column (x,y):
	if (nbuf[MAXZDIM-1] == -1) nbuf[MAXZDIM-1] = -2; //Bottom voxel must be solid
	for(z=z0;z<z1;z++)
		if (nbuf[z] != -3) mptr[z] = nbuf[z];
}

void scumfinish ()
{
	long i;

	if (scoy == 0x80000000) return;
	for(i=2;i;i--)
	{
		scumline(); scx0 = 0x7fffffff; scx1 = 0x80000000;
		scoy++; scoym3 += SCPITCH; if (scoym3 == &radar[SCPITCH*9]) scoym3 = &radar[SCPITCH*6];
	}
	scumline(); scoy = 0x80000000;
}

	//Example of how to use this code:
	//vx5.colfunc = curcolfunc; //0<x0<x1<VSID, 0<y0<y1<VSID, 0<z0<z1<256,
	//clearbuf((void *)&templongbuf[z0],z1-z0,-1); //Ex: set all voxels to air
	//for(y=y0;y<y1;y++) //MUST iterate x&y in this order, but can skip around
	//   for(x=x0;x<x1;x++)
	//      if (rand()&8) scum(x,y,z0,z1,templongbuf));
	//scumfinish(); //MUST call this when done!

void scum2line ()
{
	long i, j, k, x, y, x0, x1, *mptr, *uptr;
	char *v;

	x0 = min(scox0-1,min(scx0,scoox0)); scoox0 = scox0; scox0 = scx0;
	x1 = max(scox1+1,max(scx1,scoox1)); scoox1 = scox1; scox1 = scx1;

	uptr = &scoym3[SCPITCH]; if (uptr == &radar[SCPITCH*9]) uptr = &radar[SCPITCH*6];
	mptr = &uptr[SCPITCH];   if (mptr == &radar[SCPITCH*9]) mptr = &radar[SCPITCH*6];

	if ((x1 < sceox0) || (x0 > sceox1))
	{
		for(x=x0;x<=x1;x++) expandrle(x,scoy-2,&uptr[x*SCPITCH*3]);
	}
	else
	{
		for(x=x0;x<sceox0;x++) expandrle(x,scoy-2,&uptr[x*SCPITCH*3]);
		for(x=x1;x>sceox1;x--) expandrle(x,scoy-2,&uptr[x*SCPITCH*3]);
	}

	if ((scex1|x1) >= 0)
	{
		for(x=x1+2;x<scex0;x++) expandrle(x,scoy-1,&mptr[x*SCPITCH*3]);
		for(x=x0-2;x>scex1;x--) expandrle(x,scoy-1,&mptr[x*SCPITCH*3]);
	}
	if ((x1+1 < scex0) || (x0-1 > scex1))
	{
		for(x=x0-1;x<=x1+1;x++) expandrle(x,scoy-1,&mptr[x*SCPITCH*3]);
	}
	else
	{
		for(x=x0-1;x<scex0;x++) expandrle(x,scoy-1,&mptr[x*SCPITCH*3]);
		for(x=x1+1;x>scex1;x--) expandrle(x,scoy-1,&mptr[x*SCPITCH*3]);
	}
	sceox0 = min(x0-1,scex0);
	sceox1 = max(x1+1,scex1);

	if ((x1 < scx0) || (x0 > scx1))
	{
		for(x=x0;x<=x1;x++) expandrle(x,scoy,&scoym3[x*SCPITCH*3]);
	}
	else
	{
		for(x=x0;x<scx0;x++) expandrle(x,scoy,&scoym3[x*SCPITCH*3]);
		for(x=x1;x>scx1;x--) expandrle(x,scoy,&scoym3[x*SCPITCH*3]);
	}
	scex0 = x0;
	scex1 = x1;

	y = scoy-1; if (y&(~(VSID-1))) return;
	if (x0 < 0) x0 = 0;
	if (x1 >= VSID) x1 = VSID-1;
	i = y*VSID+x0; k = x0*SCPITCH*3;
	for(x=x0;x<=x1;x++,i++,k+=SCPITCH*3)
	{
		j = compilerle(&mptr[k],&mptr[k-SCPITCH*3],&mptr[k+SCPITCH*3],&uptr[k],&scoym3[k],
							tbuf,x,y);

		v = sptr[i]; vx5.globalmass += v[1];
		while (v[0]) { v += v[0]*4; vx5.globalmass += v[1]-v[3]; }

			//De-allocate column (x,y)  Note: Must be AFTER compilerle!
		voxdealloc(sptr[i]);

			//Allocate & copy to new column (x,y)
		sptr[i] = v = voxalloc(j); copybuf((void *)tbuf,(void *)v,j>>2);

		vx5.globalmass -= v[1];
		while (v[0]) { v += v[0]*4; vx5.globalmass += v[3]-v[1]; }
	}
}

	//x: x on voxel map
	//y: y on voxel map
	//Returns pointer to rle column (x,y)
long *scum2 (long x, long y)
{
	long *mptr;

	if ((x|y)&(~(VSID-1))) return(0);

	if (y != scoy)
	{
		if (scoy >= 0)
		{
			scum2line();
			while (scoy < y-1)
			{
				scx0 = 0x7fffffff; scx1 = 0x80000000;
				scoy++; scoym3 += SCPITCH; if (scoym3 == &radar[SCPITCH*9]) scoym3 = &radar[SCPITCH*6];
				scum2line();
			}
			scoy++; scoym3 += SCPITCH; if (scoym3 == &radar[SCPITCH*9]) scoym3 = &radar[SCPITCH*6];
		}
		else
		{
			scoox0 = scox0 = 0x7fffffff;
			sceox0 = scex0 = x+1;
			sceox1 = scex1 = x;
			scoy = y; scoym3 = &radar[SCPITCH*6];
		}
		scx0 = x;
	}
	else
	{
		while (scx1 < x-1) { scx1++; expandrle(scx1,y,&scoym3[scx1*SCPITCH*3]); }
	}

	mptr = &scoym3[x*SCPITCH*3]; scx1 = x; expandrle(x,y,mptr);
	return(mptr);
}

void scum2finish ()
{
	long i;

	if (scoy == 0x80000000) return;
	for(i=2;i;i--)
	{
		scum2line(); scx0 = 0x7fffffff; scx1 = 0x80000000;
		scoy++; scoym3 += SCPITCH; if (scoym3 == &radar[SCPITCH*9]) scoym3 = &radar[SCPITCH*6];
	}
	scum2line(); scoy = 0x80000000;
}

//------------------- editing backup / restore begins ------------------------

void voxdontrestore ()
{
	long i;

	if (backedup == 1)
	{
		for(i=(bacx1-bacx0)*(bacy1-bacy0)-1;i>=0;i--) voxdealloc(bacsptr[i]);
	}
	backedup = -1;
}

void voxrestore ()
{
	long i, j, x, y;
	char *v, *daptr;

	if (backedup == 1)
	{
		i = 0;
		for(y=bacy0;y<bacy1;y++)
		{
			j = y*VSID;
			for(x=bacx0;x<bacx1;x++)
			{
				v = sptr[j+x]; vx5.globalmass += v[1];
				while (v[0]) { v += v[0]*4; vx5.globalmass += v[1]-v[3]; }

				voxdealloc(sptr[j+x]);
				sptr[j+x] = bacsptr[i]; i++;

				v = sptr[j+x]; vx5.globalmass -= v[1];
				while (v[0]) { v += v[0]*4; vx5.globalmass += v[3]-v[1]; }
			}
		}
		if (vx5.vxlmipuse > 1) genmipvxl(bacx0,bacy0,bacx1,bacy1);
	}
	else if (backedup == 2)
	{
		daptr = (char *)bacsptr;
		for(y=bacy0;y<bacy1;y++)
		{
			j = y*VSID;
			for(x=bacx0;x<bacx1;x++)
			{
				for(v=sptr[j+x];v[0];v+=v[0]*4)
					for(i=1;i<v[0];i++) v[(i<<2)+3] = *daptr++;
				for(i=1;i<=v[2]-v[1]+1;i++) v[(i<<2)+3] = *daptr++;
			}
		}
		if (vx5.vxlmipuse > 1) genmipvxl(bacx0,bacy0,bacx1,bacy1);
	}
	backedup = -1;
}

void voxbackup (long x0, long y0, long x1, long y1, long tag)
{
	long i, j, n, x, y;
	char *v, *daptr;

	voxdontrestore();

	x0 = max(x0-2,0); y0 = max(y0-2,0);
	x1 = min(x1+2,VSID); y1 = min(y1+2,VSID);
	if ((x1-x0)*(y1-y0) > 262144) return;

	bacx0 = x0; bacy0 = y0; bacx1 = x1; bacy1 = y1; backtag = tag;

	if (tag&0x10000)
	{
		backedup = 1;
		i = 0;
		for(y=bacy0;y<bacy1;y++)
		{
			j = y*VSID;
			for(x=bacx0;x<bacx1;x++)
			{
				bacsptr[i] = v = sptr[j+x]; i++;
				n = slng(v); sptr[j+x] = voxalloc(n);

				copybuf((void *)v,(void *)sptr[j+x],n>>2);
			}
		}
	}
	else if (tag&0x20000)
	{
		backedup = 2;
			//WARNING!!! Right now, function will crash if saving more than
			//   1<<20 colors :( This needs to be addressed!!!
		daptr = (char *)bacsptr;
		for(y=bacy0;y<bacy1;y++)
		{
			j = y*VSID;
			for(x=bacx0;x<bacx1;x++)
			{
				for(v=sptr[j+x];v[0];v+=v[0]*4)
					for(i=1;i<v[0];i++) *daptr++ = v[(i<<2)+3];
				for(i=1;i<=v[2]-v[1]+1;i++) *daptr++ = v[(i<<2)+3];
			}
		}
	}
	else backedup = 0;
}

//-------------------- editing backup / restore ends -------------------------

	//WARNING! Make sure to set vx5.colfunc before calling this function!
	//This function is here for simplicity only - it is NOT optimal.
	//
	//   -1: set air
	//   -2: use vx5.colfunc
void setcube (long px, long py, long pz, long col)
{
	long bakcol, (*bakcolfunc)(lpoint3d *), *lptr;

	vx5.minx = px; vx5.maxx = px+1;
	vx5.miny = py; vx5.maxy = py+1;
	vx5.minz = pz; vx5.maxz = pz+1;
	if ((unsigned long)pz >= MAXZDIM) return;
	if ((unsigned long)col >= (unsigned long)0xfffffffe) //-1 or -2
	{
		lptr = scum2(px,py);
		if (col == -1) delslab(lptr,pz,pz+1); else insslab(lptr,pz,pz+1);
		scum2finish();
		updatebbox(vx5.minx,vx5.miny,vx5.minz,vx5.maxx,vx5.maxy,vx5.maxz,col);
		return;
	}

	bakcol = getcube(px,py,pz);
	if (bakcol == 1) return; //Unexposed solid
	if (bakcol != 0) //Not 0 (air)
		*(long *)bakcol = col;
	else
	{
		bakcolfunc = vx5.colfunc; bakcol = vx5.curcol;
		vx5.colfunc = curcolfunc; vx5.curcol = col;
		insslab(scum2(px,py),pz,pz+1); scum2finish();
		vx5.colfunc = bakcolfunc; vx5.curcol = bakcol;
	}
	updatebbox(vx5.minx,vx5.miny,vx5.minz,vx5.maxx,vx5.maxy,vx5.maxz,0);
}

//-------------------------- SETCOLUMN CODE ENDS ----------------------------

static __int64 qmulmip[8] =
{
	0x7fff7fff7fff7fff,0x4000400040004000,0x2aaa2aaa2aaa2aaa,0x2000200020002000,
	0x1999199919991999,0x1555155515551555,0x1249124912491249,0x1000100010001000
};
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

void setsphere (lpoint3d *hit, long hitrad, long dacol)
{
	void (*modslab)(long *, long, long);
	long i, x, y, xs, ys, zs, xe, ye, ze, sq;
	float f, ff;

	xs = max(hit->x-hitrad,0); xe = min(hit->x+hitrad,VSID-1);
	ys = max(hit->y-hitrad,0); ye = min(hit->y+hitrad,VSID-1);
	zs = max(hit->z-hitrad,0); ze = min(hit->z+hitrad,MAXZDIM-1);
	vx5.minx = xs; vx5.maxx = xe+1;
	vx5.miny = ys; vx5.maxy = ye+1;
	vx5.minz = zs; vx5.maxz = ze+1;
	if ((xs > xe) || (ys > ye) || (zs > ze)) return;

	if (vx5.colfunc == sphcolfunc)
	{
		vx5.cen = hit->x+hit->y+hit->z;
		vx5.daf = 1.f/(hitrad*sqrt(3.f));
	}

	if (hitrad >= SETSPHMAXRAD-1) hitrad = SETSPHMAXRAD-2;
	if (dacol == -1) modslab = delslab; else modslab = insslab;

	tempfloatbuf[0] = 0.0f;
#if 0
		//Totally unoptimized
	for(i=1;i<=hitrad;i++) tempfloatbuf[i] = pow((float)i,vx5.curpow);
#else
	tempfloatbuf[1] = 1.0f;
	for(i=2;i<=hitrad;i++)
	{
		if (!factr[i][0]) tempfloatbuf[i] = exp(logint[i]*vx5.curpow);
		else tempfloatbuf[i] = tempfloatbuf[factr[i][0]]*tempfloatbuf[factr[i][1]];
	}
#endif
	*(long *)&tempfloatbuf[hitrad+1] = 0x7f7fffff; //3.4028235e38f; //Highest float

	sq = 0; //pow(fabs(x-hit->x),vx5.curpow) + "y + "z < pow(vx5.currad,vx5.curpow)
	for(y=ys;y<=ye;y++)
	{
		ff = tempfloatbuf[hitrad]-tempfloatbuf[labs(y-hit->y)];
		if (*(long *)&ff <= 0) continue;
		for(x=xs;x<=xe;x++)
		{
			f = ff-tempfloatbuf[labs(x-hit->x)]; if (*(long *)&f <= 0) continue;
			while (*(long *)&tempfloatbuf[sq] <  *(long *)&f) sq++;
			while (*(long *)&tempfloatbuf[sq] >= *(long *)&f) sq--;
			modslab(scum2(x,y),max(hit->z-sq,zs),min(hit->z+sq+1,ze));
		}
	}
	scum2finish();
	updatebbox(vx5.minx,vx5.miny,vx5.minz,vx5.maxx,vx5.maxy,vx5.maxz,dacol);
}

void setellipsoid (lpoint3d *hit, lpoint3d *hit2, long hitrad, long dacol, long bakit)
{
	void (*modslab)(long *, long, long);
	long x, y, xs, ys, zs, xe, ye, ze;
	float a, b, c, d, e, f, g, h, r, t, u, Za, Zb, fx0, fy0, fz0, fx1, fy1, fz1;

	xs = min(hit->x,hit2->x)-hitrad; xs = max(xs,0);
	ys = min(hit->y,hit2->y)-hitrad; ys = max(ys,0);
	zs = min(hit->z,hit2->z)-hitrad; zs = max(zs,0);
	xe = max(hit->x,hit2->x)+hitrad; xe = min(xe,VSID-1);
	ye = max(hit->y,hit2->y)+hitrad; ye = min(ye,VSID-1);
	ze = max(hit->z,hit2->z)+hitrad; ze = min(ze,MAXZDIM-1);
	vx5.minx = xs; vx5.maxx = xe+1;
	vx5.miny = ys; vx5.maxy = ye+1;
	vx5.minz = zs; vx5.maxz = ze+1;
	if ((xs > xe) || (ys > ye) || (zs > ze))
		{ if (bakit) voxbackup(xs,ys,xs,ys,bakit); return; }

	fx0 = (float)hit->x; fy0 = (float)hit->y; fz0 = (float)hit->z;
	fx1 = (float)hit2->x; fy1 = (float)hit2->y; fz1 = (float)hit2->z;

	r = (fx1-fx0)*(fx1-fx0) + (fy1-fy0)*(fy1-fy0) + (fz1-fz0)*(fz1-fz0);
	r = sqrt((float)hitrad*(float)hitrad + r*.25);
	c = fz0*fz0 - fz1*fz1; d = r*r*-4; e = d*4;
	f = c*c + fz1*fz1 * e; g = c + c; h = (fz1-fz0)*2; c = c*h - fz1*e;
	Za = -h*h - e; if (Za <= 0) { if (bakit) voxbackup(xs,ys,xs,ys,bakit); return; }
	u = 1 / Za;

	if (vx5.colfunc == jitcolfunc) vx5.amount = 0x70707;

	if (dacol == -1) modslab = delslab; else modslab = insslab;

	if (bakit) voxbackup(xs,ys,xe+1,ye+1,bakit);

	for(y=ys;y<=ye;y++)
		for(x=xs;x<=xe;x++)
		{
			a = (x-fx0)*(x-fx0) + (y-fy0)*(y-fy0);
			b = (x-fx1)*(x-fx1) + (y-fy1)*(y-fy1);
			t = a-b+d; Zb = t*h + c;
			t = ((t+g)*t + b*e + f)*Za + Zb*Zb; if (t <= 0) continue;
			t = sqrt(t);
			ftol((Zb - t)*u,&zs); if (zs < 0) zs = 0;
			ftol((Zb + t)*u,&ze); if (ze > MAXZDIM) ze = MAXZDIM;
			modslab(scum2(x,y),zs,ze);
		}
	scum2finish();
	updatebbox(vx5.minx,vx5.miny,vx5.minz,vx5.maxx,vx5.maxy,vx5.maxz,dacol);
}

	//Draws a cylinder, given: 2 points, a radius, and a color
	//Code mostly optimized - original code from CYLINDER.BAS:drawcylinder
void setcylinder (lpoint3d *p0, lpoint3d *p1, long cr, long dacol, long bakit)
{
	void (*modslab)(long *, long, long);

	float t, ax, ay, az, bx, by, bz, cx, cy, cz, ux, uy, uz, vx, vy, vz;
	float Za, Zb, Zc, tcr, xxyy, rcz, rZa;
	float fx, fxi, xof, vx0, vy0, vz0, vz0i, vxo, vyo, vzo;
	long i, j, ix, iy, ix0, ix1, iz0, iz1, minx, maxx, miny, maxy;
	long x0, y0, z0, x1, y1, z1;

		//Map generic cylinder into unit space:  (0,0,0), (0,0,1), cr = 1
		//   x*x + y*y < 1, z >= 0, z < 1
	if (p0->z > p1->z)
	{
		x0 = p1->x; y0 = p1->y; z0 = p1->z;
		x1 = p0->x; y1 = p0->y; z1 = p0->z;
	}
	else
	{
		x0 = p0->x; y0 = p0->y; z0 = p0->z;
		x1 = p1->x; y1 = p1->y; z1 = p1->z;
	}

	xxyy = (float)((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0));
	t = xxyy + (float)(z1-z0)*(z1-z0);
	if ((t == 0) || (cr == 0))
	{
		vx5.minx = x0; vx5.maxx = x0+1;
		vx5.miny = y0; vx5.maxy = y0+1;
		vx5.minz = z0; vx5.maxz = z0+1;
		if (bakit) voxbackup(x0,y0,x0,y0,bakit);
		return;
	}
	t = 1 / t; cx = ((float)(x1-x0))*t; cy = ((float)(y1-y0))*t; cz = ((float)(z1-z0))*t;
	t = sqrt(t); ux = ((float)(x1-x0))*t; uy = ((float)(y1-y0))*t; uz = ((float)(z1-z0))*t;

	if (vx5.colfunc == jitcolfunc) vx5.amount = 0x70707;

	if (dacol == -1) modslab = delslab; else modslab = insslab;

	if (xxyy == 0)
	{
		iz0 = max(z0,0); iz1 = min(z1,MAXZDIM);
		minx = max(x0-cr,0); maxx = min(x0+cr,VSID-1);
		miny = max(y0-cr,0); maxy = min(y0+cr,VSID-1);

		vx5.minx = minx; vx5.maxx = maxx+1;
		vx5.miny = miny; vx5.maxy = maxy+1;
		vx5.minz = iz0; vx5.maxz = iz1;
		if (bakit) voxbackup(minx,miny,maxx+1,maxy+1,bakit);

		j = cr*cr;
		for(iy=miny;iy<=maxy;iy++)
		{
			i = j-(iy-y0)*(iy-y0);
			for(ix=minx;ix<=maxx;ix++)
				if ((ix-x0)*(ix-x0) < i) modslab(scum2(ix,iy),iz0,iz1);
		}
		scum2finish();
		updatebbox(vx5.minx,vx5.miny,vx5.minz,vx5.maxx,vx5.maxy,vx5.maxz,dacol);
		return;
	}

	if (x0 < x1) { minx = x0; maxx = x1; } else { minx = x1; maxx = x0; }
	if (y0 < y1) { miny = y0; maxy = y1; } else { miny = y1; maxy = y0; }
	tcr = cr / sqrt(xxyy); vx = fabs((float)(x1-x0))*tcr; vy = fabs((float)(y1-y0))*tcr;
	t = vx*uz + vy;
	ftol((float)minx-t,&minx); if (minx < 0) minx = 0;
	ftol((float)maxx+t,&maxx); if (maxx >= VSID) maxx = VSID-1;
	t = vy*uz + vx;
	ftol((float)miny-t,&miny); if (miny < 0) miny = 0;
	ftol((float)maxy+t,&maxy); if (maxy >= VSID) maxy = VSID-1;

	vx5.minx = minx; vx5.maxx = maxx+1;
	vx5.miny = miny; vx5.maxy = maxy+1;
	vx5.minz = z0-cr; vx5.maxz = z1+cr+1;
	if (bakit) voxbackup(minx,miny,maxx+1,maxy+1,bakit);

	vx = (fabs(ux) < fabs(uy)); vy = 1.0f-vx; vz = 0;
	ax = uy*vz - uz*vy; ay = uz*vx - ux*vz; az = ux*vy - uy*vx;
	t = 1.0 / (sqrt(ax*ax + ay*ay + az*az)*cr);
	ax *= t; ay *= t; az *= t;
	bx = ay*uz - az*uy; by = az*ux - ax*uz; bz = ax*uy - ay*ux;

	Za = az*az + bz*bz; rZa = 1.0f / Za;
	if (cz != 0) { rcz = 1.0f / cz; vz0i = -rcz*cx; }
	if (y0 != y1)
	{
		t = 1.0f / ((float)(y1-y0)); fxi = ((float)(x1-x0))*t;
		fx = ((float)miny-y0)*fxi + x0; xof = fabs(tcr*xxyy*t);
	}
	else { fx = (float)minx; fxi = 0.0; xof = (float)(maxx-minx); }

	vy = (float)(miny-y0);
	vxo = vy*ay - z0*az;
	vyo = vy*by - z0*bz;
	vzo = vy*cy - z0*cz;
	for(iy=miny;iy<=maxy;iy++)
	{
		ftol(fx-xof,&ix0); if (ix0 < minx) ix0 = minx;
		ftol(fx+xof,&ix1); if (ix1 > maxx) ix1 = maxx;
		fx += fxi;

		vx = (float)(ix0-x0);
		vx0 = vx*ax + vxo; vxo += ay;
		vy0 = vx*bx + vyo; vyo += by;
		vz0 = vx*cx + vzo; vzo += cy;

		if (cz != 0)   //(vx0 + vx1*t)ý + (vy0 + vy1*t)ý = 1
		{
			vz0 *= -rcz;
			for(ix=ix0;ix<=ix1;ix++,vx0+=ax,vy0+=bx,vz0+=vz0i)
			{
				Zb = vx0*az + vy0*bz; Zc = vx0*vx0 + vy0*vy0 - 1;
				t = Zb*Zb - Za*Zc; if (*(long *)&t <= 0) continue; t = sqrt(t);
				ftol(max((-Zb-t)*rZa,vz0    ),&iz0); if (iz0 < 0) iz0 = 0;
				ftol(min((-Zb+t)*rZa,vz0+rcz),&iz1); if (iz1 > MAXZDIM) iz1 = MAXZDIM;
				modslab(scum2(ix,iy),iz0,iz1);
			}
		}
		else
		{
			for(ix=ix0;ix<=ix1;ix++,vx0+=ax,vy0+=bx,vz0+=cx)
			{
				if (*(unsigned long *)&vz0 >= 0x3f800000) continue; //vz0<0||vz0>=1
				Zb = vx0*az + vy0*bz; Zc = vx0*vx0 + vy0*vy0 - 1;
				t = Zb*Zb - Za*Zc; if (*(long *)&t <= 0) continue; t = sqrt(t);
				ftol((-Zb-t)*rZa,&iz0); if (iz0 < 0) iz0 = 0;
				ftol((-Zb+t)*rZa,&iz1); if (iz1 > MAXZDIM) iz1 = MAXZDIM;
				modslab(scum2(ix,iy),iz0,iz1);
			}
		}
	}
	scum2finish();
	updatebbox(vx5.minx,vx5.miny,vx5.minz,vx5.maxx,vx5.maxy,vx5.maxz,dacol);
}

	//Draws a rectangle, given: 2 points as opposite corners, and a color
void setrect (lpoint3d *hit, lpoint3d *hit2, long dacol)
{
	long x, y, xs, ys, zs, xe, ye, ze;

		//WARNING: do NOT use lbound because 'c' not guaranteed to be >= 'b'
	xs = max(min(hit->x,hit2->x),0); xe = min(max(hit->x,hit2->x),VSID-1);
	ys = max(min(hit->y,hit2->y),0); ye = min(max(hit->y,hit2->y),VSID-1);
	zs = max(min(hit->z,hit2->z),0); ze = min(max(hit->z,hit2->z),MAXZDIM-1);
	vx5.minx = xs; vx5.maxx = xe+1;
	vx5.miny = ys; vx5.maxy = ye+1;
	vx5.minz = zs; vx5.maxz = ze+1;
	if ((xs > xe) || (ys > ye) || (zs > ze)) return;

	if (vx5.colfunc == jitcolfunc) vx5.amount = 0x70707;

	ze++;
	if (dacol == -1)
	{
		for(y=ys;y<=ye;y++)
			for(x=xs;x<=xe;x++)
				delslab(scum2(x,y),zs,ze);
	}
	else
	{
		for(y=ys;y<=ye;y++)
			for(x=xs;x<=xe;x++)
				insslab(scum2(x,y),zs,ze);
	}
	scum2finish();
	updatebbox(vx5.minx,vx5.miny,vx5.minz,vx5.maxx,vx5.maxy,vx5.maxz,dacol);
}

	//Does CSG using pre-sorted spanlist
void setspans (vspans *lst, long lstnum, lpoint3d *offs, long dacol)
{
	void (*modslab)(long *, long, long);
	long i, j, x, y, z0, z1, *lptr;
	char ox, oy;

	if (lstnum <= 0) return;
	if (dacol == -1) modslab = delslab; else modslab = insslab;
	vx5.minx = vx5.maxx = ((long)lst[0].x)+offs->x;
	vx5.miny = ((long)lst[       0].y)+offs->y;
	vx5.maxy = ((long)lst[lstnum-1].y)+offs->y+1;
	vx5.minz = vx5.maxz = ((long)lst[0].z0)+offs->z;

	i = 0; goto in2setlist;
	do
	{
		if ((ox != lst[i].x) || (oy != lst[i].y))
		{
in2setlist:;
			ox = lst[i].x; oy = lst[i].y;
			x = ((long)lst[i].x)+offs->x;
			y = ((long)lst[i].y)+offs->y;
				  if (x < vx5.minx) vx5.minx = x;
			else if (x > vx5.maxx) vx5.maxx = x;
			lptr = scum2(x,y);
		}
		if ((x|y)&(~(VSID-1))) { i++; continue; }
		z0 = ((long)lst[i].z0)+offs->z;   if (z0 < 0) z0 = 0;
		z1 = ((long)lst[i].z1)+offs->z+1; if (z1 > MAXZDIM) z1 = MAXZDIM;
		if (z0 < vx5.minz) vx5.minz = z0;
		if (z1 > vx5.maxz) vx5.maxz = z1;
		modslab(lptr,z0,z1);
		i++;
	} while (i < lstnum);
	vx5.maxx++; vx5.maxz++;
	if (vx5.minx < 0) vx5.minx = 0;
	if (vx5.miny < 0) vx5.miny = 0;
	if (vx5.maxx > VSID) vx5.maxx = VSID;
	if (vx5.maxy > VSID) vx5.maxy = VSID;

	scum2finish();
	updatebbox(vx5.minx,vx5.miny,vx5.minz,vx5.maxx,vx5.maxy,vx5.maxz,dacol);
}

void setheightmap (const unsigned char *hptr, long hpitch, long hxdim, long hydim,
						 long x0, long y0, long x1, long y1)
{
	long x, y, su, sv, u, v;

	if (x0 < 0) x0 = 0;
	if (y0 < 0) y0 = 0;
	if (x1 > VSID) x1 = VSID;
	if (y1 > VSID) y1 = VSID;
	vx5.minx = x0; vx5.maxx = x1;
	vx5.miny = y0; vx5.maxy = y1;
	vx5.minz = 0; vx5.maxz = MAXZDIM;
	if ((x0 >= x1) || (y0 >= y1)) return;

	su = x0%hxdim; sv = y0%hydim;
	for(y=y0,v=sv;y<y1;y++)
	{
		for(x=x0,u=su;x<x1;x++)
		{
			insslab(scum2(x,y),hptr[v*hpitch+u],MAXZDIM);
			u++; if (u >= hxdim) u = 0;
		}
		v++; if (v >= hydim) v = 0;
	}
	scum2finish();
	updatebbox(vx5.minx,vx5.miny,vx5.minz,vx5.maxx,vx5.maxy,vx5.maxz,0);
}

static long min0[VSID], max0[VSID]; //MAXY
static long min1[VSID], max1[VSID]; //MAXX
static long min2[VSID], max2[VSID]; //MAXY

static void canseerange (point3d *p0, point3d *p1)
{
	lpoint3d a, c, d, p, i;
	point3d f, g;
	long cnt, j;

	ftol(p0->x-.5,&a.x); ftol(p0->y-.5,&a.y); ftol(p0->z-.5,&a.z);
	ftol(p1->x-.5,&c.x); ftol(p1->y-.5,&c.y); ftol(p1->z-.5,&c.z);
	cnt = 0;

		  if (c.x <  a.x) { d.x = -1; f.x = p0->x-a.x;   g.x = (p0->x-p1->x)*1024; cnt += a.x-c.x; }
	else if (c.x != a.x) { d.x =  1; f.x = a.x+1-p0->x; g.x = (p1->x-p0->x)*1024; cnt += c.x-a.x; }
	else f.x = g.x = 0;
		  if (c.y <  a.y) { d.y = -1; f.y = p0->y-a.y;   g.y = (p0->y-p1->y)*1024; cnt += a.y-c.y; }
	else if (c.y != a.y) { d.y =  1; f.y = a.y+1-p0->y; g.y = (p1->y-p0->y)*1024; cnt += c.y-a.y; }
	else f.y = g.y = 0;
		  if (c.z <  a.z) { d.z = -1; f.z = p0->z-a.z;   g.z = (p0->z-p1->z)*1024; cnt += a.z-c.z; }
	else if (c.z != a.z) { d.z =  1; f.z = a.z+1-p0->z; g.z = (p1->z-p0->z)*1024; cnt += c.z-a.z; }
	else f.z = g.z = 0;

	ftol(f.x*g.z - f.z*g.x,&p.x); ftol(g.x,&i.x);
	ftol(f.y*g.z - f.z*g.y,&p.y); ftol(g.y,&i.y);
	ftol(f.y*g.x - f.x*g.y,&p.z); ftol(g.z,&i.z);
	for(;cnt;cnt--)
	{
			//use a.x, a.y, a.z
		if (a.x < min0[a.y]) min0[a.y] = a.x;
		if (a.x > max0[a.y]) max0[a.y] = a.x;
		if (a.z < min1[a.x]) min1[a.x] = a.z;
		if (a.z > max1[a.x]) max1[a.x] = a.z;
		if (a.z < min2[a.y]) min2[a.y] = a.z;
		if (a.z > max2[a.y]) max2[a.y] = a.z;

		if (((p.x|p.y) >= 0) && (a.z != c.z)) { a.z += d.z; p.x -= i.x; p.y -= i.y; }
		else if ((p.z >= 0) && (a.x != c.x))  { a.x += d.x; p.x += i.z; p.z -= i.y; }
		else                                  { a.y += d.y; p.y += i.z; p.z += i.x; }
	}
}

void settri (point3d *p0, point3d *p1, point3d *p2, long bakit)
{
	point3d n;
	float f, x0, y0, z0, x1, y1, z1, rx, ry, k0, k1;
	long i, x, y, z, iz0, iz1, minx, maxx, miny, maxy;

	if (p0->x < p1->x) { x0 = p0->x; x1 = p1->x; } else { x0 = p1->x; x1 = p0->x; }
	if (p2->x < x0) x0 = p2->x;
	if (p2->x > x1) x1 = p2->x;
	if (p0->y < p1->y) { y0 = p0->y; y1 = p1->y; } else { y0 = p1->y; y1 = p0->y; }
	if (p2->y < y0) y0 = p2->y;
	if (p2->y > y1) y1 = p2->y;
	if (p0->z < p1->z) { z0 = p0->z; z1 = p1->z; } else { z0 = p1->z; z1 = p0->z; }
	if (p2->z < z0) z0 = p2->z;
	if (p2->z > z1) z1 = p2->z;

	ftol(x0-.5,&minx); ftol(y0-.5,&miny);
	ftol(x1-.5,&maxx); ftol(y1-.5,&maxy);
	vx5.minx = minx; vx5.maxx = maxx+1;
	vx5.miny = miny; vx5.maxy = maxy+1;
	ftol(z0-.5,&vx5.minz); ftol(z1+.5,&vx5.maxz);
	if (bakit) voxbackup(minx,miny,maxx+1,maxy+1,bakit);

	for(i=miny;i<=maxy;i++) { min0[i] = 0x7fffffff; max0[i] = 0x80000000; }
	for(i=minx;i<=maxx;i++) { min1[i] = 0x7fffffff; max1[i] = 0x80000000; }
	for(i=miny;i<=maxy;i++) { min2[i] = 0x7fffffff; max2[i] = 0x80000000; }

	canseerange(p0,p1);
	canseerange(p1,p2);
	canseerange(p2,p0);

	n.x = (p1->z-p0->z)*(p2->y-p1->y) - (p1->y-p0->y) * (p2->z-p1->z);
	n.y = (p1->x-p0->x)*(p2->z-p1->z) - (p1->z-p0->z) * (p2->x-p1->x);
	n.z = (p1->y-p0->y)*(p2->x-p1->x) - (p1->x-p0->x) * (p2->y-p1->y);
	f = 1.0 / sqrt(n.x*n.x + n.y*n.y + n.z*n.z); if (n.z < 0) f = -f;
	n.x *= f; n.y *= f; n.z *= f;

	if (n.z > .01)
	{
		f = -1.0 / n.z; rx = n.x*f; ry = n.y*f;
		k0 = ((n.x>=0)-p0->x)*rx + ((n.y>=0)-p0->y)*ry - ((n.z>=0)-p0->z) + .5;
		k1 = ((n.x< 0)-p0->x)*rx + ((n.y< 0)-p0->y)*ry - ((n.z< 0)-p0->z) - .5;
	}
	else { rx = 0; ry = 0; k0 = -2147000000.0; k1 = 2147000000.0; }

	for(y=miny;y<=maxy;y++)
		for(x=min0[y];x<=max0[y];x++)
		{
			f = (float)x*rx + (float)y*ry; ftol(f+k0,&iz0); ftol(f+k1,&iz1);
			if (iz0 < min1[x]) iz0 = min1[x];
			if (iz1 > max1[x]) iz1 = max1[x];
			if (iz0 < min2[y]) iz0 = min2[y];
			if (iz1 > max2[y]) iz1 = max2[y];

				//set: (x,y,iz0) to (x,y,iz1) (inclusive)
			insslab(scum2(x,y),iz0,iz1+1);
	}
	scum2finish();
	updatebbox(vx5.minx,vx5.miny,vx5.minz,vx5.maxx,vx5.maxy,vx5.maxz,0);
}

// ------------------------ CONVEX 3D HULL CODE BEGINS ------------------------

char umost[VSID*VSID];

// ------------------------- CONVEX 3D HULL CODE ENDS -------------------------

	//Old&Slow sector code, but only this one supports the 3D bumpmapping :(
static void setsectorb (point3d *p, long *point2, long n, float thick, long dacol, long bakit, long bumpmap)
{
	point3d norm, p2;
	float d, f, x0, y0, x1, y1;
	long i, j, k, got, x, y, z, xs, ys, zs, xe, ye, ze, maxis, ndacol;

	norm.x = 0; norm.y = 0; norm.z = 0;
	for(i=0;i<n;i++)
	{
		j = point2[i]; k = point2[j];
		norm.x += (p[i].y-p[j].y)*(p[k].z-p[j].z) - (p[i].z-p[j].z)*(p[k].y-p[j].y);
		norm.y += (p[i].z-p[j].z)*(p[k].x-p[j].x) - (p[i].x-p[j].x)*(p[k].z-p[j].z);
		norm.z += (p[i].x-p[j].x)*(p[k].y-p[j].y) - (p[i].y-p[j].y)*(p[k].x-p[j].x);
	}
	f = 1.0 / sqrt(norm.x*norm.x + norm.y*norm.y + norm.z*norm.z);
	norm.x *= f; norm.y *= f; norm.z *= f;

	if ((fabs(norm.z) >= fabs(norm.x)) && (fabs(norm.z) >= fabs(norm.y)))
		maxis = 2;
	else if (fabs(norm.y) > fabs(norm.x))
		maxis = 1;
	else
		maxis = 0;

	xs = xe = p[0].x;
	ys = ye = p[0].y;
	zs = ze = p[0].z;
	for(i=n-1;i;i--)
	{
		if (p[i].x < xs) xs = p[i].x;
		if (p[i].y < ys) ys = p[i].y;
		if (p[i].z < zs) zs = p[i].z;
		if (p[i].x > xe) xe = p[i].x;
		if (p[i].y > ye) ye = p[i].y;
		if (p[i].z > ze) ze = p[i].z;
	}
	xs = max(xs-thick-bumpmap,0); xe = min(xe+thick+bumpmap,VSID-1);
	ys = max(ys-thick-bumpmap,0); ye = min(ye+thick+bumpmap,VSID-1);
	zs = max(zs-thick-bumpmap,0); ze = min(ze+thick+bumpmap,MAXZDIM-1);
	vx5.minx = xs; vx5.maxx = xe+1;
	vx5.miny = ys; vx5.maxy = ye+1;
	vx5.minz = zs; vx5.maxz = ze+1;
	if ((xs > xe) || (ys > ye) || (zs > ze)) return;
	if (bakit) voxbackup(xs,ys,xe+1,ye+1,bakit);

	clearbuf((void *)&templongbuf[zs],ze-zs+1,-3);

	ndacol = (dacol==-1)-2;

	for(y=ys;y<=ye;y++)
		for(x=xs;x<=xe;x++)
		{
			got = 0;
			d = ((float)x-p[0].x)*norm.x + ((float)y-p[0].y)*norm.y + ((float)zs-p[0].z)*norm.z;
			for(z=zs;z<=ze;z++,d+=norm.z)
			{
				if (bumpmap)
				{
					if (d < -thick) continue;
					p2.x = (float)x - d*norm.x;
					p2.y = (float)y - d*norm.y;
					p2.z = (float)z - d*norm.z;
					if (d > (float)hpngcolfunc(&p2)+thick) continue;
				}
				else
				{
					if (fabs(d) > thick) continue;
					p2.x = (float)x - d*norm.x;
					p2.y = (float)y - d*norm.y;
					p2.z = (float)z - d*norm.z;
				}

				k = 0;
				for(i=n-1;i>=0;i--)
				{
					j = point2[i];
					switch(maxis)
					{
						case 0: x0 = p[i].z-p2.z; x1 = p[j].z-p2.z;
								  y0 = p[i].y-p2.y; y1 = p[j].y-p2.y; break;
						case 1: x0 = p[i].x-p2.x; x1 = p[j].x-p2.x;
								  y0 = p[i].z-p2.z; y1 = p[j].z-p2.z; break;
						case 2: x0 = p[i].x-p2.x; x1 = p[j].x-p2.x;
								  y0 = p[i].y-p2.y; y1 = p[j].y-p2.y; break;
						default: __assume(0); //tells MSVC default can't be reached
					}
					if (((*(long *)&y0)^(*(long *)&y1)) < 0)
					{
						if (((*(long *)&x0)^(*(long *)&x1)) >= 0) k ^= (*(long *)&x0);
						else { f = (x0*y1-x1*y0); k ^= (*(long *)&f)^(*(long *)&y1); }
					}
				}
				if (k >= 0) continue;

				templongbuf[z] = ndacol; got = 1;
			}
			if (got)
			{
				scum(x,y,zs,ze+1,templongbuf);
				clearbuf((void *)&templongbuf[zs],ze-zs+1,-3);
			}
		}
	scumfinish();
	updatebbox(vx5.minx,vx5.miny,vx5.minz,vx5.maxx,vx5.maxy,vx5.maxz,dacol);
}

	//This is for ordfillpolygon&splitpoly
typedef struct { long p, i, t; } raster;
#define MAXCURS 100 //THIS IS VERY EVIL... FIX IT!!!
static raster rst[MAXCURS];
static long slist[MAXCURS];

	//Code taken from POLYOLD\POLYSPLI.BAS:splitpoly (06/09/2001)
void splitpoly (float *px, float *py, long *point2, long *bakn,
					 float x0, float y0, float dx, float dy)
{
	long i, j, s2, n, sn, splcnt, z0, z1, z2, z3;
	float t, t1;

	n = (*bakn); if (n < 3) return;
	i = 0; s2 = sn = n; splcnt = 0;
	do
	{
		t1 = (px[i]-x0)*dy - (py[i]-y0)*dx;
		do
		{
			j = point2[i]; point2[i] |= 0x80000000;
			t = t1; t1 = (px[j]-x0)*dy - (py[j]-y0)*dx;
			if ((*(long *)&t) < 0)
				{ px[n] = px[i]; py[n] = py[i]; point2[n] = n+1; n++; }
			if (((*(long *)&t) ^ (*(long *)&t1)) < 0)
			{
				if ((*(long *)&t) < 0) slist[splcnt++] = n;
				t /= (t-t1);
				px[n] = (px[j]-px[i])*t + px[i];
				py[n] = (py[j]-py[i])*t + py[i];
				point2[n] = n+1; n++;
			}
			i = j;
		} while (point2[i] >= 0);
		if (n > s2) { point2[n-1] = s2; s2 = n; }
		for(i=sn-1;(i) && (point2[i] < 0);i--);
	} while (i > 0);

	if (fabs(dx) > fabs(dy))
	{
		for(i=1;i<splcnt;i++)
		{
			z0 = slist[i];
			for(j=0;j<i;j++)
			{
				z1 = point2[z0]; z2 = slist[j]; z3 = point2[z2];
				if (fabs(px[z0]-px[z3])+fabs(px[z2]-px[z1]) < fabs(px[z0]-px[z1])+fabs(px[z2]-px[z3]))
					{ point2[z0] = z3; point2[z2] = z1; }
			}
		}
	}
	else
	{
		for(i=1;i<splcnt;i++)
		{
			z0 = slist[i];
			for(j=0;j<i;j++)
			{
				z1 = point2[z0]; z2 = slist[j]; z3 = point2[z2];
				if (fabs(py[z0]-py[z3])+fabs(py[z2]-py[z1]) < fabs(py[z0]-py[z1])+fabs(py[z2]-py[z3]))
					{ point2[z0] = z3; point2[z2] = z1; }
			}
		}
	}

	for(i=sn;i<n;i++)
		{ px[i-sn] = px[i]; py[i-sn] = py[i]; point2[i-sn] = point2[i]-sn; }
	(*bakn) = n-sn;
}

void ordfillpolygon (float *px, float *py, long *point2, long n, long day, long xs, long xe, void (*modslab)(long *, long, long))
{
	float f;
	long k, i, z, zz, z0, z1, zx, sx0, sy0, sx1, sy1, sy, nsy, gap, numrst;
	long np, ni;

	if (n < 3) return;

	for(z=0;z<n;z++) slist[z] = z;

		//Sort points by y's
	for(gap=(n>>1);gap;gap>>=1)
		for(z=0;z<n-gap;z++)
			for(zz=z;zz>=0;zz-=gap)
			{
				if (py[point2[slist[zz]]] <= py[point2[slist[zz+gap]]]) break;
				z0 = slist[zz]; slist[zz] = slist[zz+gap]; slist[zz+gap] = z0;
			}

	ftol(py[point2[slist[0]]]+.5,&sy); if (sy < xs) sy = xs;

	numrst = 0; z = 0; n--; //Note: n is local variable!
	while (z < n)
	{
		z1 = slist[z]; z0 = point2[z1];
		for(zx=0;zx<2;zx++)
		{
			ftol(py[z0]+.5,&sy0); ftol(py[z1]+.5,&sy1);
			if (sy1 > sy0) //Insert raster (z0,z1)
			{
				f = (px[z1]-px[z0]) / (py[z1]-py[z0]);
				ftol(((sy-py[z0])*f + px[z0])*65536.0 + 65535.0,&np);
				if (sy1-sy0 >= 2) ftol(f*65536.0,&ni); else ni = 0;
				k = (np<<1)+ni;
				for(i=numrst;i>0;i--)
				{
					if ((rst[i-1].p<<1)+rst[i-1].i < k) break;
					rst[i] = rst[i-1];
				}
				rst[i].i = ni; rst[i].p = np; rst[i].t = (z0<<16)+z1;
				numrst++;
			}
			else if (sy1 < sy0) //Delete raster (z1,z0)
			{
				numrst--;
				k = (z1<<16)+z0; i = 0;
				while ((i < numrst) && (rst[i].t != k)) i++;
				while (i < numrst) { rst[i] = rst[i+1]; i++; }
			}
			z1 = point2[z0];
		}

		z++;
		ftol(py[point2[slist[z]]]+.5,&nsy); if (nsy > xe) nsy = xe;
		for(;sy<nsy;sy++)
			for(i=0;i<numrst;i+=2)
			{
				modslab(scum2(sy,day),max(rst[i].p>>16,0),min(rst[i+1].p>>16,MAXZDIM));
				rst[i].p += rst[i].i; rst[i+1].p += rst[i+1].i;
			}
	}
}

	//Draws a flat polygon
	//given: p&point2: 3D points, n: # points, thick: thickness, dacol: color
static float ppx[MAXCURS*4], ppy[MAXCURS*4];
static long npoint2[MAXCURS*4];
void setsector (point3d *p, long *point2, long n, float thick, long dacol, long bakit)
{
	void (*modslab)(long *, long, long);
	point3d norm;
	float f, rnormy, xth, zth, dax, daz, t, t1;
	long i, j, k, x, y, z, sn, s2, nn, xs, ys, zs, xe, ye, ze;

	norm.x = 0; norm.y = 0; norm.z = 0;
	for(i=0;i<n;i++)
	{
		j = point2[i]; k = point2[j];
		norm.x += (p[i].y-p[j].y)*(p[k].z-p[j].z) - (p[i].z-p[j].z)*(p[k].y-p[j].y);
		norm.y += (p[i].z-p[j].z)*(p[k].x-p[j].x) - (p[i].x-p[j].x)*(p[k].z-p[j].z);
		norm.z += (p[i].x-p[j].x)*(p[k].y-p[j].y) - (p[i].y-p[j].y)*(p[k].x-p[j].x);
	}
	f = 1.0 / sqrt(norm.x*norm.x + norm.y*norm.y + norm.z*norm.z);
	norm.x *= f; norm.y *= f; norm.z *= f;

	if (vx5.colfunc == jitcolfunc) vx5.amount = 0x70707;
	else if ((vx5.colfunc == pngcolfunc) && (vx5.pic) && (vx5.xsiz > 0) && (vx5.ysiz > 0) && (vx5.picmode == 3))
	{
			//Find biggest height offset to minimize bounding box size
		j = k = vx5.pic[0];
		for(y=vx5.ysiz-1;y>=0;y--)
		{
			i = y*(vx5.bpl>>2);
			for(x=vx5.xsiz-1;x>=0;x--)
			{
				if (vx5.pic[i+x] < j) j = vx5.pic[i+x];
				if (vx5.pic[i+x] > k) k = vx5.pic[i+x];
			}
		}
		if ((j^k)&0xff000000) //If high bytes are !=, then use bumpmapping
		{
			setsectorb(p,point2,n,thick,dacol,bakit,max(labs(j>>24),labs(k>>24)));
			return;
		}
	}

	xs = xe = p[0].x;
	ys = ye = p[0].y;
	zs = ze = p[0].z;
	for(i=n-1;i;i--)
	{
		if (p[i].x < xs) xs = p[i].x;
		if (p[i].y < ys) ys = p[i].y;
		if (p[i].z < zs) zs = p[i].z;
		if (p[i].x > xe) xe = p[i].x;
		if (p[i].y > ye) ye = p[i].y;
		if (p[i].z > ze) ze = p[i].z;
	}
	xs = max(xs-thick,0); xe = min(xe+thick,VSID-1);
	ys = max(ys-thick,0); ye = min(ye+thick,VSID-1);
	zs = max(zs-thick,0); ze = min(ze+thick,MAXZDIM-1);
	vx5.minx = xs; vx5.maxx = xe+1;
	vx5.miny = ys; vx5.maxy = ye+1;
	vx5.minz = zs; vx5.maxz = ze+1;
	if ((xs > xe) || (ys > ye) || (zs > ze)) return;
	if (bakit) voxbackup(xs,ys,xe+1,ye+1,bakit);

	if (dacol == -1) modslab = delslab; else modslab = insslab;

	if (fabs(norm.y) >= .001)
	{
		rnormy = 1.0 / norm.y;
		for(y=ys;y<=ye;y++)
		{
			nn = n;
			for(i=0;i<n;i++)
			{
				f = ((float)y-p[i].y) * rnormy;
				ppx[i] = norm.z*f + p[i].z;
				ppy[i] = norm.x*f + p[i].x;
				npoint2[i] = point2[i];
			}
			if (fabs(norm.x) > fabs(norm.z))
			{
				splitpoly(ppx,ppy,npoint2,&nn,p[0].z,((p[0].y-(float)y)*norm.y-thick)/norm.x+p[0].x,norm.x,-norm.z);
				splitpoly(ppx,ppy,npoint2,&nn,p[0].z,((p[0].y-(float)y)*norm.y+thick)/norm.x+p[0].x,-norm.x,norm.z);
			}
			else
			{
				splitpoly(ppx,ppy,npoint2,&nn,((p[0].y-(float)y)*norm.y-thick)/norm.z+p[0].z,p[0].x,norm.x,-norm.z);
				splitpoly(ppx,ppy,npoint2,&nn,((p[0].y-(float)y)*norm.y+thick)/norm.z+p[0].z,p[0].x,-norm.x,norm.z);
			}
			ordfillpolygon(ppx,ppy,npoint2,nn,y,xs,xe,modslab);
		}
	}
	else
	{
		xth = norm.x*thick; zth = norm.z*thick;
		for(y=ys;y<=ye;y++)
		{
			for(z=0;z<n;z++) slist[z] = 0;
			nn = 0; i = 0; sn = n;
			do
			{
				s2 = nn; t1 = p[i].y-(float)y;
				do
				{
					j = point2[i]; slist[i] = 1; t = t1; t1 = p[j].y-(float)y;
					if (((*(long *)&t) ^ (*(long *)&t1)) < 0)
					{
						k = ((*(unsigned long *)&t)>>31); t /= (t-t1);
						daz = (p[j].z-p[i].z)*t + p[i].z;
						dax = (p[j].x-p[i].x)*t + p[i].x;
						ppx[nn+k] = daz+zth; ppx[nn+1-k] = daz-zth;
						ppy[nn+k] = dax+xth; ppy[nn+1-k] = dax-xth;
						npoint2[nn] = nn+1; npoint2[nn+1] = nn+2; nn += 2;
					}
					i = j;
				} while (!slist[i]);
				if (nn > s2) { npoint2[nn-1] = s2; s2 = nn; }
				for(i=sn-1;(i) && (slist[i]);i--);
			} while (i);
			ordfillpolygon(ppx,ppy,npoint2,nn,y,xs,xe,modslab);
		}
	}
	scum2finish();
	updatebbox(vx5.minx,vx5.miny,vx5.minz,vx5.maxx,vx5.maxy,vx5.maxz,dacol);
}

//FLOODFILL3D begins --------------------------------------------------------

#define FILLBUFSIZ 16384 //Use realloc instead!
typedef struct { unsigned short x, y, z0, z1; } spoint4d; //128K
static spoint4d fbuf[FILLBUFSIZ];

long dntil0 (long x, long y, long z)
{
	char *v = sptr[y*VSID+x];
	while (1)
	{
		if (z < v[1]) break;
		if (!v[0]) return(MAXZDIM);
		v += v[0]*4;
		if (z < v[3]) return(v[3]);
	}
	return(z);
}

long dntil1 (long x, long y, long z)
{
	char *v = sptr[y*VSID+x];
	while (1)
	{
		if (z <= v[1]) return(v[1]);
		if (!v[0]) break;
		v += v[0]*4;
		if (z < v[3]) break;
	}
	return(z);
}

long uptil1 (long x, long y, long z)
{
	char *v = sptr[y*VSID+x];
	if (z < v[1]) return(0);
	while (v[0])
	{
		v += v[0]*4;
		if (z < v[3]) break;
		if (z < v[1]) return(v[3]);
	}
	return(z);
}

void hollowfillstart (long x, long y, long z)
{
	spoint4d a;
	char *v;
	long i, j, z0, z1, i0, i1;

	a.x = x; a.y = y;

	v = sptr[y*VSID+x]; j = ((((long)v)-(long)vbuf)>>2); a.z0 = 0;
	while (1)
	{
		a.z1 = (long)(v[1]);
		if ((a.z0 <= z) && (z < a.z1) && (!(vbit[j>>5]&(1<<j)))) break;
		if (!v[0]) return;
		v += v[0]*4; j += 2;
		a.z0 = (long)(v[3]);
	}
	vbit[j>>5] |= (1<<j); //fill a.x,a.y,a.z0<=?<a.z1

	i0 = i1 = 0; goto floodfill3dskip2;
	do
	{
		a = fbuf[i0]; i0 = ((i0+1)&(FILLBUFSIZ-1));
floodfill3dskip2:;
		for(i=3;i>=0;i--)
		{
			if (i&1) { x = a.x+(i&2)-1; if ((unsigned long)x >= VSID) continue; y = a.y; }
				 else { y = a.y+(i&2)-1; if ((unsigned long)y >= VSID) continue; x = a.x; }

			v = sptr[y*VSID+x]; j = ((((long)v)-(long)vbuf)>>2); z0 = 0;
			while (1)
			{
				z1 = (long)(v[1]);
				if ((z0 < a.z1) && (a.z0 < z1) && (!(vbit[j>>5]&(1<<j))))
				{
					fbuf[i1].x = x; fbuf[i1].y = y;
					fbuf[i1].z0 = z0; fbuf[i1].z1 = z1;
					i1 = ((i1+1)&(FILLBUFSIZ-1));
					if (i0 == i1) return; //floodfill stack overflow!
					vbit[j>>5] |= (1<<j); //fill x,y,z0<=?<z1
				}
				if (!v[0]) break;
				v += v[0]*4; j += 2;
				z0 = (long)(v[3]);
			}
		}
	} while (i0 != i1);
}

//FLOODFILL3D ends ----------------------------------------------------------

#define LPATBUFSIZ 14
static lpoint2d *patbuf;
#define LPATHASHSIZ 12
static lpoint3d *pathashdat;
static long *pathashead, pathashcnt, pathashmax;

static void initpathash ()
{
	patbuf = (lpoint2d *)radar;
	pathashead = (long *)(((long)patbuf)+(1<<LPATBUFSIZ)*sizeof(lpoint2d));
	pathashdat = (lpoint3d *)(((long)pathashead)+((1<<LPATHASHSIZ)*4));
	pathashmax = ((max((MAXXDIM*MAXYDIM*27)>>1,(VSID+4)*3*256*4)-((1<<LPATBUFSIZ)*sizeof(lpoint2d))-(1<<LPATHASHSIZ)*4)/12);
	memset(pathashead,-1,(1<<LPATHASHSIZ)*4);
	pathashcnt = 0;
}

static long readpathash (long i)
{
	long j = (((i>>LPATHASHSIZ)-i) & ((1<<LPATHASHSIZ)-1));
	for(j=pathashead[j];j>=0;j=pathashdat[j].x)
		if (pathashdat[j].y == i) return(pathashdat[j].z);
	return(-1);
}

static void writepathash (long i, long v)
{
	long k, j = (((i>>LPATHASHSIZ)-i) & ((1<<LPATHASHSIZ)-1));
	for(k=pathashead[j];k>=0;k=pathashdat[k].x)
		if (pathashdat[k].y == i) { pathashdat[k].z = v; return; }
	pathashdat[pathashcnt].x = pathashead[j]; pathashead[j] = pathashcnt;
	pathashdat[pathashcnt].y = i;
	pathashdat[pathashcnt].z = v;
	pathashcnt++;
}

static signed char cdir[26*4] = //sqrt(2) =~ 58/41, sqrt(3) =~ 71/41;
{
	-1, 0, 0,41,  1, 0, 0,41,  0,-1, 0,41,  0, 1, 0,41,  0, 0,-1,41,  0, 0, 1,41,
	-1,-1, 0,58, -1, 1, 0,58, -1, 0,-1,58, -1, 0, 1,58,  0,-1,-1,58,  0,-1, 1,58,
	 1,-1, 0,58,  1, 1, 0,58,  1, 0,-1,58,  1, 0, 1,58,  0, 1,-1,58,  0, 1, 1,58,
	-1,-1,-1,71, -1,-1, 1,71, -1, 1,-1,71, -1, 1, 1,71,
	 1,-1,-1,71,  1,-1, 1,71,  1, 1,-1,71,  1, 1, 1,71,
};

long findpath (long *pathpos, long pathmax, lpoint3d *p1, lpoint3d *p0)
{
	long i, j, k, x, y, z, c, nc, xx, yy, zz, bufr, bufw, pcnt;

	if (!(getcube(p0->x,p0->y,p0->z)&~1))
	{
		for(i=5;i>=0;i--)
		{
			x = p0->x+(long)cdir[i*4]; y = p0->y+(long)cdir[i*4+1]; z = p0->z+(long)cdir[i*4+2];
			if (getcube(x,y,z)&~1) { p0->x = x; p0->y = y; p0->z = z; break; }
		}
		if (i < 0) return(0);
	}
	if (!(getcube(p1->x,p1->y,p1->z)&~1))
	{
		for(i=5;i>=0;i--)
		{
			x = p1->x+(long)cdir[i*4]; y = p1->y+(long)cdir[i*4+1]; z = p1->z+(long)cdir[i*4+2];
			if (getcube(x,y,z)&~1) { p1->x = x; p1->y = y; p1->z = z; break; }
		}
		if (i < 0) return(0);
	}

	initpathash();
	j = (p0->x*VSID + p0->y)*MAXZDIM+p0->z;
	patbuf[0].x = j; patbuf[0].y = 0; bufr = 0; bufw = 1;
	writepathash(j,0);
	do
	{
		j = patbuf[bufr&((1<<LPATBUFSIZ)-1)].x;
		x = j/(VSID*MAXZDIM); y = ((j/MAXZDIM)&(VSID-1)); z = (j&(MAXZDIM-1));
		c = patbuf[bufr&((1<<LPATBUFSIZ)-1)].y; bufr++;
		for(i=0;i<26;i++)
		{
			xx = x+(long)cdir[i*4]; yy = y+(long)cdir[i*4+1]; zz = z+(long)cdir[i*4+2];
			j = (xx*VSID + yy)*MAXZDIM+zz;

			//nc = c+(long)cdir[i*4+3]; //More accurate but lowers max distance a lot!
			//if (((k = getcube(xx,yy,zz))&~1) && ((unsigned long)nc < (unsigned long)readpathash(j)))

			if (((k = getcube(xx,yy,zz))&~1) && (readpathash(j) < 0))
			{
				nc = c+(long)cdir[i*4+3];
				if ((xx == p1->x) && (yy == p1->y) && (zz == p1->z)) { c = nc; goto pathfound; }
				writepathash(j,nc);
				if (pathashcnt >= pathashmax) return(0);
				patbuf[bufw&((1<<LPATBUFSIZ)-1)].x = (xx*VSID + yy)*MAXZDIM+zz;
				patbuf[bufw&((1<<LPATBUFSIZ)-1)].y = nc; bufw++;
			}
		}
	} while (bufr != bufw);

pathfound:
	if (pathmax <= 0) return(0);
	pathpos[0] = (p1->x*VSID + p1->y)*MAXZDIM+p1->z; pcnt = 1;
	x = p1->x; y = p1->y; z = p1->z;
	do
	{
		for(i=0;i<26;i++)
		{
			xx = x+(long)cdir[i*4]; yy = y+(long)cdir[i*4+1]; zz = z+(long)cdir[i*4+2];
			nc = c-(long)cdir[i*4+3];
			if (readpathash((xx*VSID + yy)*MAXZDIM+zz) == nc)
			{
				if (pcnt >= pathmax) return(0);
				pathpos[pcnt] = (xx*VSID + yy)*MAXZDIM+zz; pcnt++;
				x = xx; y = yy; z = zz; c = nc; break;
			}
		}
	} while (i < 26);
	if (pcnt >= pathmax) return(0);
	pathpos[pcnt] = (p0->x*VSID + p0->y)*MAXZDIM+p0->z;
	return(pcnt+1);
}

//---------------------------------------------------------------------

	//This here for game programmer only. I would never use it!
void drawpoint2d (long sx, long sy, long col)
{
	if ((unsigned long)sx >= (unsigned long)xres) return;
	if ((unsigned long)sy >= (unsigned long)yres) return;
	*(long *)(ylookup[sy]+(sx<<2)+frameplace) = col;
}

	//This here for game programmer only. I would never use it!
void drawpoint3d (float x0, float y0, float z0, long col)
{
	float ox, oy, oz, r;
	long x, y;

	ox = x0-gipos.x; oy = y0-gipos.y; oz = z0-gipos.z;
	z0 = ox*gifor.x + oy*gifor.y + oz*gifor.z; if (z0 < SCISDIST) return;
	r = 1.0f / z0;
	x0 = (ox*gistr.x + oy*gistr.y + oz*gistr.z)*gihz;
	y0 = (ox*gihei.x + oy*gihei.y + oz*gihei.z)*gihz;

	ftol(x0*r + gihx-.5f,&x); if ((unsigned long)x >= (unsigned long)xres) return;
	ftol(y0*r + gihy-.5f,&y); if ((unsigned long)y >= (unsigned long)yres) return;
	*(long *)(ylookup[y]+(x<<2)+frameplace) = col;
}

	//returns 1 if visible
long project2d (float x, float y, float z, float *px, float *py, float *sx)
{
	float ox, oy, oz;

	ox = x-gipos.x; oy = y-gipos.y; oz = z-gipos.z;
	z = ox*gifor.x + oy*gifor.y + oz*gifor.z; if (z < SCISDIST) return(0);
	
	z = gihz / z;
	*px = (ox*gistr.x + oy*gistr.y + oz*gistr.z)*z + gihx;
	*py = (ox*gihei.x + oy*gihei.y + oz*gihei.z)*z + gihy;
	*sx = z;   
	return(1);
}

static __int64 mskp255 = 0x00ff00ff00ff00ff;
static __int64 mskn255 = 0xff01ff01ff01ff01;
static __int64 rgbmask64 = 0xffffff00ffffff;

void drawline2d (float x1, float y1, float x2, float y2, long col)
{
	float dx, dy, fxresm1, fyresm1;
	long i, j, incr, ie;

	dx = x2-x1; dy = y2-y1; if ((dx == 0) && (dy == 0)) return;
	fxresm1 = (float)xres-.5; fyresm1 = (float)yres-.5;
	if (x1 >= fxresm1) { if (x2 >= fxresm1) return; y1 += (fxresm1-x1)*dy/dx; x1 = fxresm1; }
	else if (x1 < 0) { if (x2 < 0) return; y1 += (0-x1)*dy/dx; x1 = 0; }
	if (x2 >= fxresm1) { y2 += (fxresm1-x2)*dy/dx; x2 = fxresm1; }
	else if (x2 < 0) { y2 += (0-x2)*dy/dx; x2 = 0; }
	if (y1 >= fyresm1) { if (y2 >= fyresm1) return; x1 += (fyresm1-y1)*dx/dy; y1 = fyresm1; }
	else if (y1 < 0) { if (y2 < 0) return; x1 += (0-y1)*dx/dy; y1 = 0; }
	if (y2 >= fyresm1) { x2 += (fyresm1-y2)*dx/dy; y2 = fyresm1; }
	else if (y2 < 0) { x2 += (0-y2)*dx/dy; y2 = 0; }

	if (fabs(dx) >= fabs(dy))
	{
		if (x2 > x1) { ftol(x1,&i); ftol(x2,&ie); } else { ftol(x2,&i); ftol(x1,&ie); }
		if (i < 0) i = 0; if (ie >= xres) ie = xres-1;
		ftol(1048576.0*dy/dx,&incr); ftol(y1*1048576.0+((float)i+.5f-x1)*incr,&j);
		for(;i<=ie;i++,j+=incr)
			if ((unsigned long)(j>>20) < (unsigned long)yres)
				*(long *)(ylookup[j>>20]+(i<<2)+frameplace) = col;
	}
	else
	{
		if (y2 > y1) { ftol(y1,&i); ftol(y2,&ie); } else { ftol(y2,&i); ftol(y1,&ie); }
		if (i < 0) i = 0; if (ie >= yres) ie = yres-1;
		ftol(1048576.0*dx/dy,&incr); ftol(x1*1048576.0+((float)i+.5f-y1)*incr,&j);
		for(;i<=ie;i++,j+=incr)
			if ((unsigned long)(j>>20) < (unsigned long)xres)
				*(long *)(ylookup[i]+((j>>18)&~3)+frameplace) = col;
	}
}

#if (USEZBUFFER == 1)
void drawline2dclip (float x1, float y1, float x2, float y2, float rx0, float ry0, float rz0, float rx1, float ry1, float rz1, long col)
{
	float dx, dy, fxresm1, fyresm1, Za, Zb, Zc, z;
	long i, j, incr, ie, p;

	dx = x2-x1; dy = y2-y1; if ((dx == 0) && (dy == 0)) return;
	fxresm1 = (float)xres-.5; fyresm1 = (float)yres-.5;
	if (x1 >= fxresm1) { if (x2 >= fxresm1) return; y1 += (fxresm1-x1)*dy/dx; x1 = fxresm1; }
	else if (x1 < 0) { if (x2 < 0) return; y1 += (0-x1)*dy/dx; x1 = 0; }
	if (x2 >= fxresm1) { y2 += (fxresm1-x2)*dy/dx; x2 = fxresm1; }
	else if (x2 < 0) { y2 += (0-x2)*dy/dx; x2 = 0; }
	if (y1 >= fyresm1) { if (y2 >= fyresm1) return; x1 += (fyresm1-y1)*dx/dy; y1 = fyresm1; }
	else if (y1 < 0) { if (y2 < 0) return; x1 += (0-y1)*dx/dy; y1 = 0; }
	if (y2 >= fyresm1) { x2 += (fyresm1-y2)*dx/dy; y2 = fyresm1; }
	else if (y2 < 0) { x2 += (0-y2)*dx/dy; y2 = 0; }

	if (fabs(dx) >= fabs(dy))
	{
			//Original equation: (rz1*t+rz0) / (rx1*t+rx0) = gihz/(sx-gihx)
		Za = gihz*(rx0*rz1 - rx1*rz0); Zb = rz1; Zc = -gihx*rz1 - gihz*rx1;

		if (x2 > x1) { ftol(x1,&i); ftol(x2,&ie); } else { ftol(x2,&i); ftol(x1,&ie); }
		if (i < 0) i = 0; if (ie >= xres) ie = xres-1;
		ftol(1048576.0*dy/dx,&incr); ftol(y1*1048576.0+((float)i+.5f-x1)*incr,&j);
		for(;i<=ie;i++,j+=incr)
			if ((unsigned long)(j>>20) < (unsigned long)yres)
			{
				p = ylookup[j>>20]+(i<<2)+frameplace;
				z = Za / ((float)i*Zb + Zc);
				if (*(long *)&z >= *(long *)(p+zbufoff)) continue;
				*(long *)(p+zbufoff) = *(long *)&z;
				*(long *)p = col;
			}
	}
	else
	{
		Za = gihz*(ry0*rz1 - ry1*rz0); Zb = rz1; Zc = -gihy*rz1 - gihz*ry1;

		if (y2 > y1) { ftol(y1,&i); ftol(y2,&ie); } else { ftol(y2,&i); ftol(y1,&ie); }
		if (i < 0) i = 0; if (ie >= yres) ie = yres-1;
		ftol(1048576.0*dx/dy,&incr); ftol(x1*1048576.0+((float)i+.5f-y1)*incr,&j);
		for(;i<=ie;i++,j+=incr)
			if ((unsigned long)(j>>20) < (unsigned long)xres)
			{
				p = ylookup[i]+((j>>18)&~3)+frameplace;
				z = Za / ((float)i*Zb + Zc);
				if (*(long *)&z >= *(long *)(p+zbufoff)) continue;
				*(long *)(p+zbufoff) = *(long *)&z;
				*(long *)p = col;
			}
	}
}
#endif

void drawline3d (float x0, float y0, float z0, float x1, float y1, float z1, long col)
{
	float ox, oy, oz, r;

	ox = x0-gipos.x; oy = y0-gipos.y; oz = z0-gipos.z;
	x0 = ox*gistr.x + oy*gistr.y + oz*gistr.z;
	y0 = ox*gihei.x + oy*gihei.y + oz*gihei.z;
	z0 = ox*gifor.x + oy*gifor.y + oz*gifor.z;

	ox = x1-gipos.x; oy = y1-gipos.y; oz = z1-gipos.z;
	x1 = ox*gistr.x + oy*gistr.y + oz*gistr.z;
	y1 = ox*gihei.x + oy*gihei.y + oz*gihei.z;
	z1 = ox*gifor.x + oy*gifor.y + oz*gifor.z;

	if (z0 < SCISDIST)
	{
		if (z1 < SCISDIST) return;
		r = (SCISDIST-z0)/(z1-z0); z0 = SCISDIST;
		x0 += (x1-x0)*r; y0 += (y1-y0)*r;
	}
	else if (z1 < SCISDIST)
	{
		r = (SCISDIST-z1)/(z1-z0); z1 = SCISDIST;
		x1 += (x1-x0)*r; y1 += (y1-y0)*r;
	}

	ox = gihz/z0;
	oy = gihz/z1;

#if (USEZBUFFER == 1)
	if (!(col&0xff000000))
		drawline2dclip(x0*ox+gihx,y0*ox+gihy,x1*oy+gihx,y1*oy+gihy,x0,y0,z0,x1-x0,y1-y0,z1-z0,col);
	else
		drawline2d(x0*ox+gihx,y0*ox+gihy,x1*oy+gihx,y1*oy+gihy,col&0xffffff);
#else
	drawline2d(x0*ox+gihx,y0*ox+gihy,x1*oy+gihx,y1*oy+gihy,col);
#endif
}

void drawpicinquad (long rpic, long rbpl, long rxsiz, long rysiz,
						  long wpic, long wbpl, long wxsiz, long wysiz,
						  float x0, float y0, float x1, float y1,
						  float x2, float y2, float x3, float y3)
{
	float px[4], py[4], k0, k1, k2, k3, k4, k5, k6, k7, k8;
	float t, u, v, dx, dy, l0, l1, m0, m1, m2, n0, n1, n2, r;
	long i, j, k, l, imin, imax, sx, sxe, sy, sy1, dd, uu, vv, ddi, uui, vvi;
	long x, xi, *p, *pe, uvmax, iu, iv;

	px[0] = x0; px[1] = x1; px[2] = x2; px[3] = x3;
	py[0] = y0; py[1] = y1; py[2] = y2; py[3] = y3;

		//This code projects 4 point2D's into a t,u,v screen-projection matrix
		//
		//Derivation: (given 4 known (sx,sy,kt,ku,kv) pairs, solve for k0-k8)
		//   kt = k0*sx + k1*sy + k2
		//   ku = k3*sx + k4*sy + k5
		//   kv = k6*sx + k7*sy + k8
		//0 = (k3*x0 + k4*y0 + k5) / (k0*x0 + k1*y0 + k2) / rxsiz
		//0 = (k6*x0 + k7*y0 + k8) / (k0*x0 + k1*y0 + k2) / rysiz
		//1 = (k3*x1 + k4*y1 + k5) / (k0*x1 + k1*y1 + k2) / rxsiz
		//0 = (k6*x1 + k7*y1 + k8) / (k0*x1 + k1*y1 + k2) / rysiz
		//1 = (k3*x2 + k4*y2 + k5) / (k0*x2 + k1*y2 + k2) / rxsiz
		//1 = (k6*x2 + k7*y2 + k8) / (k0*x2 + k1*y2 + k2) / rysiz
		//0 = (k3*x3 + k4*y3 + k5) / (k0*x3 + k1*y3 + k2) / rxsiz
		//1 = (k6*x3 + k7*y3 + k8) / (k0*x3 + k1*y3 + k2) / rysiz
		//   40*, 28+, 1~, 30W
	k3 = y3 - y0; k4 = x0 - x3; k5 = x3*y0 - x0*y3;
	k6 = y0 - y1; k7 = x1 - x0; k8 = x0*y1 - x1*y0;
	n0 = x2*y3 - x3*y2; n1 = x3*y1 - x1*y3; n2 = x1*y2 - x2*y1;
	l0 = k6*x2 + k7*y2 + k8;
	l1 = k3*x2 + k4*y2 + k5;
	t = n0 + n1 + n2; dx = (float)rxsiz*t*l0; dy = (float)rysiz*t*l1;
	t = l0*l1;
	l0 *= (k3*x1 + k4*y1 + k5);
	l1 *= (k6*x3 + k7*y3 + k8);
	m0 = l1 - t; m1 = l0 - l1; m2 = t - l0;
	k0 = m0*y1 + m1*y2 + m2*y3;
	k1 = -(m0*x1 + m1*x2 + m2*x3);
	k2 = n0*l0 + n1*t + n2*l1;
	k3 *= dx; k4 *= dx; k5 *= dx;
	k6 *= dy; k7 *= dy; k8 *= dy;

		//Make sure k's are in good range for conversion to integers...
	t = fabs(k0);
	if (fabs(k1) > t) t = fabs(k1);
	if (fabs(k2) > t) t = fabs(k2);
	if (fabs(k3) > t) t = fabs(k3);
	if (fabs(k4) > t) t = fabs(k4);
	if (fabs(k5) > t) t = fabs(k5);
	if (fabs(k6) > t) t = fabs(k6);
	if (fabs(k7) > t) t = fabs(k7);
	if (fabs(k8) > t) t = fabs(k8);
	t = -268435456.0 / t;
	k0 *= t; k1 *= t; k2 *= t;
	k3 *= t; k4 *= t; k5 *= t;
	k6 *= t; k7 *= t; k8 *= t;
	ftol(k0,&ddi);

	imin = 0; imax = 0;
	for(i=1;i<4;i++)
	{
		if (py[i] < py[imin]) imin = i;
		if (py[i] > py[imax]) imax = i;
	}

	uvmax = (rysiz-1)*rbpl + (rxsiz<<2);

	i = imax;
	do
	{
		j = ((i+1)&3);
			//offset would normally be -.5, but need to bias by +1.0
		ftol(py[j]+.5,&sy); if (sy < 0) sy = 0;
		ftol(py[i]+.5,&sy1); if (sy1 > wysiz) sy1 = wysiz;
		if (sy1 > sy)
		{
			ftol((px[i]-px[j])*4096.0/(py[i]-py[j]),&xi);
			ftol(((float)sy-py[j])*(float)xi + px[j]*4096.0 + 4096.0,&x);
			for(;sy<sy1;sy++,x+=xi) lastx[sy] = (x>>12);
		}
		i = j;
	} while (i != imin);
	do
	{
		j = ((i+1)&3);
			//offset would normally be -.5, but need to bias by +1.0
		ftol(py[i]+.5,&sy); if (sy < 0) sy = 0;
		ftol(py[j]+.5,&sy1); if (sy1 > wysiz) sy1 = wysiz;
		if (sy1 > sy)
		{
			ftol((px[j]-px[i])*4096.0/(py[j]-py[i]),&xi);
			ftol(((float)sy-py[i])*(float)xi + px[i]*4096.0 + 4096.0,&x);
			for(;sy<sy1;sy++,x+=xi)
			{
				sx = lastx[sy]; if (sx < 0) sx = 0;
				sxe = (x>>12); if (sxe > wxsiz) sxe = wxsiz;
				if (sx >= sxe) continue;
				t = k0*(float)sx + k1*(float)sy + k2; r = 1.0 / t;
				u = k3*(float)sx + k4*(float)sy + k5; ftol(u*r-.5,&iu);
				v = k6*(float)sx + k7*(float)sy + k8; ftol(v*r-.5,&iv);
				ftol(t,&dd);
				ftol((float)iu*k0 - k3,&uui); ftol((float)iu*t - u,&uu);
				ftol((float)iv*k0 - k6,&vvi); ftol((float)iv*t - v,&vv);
				if (k3*t < u*k0) k =    -4; else { uui = -(uui+ddi); uu = -(uu+dd); k =    4; }
				if (k6*t < v*k0) l = -rbpl; else { vvi = -(vvi+ddi); vv = -(vv+dd); l = rbpl; }
				iu = iv*rbpl + (iu<<2);
				p  = (long *)(sy*wbpl+(sx<<2)+wpic);
				pe = (long *)(sy*wbpl+(sxe<<2)+wpic);
				do
				{
					if ((unsigned long)iu < uvmax) p[0] = *(long *)(rpic+iu);
					dd += ddi;
					uu += uui; while (uu < 0) { iu += k; uui -= ddi; uu -= dd; }
					vv += vvi; while (vv < 0) { iu += l; vvi -= ddi; vv -= dd; }
					p++;
				} while (p < pe);
			}
		}
		i = j;
	} while (i != imax);
}

__declspec(align(16)) static float dpqdistlut[MAXXDIM];
__declspec(align(16)) static float dpqmulval[4] = {0,1,2,3}, dpqfour[4] = {4,4,4,4};
__declspec(align(8)) static float dpq3dn[4];
void drawpolyquad (long rpic, long rbpl, long rxsiz, long rysiz,
						 float x0, float y0, float z0, float u0, float v0,
						 float x1, float y1, float z1, float u1, float v1,
						 float x2, float y2, float z2, float u2, float v2,
						 float x3, float y3, float z3)
{
	point3d fp, fp2;
	float px[6], py[6], pz[6], pu[6], pv[6], px2[4], py2[4], pz2[4], pu2[4], pv2[4];
	float f, t, u, v, r, nx, ny, nz, ox, oy, oz, scaler;
	float dx, dy, db, ux, uy, ub, vx, vy, vb;
	long i, j, k, l, imin, imax, sx, sxe, sy, sy1;
	long x, xi, *p, *pe, uvmax, iu, iv, n;
	long dd, uu, vv, ddi, uui, vvi, distlutoffs;

	px2[0] = x0; py2[0] = y0; pz2[0] = z0; pu2[0] = u0; pv2[0] = v0;
	px2[1] = x1; py2[1] = y1; pz2[1] = z1; pu2[1] = u1; pv2[1] = v1;
	px2[2] = x2; py2[2] = y2; pz2[2] = z2; pu2[2] = u2; pv2[2] = v2;
	px2[3] = x3; py2[3] = y3; pz2[3] = z3;

		//Calculate U-V coordinate of 4th point on quad (based on 1st 3)
	nx = (y1-y0)*(z2-z0) - (z1-z0)*(y2-y0);
	ny = (z1-z0)*(x2-x0) - (x1-x0)*(z2-z0);
	nz = (x1-x0)*(y2-y0) - (y1-y0)*(x2-x0);
	if ((fabs(nx) > fabs(ny)) && (fabs(nx) > fabs(nz)))
	{     //(y1-y0)*u + (y2-y0)*v = (y3-y0)
			//(z1-z0)*u + (z2-z0)*v = (z3-z0)
		f = 1/nx;
		u = ((y3-y0)*(z2-z0) - (z3-z0)*(y2-y0))*f;
		v = ((y1-y0)*(z3-z0) - (z1-z0)*(y3-y0))*f;
	}
	else if (fabs(ny) > fabs(nz))
	{     //(x1-x0)*u + (x2-x0)*v = (x3-x0)
			//(z1-z0)*u + (z2-z0)*v = (z3-z0)
		f = -1/ny;
		u = ((x3-x0)*(z2-z0) - (z3-z0)*(x2-x0))*f;
		v = ((x1-x0)*(z3-z0) - (z1-z0)*(x3-x0))*f;
	}
	else
	{     //(x1-x0)*u + (x2-x0)*v = (x3-x0)
			//(y1-y0)*u + (y2-y0)*v = (y3-y0)
		f = 1/nz;
		u = ((x3-x0)*(y2-y0) - (y3-y0)*(x2-x0))*f;
		v = ((x1-x0)*(y3-y0) - (y1-y0)*(x3-x0))*f;
	}
	pu2[3] = (u1-u0)*u + (u2-u0)*v + u0;
	pv2[3] = (v1-v0)*u + (v2-v0)*v + v0;


	for(i=4-1;i>=0;i--) //rotation
	{
		fp.x = px2[i]-gipos.x; fp.y = py2[i]-gipos.y; fp.z = pz2[i]-gipos.z;
		px2[i] = fp.x*gistr.x + fp.y*gistr.y + fp.z*gistr.z;
		py2[i] = fp.x*gihei.x + fp.y*gihei.y + fp.z*gihei.z;
		pz2[i] = fp.x*gifor.x + fp.y*gifor.y + fp.z*gifor.z;
	}

		//Clip to SCISDIST plane
	n = 0;
	for(i=0;i<4;i++)
	{
		j = ((i+1)&3);
		if (pz2[i] >= SCISDIST) { px[n] = px2[i]; py[n] = py2[i]; pz[n] = pz2[i]; pu[n] = pu2[i]; pv[n] = pv2[i]; n++; }
		if ((pz2[i] >= SCISDIST) != (pz2[j] >= SCISDIST))
		{
			f = (SCISDIST-pz2[i])/(pz2[j]-pz2[i]);
			px[n] = (px2[j]-px2[i])*f + px2[i];
			py[n] = (py2[j]-py2[i])*f + py2[i];
			pz[n] = SCISDIST;
			pu[n] = (pu2[j]-pu2[i])*f + pu2[i];
			pv[n] = (pv2[j]-pv2[i])*f + pv2[i]; n++;
		}
	}
	if (n < 3) return;

	for(i=n-1;i>=0;i--) //projection
	{
		pz[i] = 1/pz[i]; f = pz[i]*gihz;
		px[i] = px[i]*f + gihx;
		py[i] = py[i]*f + gihy;
	}

		//General equations:
		//pz[i] = (px[i]*gdx + py[i]*gdy + gdo)
		//pu[i] = (px[i]*gux + py[i]*guy + guo)/pz[i]
		//pv[i] = (px[i]*gvx + py[i]*gvy + gvo)/pz[i]
		//
		//px[0]*gdx + py[0]*gdy + 1*gdo = pz[0]
		//px[1]*gdx + py[1]*gdy + 1*gdo = pz[1]
		//px[2]*gdx + py[2]*gdy + 1*gdo = pz[2]
		//
		//px[0]*gux + py[0]*guy + 1*guo = pu[0]*pz[0] (pu[i] premultiplied by pz[i] above)
		//px[1]*gux + py[1]*guy + 1*guo = pu[1]*pz[1]
		//px[2]*gux + py[2]*guy + 1*guo = pu[2]*pz[2]
		//
		//px[0]*gvx + py[0]*gvy + 1*gvo = pv[0]*pz[0] (pv[i] premultiplied by pz[i] above)
		//px[1]*gvx + py[1]*gvy + 1*gvo = pv[1]*pz[1]
		//px[2]*gvx + py[2]*gvy + 1*gvo = pv[2]*pz[2]
	pu[0] *= pz[0]; pu[1] *= pz[1]; pu[2] *= pz[2];
	pv[0] *= pz[0]; pv[1] *= pz[1]; pv[2] *= pz[2];
	ox = py[1]-py[2]; oy = py[2]-py[0]; oz = py[0]-py[1];
	r = 1.0 / (ox*px[0] + oy*px[1] + oz*px[2]);
	dx = (ox*pz[0] + oy*pz[1] + oz*pz[2])*r;
	ux = (ox*pu[0] + oy*pu[1] + oz*pu[2])*r;
	vx = (ox*pv[0] + oy*pv[1] + oz*pv[2])*r;
	ox = px[2]-px[1]; oy = px[0]-px[2]; oz = px[1]-px[0];
	dy = (ox*pz[0] + oy*pz[1] + oz*pz[2])*r;
	uy = (ox*pu[0] + oy*pu[1] + oz*pu[2])*r;
	vy = (ox*pv[0] + oy*pv[1] + oz*pv[2])*r;
	db = pz[0] - px[0]*dx - py[0]*dy;
	ub = pu[0] - px[0]*ux - py[0]*uy;
	vb = pv[0] - px[0]*vx - py[0]*vy;

#if 1
		//Make sure k's are in good range for conversion to integers...
	t = fabs(ux);
	if (fabs(uy) > t) t = fabs(uy);
	if (fabs(ub) > t) t = fabs(ub);
	if (fabs(vx) > t) t = fabs(vx);
	if (fabs(vy) > t) t = fabs(vy);
	if (fabs(vb) > t) t = fabs(vb);
	if (fabs(dx) > t) t = fabs(dx);
	if (fabs(dy) > t) t = fabs(dy);
	if (fabs(db) > t) t = fabs(db);
	scaler = -268435456.0 / t;
	ux *= scaler; uy *= scaler; ub *= scaler;
	vx *= scaler; vy *= scaler; vb *= scaler;
	dx *= scaler; dy *= scaler; db *= scaler;
	ftol(dx,&ddi);
	uvmax = (rysiz-1)*rbpl + (rxsiz<<2);

	scaler = 1.f/scaler; t = dx*scaler;
	if (cputype&(1<<25))
	{
		_asm //SSE
		{
			movss xmm6, t         ;xmm6: -,-,-,dx*scaler
			shufps xmm6, xmm6, 0  ;xmm6: dx*scaler,dx*scaler,dx*scaler,dx*scaler
			movaps xmm7, xmm6     ;xmm7: dx*scaler,dx*scaler,dx*scaler,dx*scaler
			mulps xmm6, dpqmulval ;xmm6: dx*scaler*3,dx*scaler*2,dx*scaler*1,0
			mulps xmm7, dpqfour   ;xmm7: dx*scaler*4,dx*scaler*4,dx*scaler*4,dx*scaler*4
		}
	}
	else { dpq3dn[0] = 0; dpq3dn[1] = t; dpq3dn[2] = dpq3dn[3] = t+t; } //3DNow!
#endif

	imin = (py[1]<py[0]); imax = 1-imin;
	for(i=n-1;i>1;i--)
	{
		if (py[i] < py[imin]) imin = i;
		if (py[i] > py[imax]) imax = i;
	}

	i = imax;
	do
	{
		j = i+1; if (j >= n) j = 0;
			//offset would normally be -.5, but need to bias by +1.0
		ftol(py[j]+.5,&sy); if (sy < 0) sy = 0;
		ftol(py[i]+.5,&sy1); if (sy1 > yres) sy1 = yres;
		if (sy1 > sy)
		{
			ftol((px[i]-px[j])*4096.0/(py[i]-py[j]),&xi);
			ftol(((float)sy-py[j])*(float)xi + px[j]*4096.0 + 4096.0,&x);
			for(;sy<sy1;sy++,x+=xi) lastx[sy] = (x>>12);
		}
		i = j;
	} while (i != imin);
	do
	{
		j = i+1; if (j >= n) j = 0;
			//offset would normally be -.5, but need to bias by +1.0
		ftol(py[i]+.5,&sy); if (sy < 0) sy = 0;
		ftol(py[j]+.5,&sy1); if (sy1 > yres) sy1 = yres;
		if (sy1 > sy)
		{
			ftol((px[j]-px[i])*4096.0/(py[j]-py[i]),&xi);
			ftol(((float)sy-py[i])*(float)xi + px[i]*4096.0 + 4096.0,&x);
			for(;sy<sy1;sy++,x+=xi)
			{
				sx = lastx[sy]; if (sx < 0) sx = 0;
				sxe = (x>>12); if (sxe > xres) sxe = xres;
				if (sx >= sxe) continue;
				p  = (long *)(sy*bytesperline+(sx<<2)+frameplace);
				pe = (long *)(sy*bytesperline+(sxe<<2)+frameplace);
#if 0
					//Brute force
				do
				{
					f = 1.f/(dx*(float)sx + dy*(float)sy + db);
					if (f < *(float *)(((long)p)+zbufoff))
					{
						*(float *)(((long)p)+zbufoff) = f;
						ftol((ux*(float)sx + uy*(float)sy + ub)*f-.5,&iu);
						ftol((vx*(float)sx + vy*(float)sy + vb)*f-.5,&iv);
						if ((unsigned long)iu >= rxsiz) iu = 0;
						if ((unsigned long)iv >= rysiz) iv = 0;
						p[0] = *(long *)(iv*rbpl+(iu<<2)+rpic);
					}
					p++; sx++;
				} while (p < pe);
#else
					//Optimized (in C) hyperbolic texture-mapping (Added Z-buffer using SSE/3DNow! for recip's)
				t = dx*(float)sx + dy*(float)sy + db; r = 1.0 / t;
				u = ux*(float)sx + uy*(float)sy + ub; ftol(u*r-.5,&iu);
				v = vx*(float)sx + vy*(float)sy + vb; ftol(v*r-.5,&iv);
				ftol(t,&dd);
				ftol((float)iu*dx - ux,&uui); ftol((float)iu*t - u,&uu);
				ftol((float)iv*dx - vx,&vvi); ftol((float)iv*t - v,&vv);
				if (ux*t < u*dx) k =    -4; else { uui = -(uui+ddi); uu = -(uu+dd); k =    4; }
				if (vx*t < v*dx) l = -rbpl; else { vvi = -(vvi+ddi); vv = -(vv+dd); l = rbpl; }
				iu = iv*rbpl + (iu<<2);

				t *= scaler;
				_asm
				{
					mov ecx, sxe
					sub ecx, sx
					xor eax, eax
					lea ecx, [ecx*4]
					sub eax, ecx
					add ecx, offset dpqdistlut

					test cputype, 1 shl 25
					jz short dpqpre3dn

					movss xmm0, t ;dd+ddi*3 dd+ddi*2 dd+ddi*1 dd+ddi*0
					shufps xmm0, xmm0, 0
					addps xmm0, xmm6
	 dpqbegsse: rcpps xmm1, xmm0
					addps xmm0, xmm7
					movaps [eax+ecx], xmm1
					add eax, 16
					jl short dpqbegsse
					jmp short dpqendit

	 dpqpre3dn: movd mm0, t ;dd+ddi*1 dd+ddi*0
					punpckldq mm0, mm0
					pfadd mm0, dpq3dn[0]
					movq mm7, dpq3dn[8]
	 dpqbeg3dn: pswapd mm2, mm0
					pfrcp mm1, mm0     ;mm1: 1/mm0l 1/mm0l
					pfrcp mm2, mm2     ;mm2: 1/mm0h 1/mm0h
					punpckldq mm1, mm2 ;mm1: 1/mm0h 1/mm0l
					pfadd mm0, mm7
					movq [eax+ecx], mm1
					add eax, 8
					jl short dpqbeg3dn
					femms
	  dpqendit:
				}
				distlutoffs = ((long)dpqdistlut)-((long)p);
				do
				{
#if (USEZBUFFER != 0)
					if (*(long *)(((long)p)+zbufoff) > *(long *)(((long)p)+distlutoffs))
					{
						*(long *)(((long)p)+zbufoff) = *(long *)(((long)p)+distlutoffs);
#endif
						if ((unsigned long)iu < uvmax) p[0] = *(long *)(rpic+iu);
#if (USEZBUFFER != 0)
					}
#endif
					dd += ddi;
					uu += uui; while (uu < 0) { iu += k; uui -= ddi; uu -= dd; }
					vv += vvi; while (vv < 0) { iu += l; vvi -= ddi; vv -= dd; }
					p++;
				} while (p < pe);
#endif
			}
		}
		i = j;
	} while (i != imax);
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

//--------------------------  Name hash code begins --------------------------

	//khashbuf format: (used by getkv6/getkfa to avoid duplicate loads)
	//[long index to next hash or -1][pointer to struct][char type]string[\0]
	//[long index to next hash or -1][pointer to struct][chat type]string[\0]
	//...
	//type:0 = kv6data
	//type:1 = kfatype
#define KHASHINITSIZE 8192
static char *khashbuf = 0;
static long khashead[256], khashpos = 0, khashsiz = 0;

char *getkfilname (long namoff) { return(&khashbuf[namoff]); }

	//Returns: 0,retptr=-1: Error! (bad filename or out of memory)
	//         0,retptr>=0: Not in hash; new name allocated, valid index
	//         1,retptr>=0: Already in hash, valid index
	//   Uses a 256-entry hash to compare names very quickly.
static long inkhash (const char *filnam, long *retind)
{
	long i, j, hashind;

	(*retind) = -1;

	if (!filnam) return(0);
	j = strlen(filnam); if (!j) return(0);
	j += 10;
	if (khashpos+j > khashsiz) //Make sure string fits in khashbuf
	{
		i = khashsiz; do { i <<= 1; } while (khashpos+j > i);
		if (!(khashbuf = (char *)realloc(khashbuf,i))) return(0);
		khashsiz = i;
	}

		//Copy filename to avoid destroying original string
		//Also, calculate hash index (which hopefully is uniformly random :)
	strcpy(&khashbuf[khashpos+9],filnam);
	for(i=khashpos+9,hashind=0;khashbuf[i];i++)
	{
		if ((khashbuf[i] >= 'a') && (khashbuf[i] <= 'z')) khashbuf[i] -= 32;
		if (khashbuf[i] == '/') khashbuf[i] = '\\';
		hashind = (khashbuf[i] - hashind*3);
	}
	hashind %= (sizeof(khashead)/sizeof(khashead[0]));

		//Find if string is already in hash...
	for(i=khashead[hashind];i>=0;i=(*(long *)&khashbuf[i]))
		if (!strcmp(&khashbuf[i+9],&khashbuf[khashpos+9]))
			{ (*retind) = i; return(1); } //Early out: already in hash

	(*retind) = khashpos;
	*(long *)&khashbuf[khashpos] = khashead[hashind];
	*(long *)&khashbuf[khashpos+4] = 0; //Set to 0 just in case load fails
	khashead[hashind] = khashpos; khashpos += j;
	return(0);
}

//-------------------------- KV6 sprite code begins --------------------------

//EQUIVEC code begins -----------------------------------------------------
point3d univec[256];
__declspec(align(8)) short iunivec[256][4];

typedef struct
{
	float fibx[45], fiby[45];
	float azval[20], zmulk, zaddk;
	long fib[47], aztop, npoints;
} equivectyp;
static equivectyp equivec;

#ifdef _MSC_VER

static _inline long dmulshr0 (long a, long d, long s, long t)
{
	_asm
	{
		mov eax, a
		imul d
		mov ecx, eax
		mov eax, s
		imul t
		add eax, ecx
	}
}

#endif

void equiind2vec (long i, float *x, float *y, float *z)
{
	float r;
	(*z) = (float)i*equivec.zmulk + equivec.zaddk; r = sqrt(1.f - (*z)*(*z));
	fcossin((float)i*(GOLDRAT*PI*2),x,y); (*x) *= r; (*y) *= r;
}

	//Very fast; good quality
long equivec2indmem (float x, float y, float z)
{
	long b, i, j, k, bestc;
	float xy, zz, md, d;

	xy = atan2(y,x); //atan2 is 150 clock cycles!
	j = ((*(long *)&z)&0x7fffffff);
	bestc = equivec.aztop;
	do
	{
		if (j < *(long *)&equivec.azval[bestc]) break;
		bestc--;
	} while (bestc);

	zz = z + 1.f;
	ftol(equivec.fibx[bestc]*xy + equivec.fiby[bestc]*zz - .5,&i);
	bestc++;
	ftol(equivec.fibx[bestc]*xy + equivec.fiby[bestc]*zz - .5,&j);

	k = dmulshr0(equivec.fib[bestc+2],i,equivec.fib[bestc+1],j);
	if ((unsigned long)k < equivec.npoints)
	{
		md = univec[k].x*x + univec[k].y*y + univec[k].z*z;
		j = k;
	} else md = -2.f;
	b = bestc+3;
	do
	{
		i = equivec.fib[b] + k;
		if ((unsigned long)i < equivec.npoints)
		{
			d = univec[i].x*x + univec[i].y*y + univec[i].z*z;
			if (*(long *)&d > *(long *)&md) { md = d; j = i; }
		}
		b--;
	} while (b != bestc);
	return(j);
}

void equivecinit (long n)
{
	float t0, t1;
	long z;

		//Init constants for ind2vec
	equivec.npoints = n;
	equivec.zmulk = 2 / (float)n; equivec.zaddk = equivec.zmulk*.5 - 1.0;

		//equimemset
	for(z=n-1;z>=0;z--)
		equiind2vec(z,&univec[z].x,&univec[z].y,&univec[z].z);
	if (n&1) //Hack for when n=255 and want a <0,0,0> vector
		{ univec[n].x = univec[n].y = univec[n].z = 0; }

		//Init Fibonacci table
	equivec.fib[0] = 0; equivec.fib[1] = 1;
	for(z=2;z<47;z++) equivec.fib[z] = equivec.fib[z-2]+equivec.fib[z-1];

		//Init fibx/y LUT
	t0 = .5 / PI; t1 = (float)n * -.5;
	for(z=0;z<45;z++)
	{
		t0 = -t0; equivec.fibx[z] = (float)equivec.fib[z+2]*t0;
		t1 = -t1; equivec.fiby[z] = ((float)equivec.fib[z+2]*GOLDRAT - (float)equivec.fib[z])*t1;
	}

	t0 = 1 / ((float)n * PI);
	for(equivec.aztop=0;equivec.aztop<20;equivec.aztop++)
	{
		t1 = 1 - (float)equivec.fib[(equivec.aztop<<1)+6]*t0; if (t1 < 0) break;
		equivec.azval[equivec.aztop+1] = sqrt(t1);
	}
}

//EQUIVEC code ends -------------------------------------------------------

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

static __declspec(align(8)) short lightlist[MAXLIGHTS+1][4];
static __int64 all32767 = 0x7fff7fff7fff7fff;

#endif

//-------------------------- KFA sprite code begins --------------------------

static kv6voxtype *getvptr (kv6data *kv, long x, long y)
{
	kv6voxtype *v;
	long i, j;

	v = kv->vox;
	if ((x<<1) < kv->xsiz) { for(i=0         ;i< x;i++) v += kv->xlen[i]; }
	else { v += kv->numvoxs; for(i=kv->xsiz-1;i>=x;i--) v -= kv->xlen[i]; }
	j = x*kv->ysiz;
	if ((y<<1) < kv->ysiz) { for(i=0         ;i< y;i++) v += kv->ylen[j+i]; }
	else { v += kv->xlen[x]; for(i=kv->ysiz-1;i>=y;i--) v -= kv->ylen[j+i]; }
	return(v);
}

#define VFIFSIZ 16384 //SHOULDN'T BE STATIC ALLOCATION!!!
static long vfifo[VFIFSIZ];
static void floodsucksprite (vx5sprite *spr, kv6data *kv, long ox, long oy,
									  kv6voxtype *v0, kv6voxtype *v1)
{
	kv6voxtype *v, *ve, *ov, *v2, *v3;
	kv6data *kv6;
	long i, j, x, y, z, x0, y0, z0, x1, y1, z1, n, vfif0, vfif1;

	x0 = x1 = ox; y0 = y1 = oy; z0 = v0->z; z1 = v1->z;

	n = (((long)v1)-((long)v0))/sizeof(kv6voxtype)+1;
	v1->vis &= ~64;

	vfifo[0] = ox; vfifo[1] = oy;
	vfifo[2] = (long)v0; vfifo[3] = (long)v1;
	vfif0 = 0; vfif1 = 4;

	while (vfif0 < vfif1)
	{
		i = (vfif0&(VFIFSIZ-1)); vfif0 += 4;
		ox = vfifo[i]; oy = vfifo[i+1];
		v0 = (kv6voxtype *)vfifo[i+2]; v1 = (kv6voxtype *)vfifo[i+3];

		if (ox < x0) x0 = ox;
		if (ox > x1) x1 = ox;
		if (oy < y0) y0 = oy;
		if (oy > y1) y1 = oy;
		if (v0->z < z0) z0 = v0->z;
		if (v1->z > z1) z1 = v1->z;
		for(v=v0;v<=v1;v++) v->vis |= 128; //Mark as part of current piece

		for(j=0;j<4;j++)
		{
			switch(j)
			{
				case 0: x = ox-1; y = oy; break;
				case 1: x = ox+1; y = oy; break;
				case 2: x = ox; y = oy-1; break;
				case 3: x = ox; y = oy+1; break;
				default: __assume(0); //tells MSVC default can't be reached
			}
			if ((unsigned long)x >= kv->xsiz) continue;
			if ((unsigned long)y >= kv->ysiz) continue;

			v = getvptr(kv,x,y);
			for(ve=&v[kv->ylen[x*kv->ysiz+y]];v<ve;v++)
			{
				if (v->vis&16) ov = v;
				if (((v->vis&(64+32)) == 64+32) && (v0->z <= v->z) && (v1->z >= ov->z))
				{
					i = (vfif1&(VFIFSIZ-1)); vfif1 += 4;
					if (vfif1-vfif0 >= VFIFSIZ) //FIFO Overflow... make entire object 1 piece :/
					{
						for(i=kv->numvoxs-1;i>=0;i--)
						{
							if ((kv->vox[i].vis&(64+32)) == 64+32) { v1 = &kv->vox[i]; v1->vis &= ~64; }
							if (kv->vox[i].vis&16) for(v=&kv->vox[i];v<=v1;v++) kv->vox[i].vis |= 128;
						}
						x0 = y0 = z0 = 0; x1 = kv->xsiz; y1 = kv->ysiz; z1 = kv->zsiz; n = kv->numvoxs;
						goto floodsuckend;
					}
					vfifo[i] = x; vfifo[i+1] = y;
					vfifo[i+2] = (long)ov; vfifo[i+3] = (long)v;
					n += (((long)v)-((long)ov))/sizeof(kv6voxtype)+1;
					v->vis &= ~64;
				}
			}
		}
	}
	x1++; y1++; z1++;
floodsuckend:;

	i = sizeof(kv6data) + n*sizeof(kv6voxtype) + (x1-x0)*4 + (x1-x0)*(y1-y0)*2;
	if (!(kv6 = (kv6data *)malloc(i))) return;
	kv6->leng = i;
	kv6->xsiz = x1-x0;
	kv6->ysiz = y1-y0;
	kv6->zsiz = z1-z0;
	kv6->xpiv = 0; //Set limb pivots to 0 - don't need it twice!
	kv6->ypiv = 0;
	kv6->zpiv = 0;
	kv6->numvoxs = n;
	kv6->namoff = 0;
	kv6->lowermip = 0;
	kv6->vox = (kv6voxtype *)(((long)kv6)+sizeof(kv6data));
	kv6->xlen = (unsigned long *)(((long)kv6->vox)+n*sizeof(kv6voxtype));
	kv6->ylen = (unsigned short *)(((long)kv6->xlen)+(x1-x0)*4);

		//Extract sub-KV6 to newly allocated kv6data
	v3 = kv6->vox; n = 0;
	for(x=0,v=kv->vox;x<x0;x++) v += kv->xlen[x];
	for(;x<x1;x++)
	{
		v2 = v; ox = n;
		for(y=0;y<y0;y++) v += kv->ylen[x*kv->ysiz+y];
		for(;y<y1;y++)
		{
			oy = n;
			for(ve=&v[kv->ylen[x*kv->ysiz+y]];v<ve;v++)
				if (v->vis&128)
				{
					v->vis &= ~128;
					(*v3) = (*v);
					v3->z -= z0;
					v3++; n++;
				}
			kv6->ylen[(x-x0)*(y1-y0)+(y-y0)] = n-oy;
		}
		kv6->xlen[x-x0] = n-ox;
		v = v2+kv->xlen[x];
	}

	spr->p.x = x0-kv->xpiv;
	spr->p.y = y0-kv->ypiv;
	spr->p.z = z0-kv->zpiv;
	spr->s.x = 1; spr->s.y = 0; spr->s.z = 0;
	spr->h.x = 0; spr->h.y = 1; spr->h.z = 0;
	spr->f.x = 0; spr->f.y = 0; spr->f.z = 1;
	spr->voxnum = kv6;
	spr->flags = 0;
}

static char *stripdir (char *filnam)
{
	long i, j;
	for(i=0,j=-1;filnam[i];i++)
		if ((filnam[i] == '/') || (filnam[i] == '\\')) j = i;
	return(&filnam[j+1]);
}

static void kfasorthinge (hingetype *h, long nh, long *hsort)
{
	long i, j, n;

		//First pass: stick hinges with parent=-1 at end
	n = nh; j = 0;
	for(i=n-1;i>=0;i--)
	{
		if (h[i].parent < 0) hsort[--n] = i;
							 else hsort[j++] = i;
	}
		//Finish accumulation (n*log(n) if tree is perfectly balanced)
	while (n > 0)
	{
		i--; if (i < 0) i = n-1;
		j = hsort[i];
		if (h[h[j].parent].parent < 0)
		{
			h[j].parent = -2-h[j].parent; n--;
			hsort[i] = hsort[n]; hsort[n] = j;
		}
	}
		//Restore parents to original values
	for(i=nh-1;i>=0;i--) h[i].parent = -2-h[i].parent;
}

	//Given vector a, returns b&c that makes (a,b,c) orthonormal
void genperp (point3d *a, point3d *b, point3d *c)
{
	float t;

	if ((a->x == 0) && (a->y == 0) && (a->z == 0))
		{ b->x = 0; b->y = 0; b->z = 0; return; }
	if ((fabs(a->x) < fabs(a->y)) && (fabs(a->x) < fabs(a->z)))
		{ t = 1.0 / sqrt(a->y*a->y + a->z*a->z); b->x = 0; b->y = a->z*t; b->z = -a->y*t; }
	else if (fabs(a->y) < fabs(a->z))
		{ t = 1.0 / sqrt(a->x*a->x + a->z*a->z); b->x = -a->z*t; b->y = 0; b->z = a->x*t; }
	else
		{ t = 1.0 / sqrt(a->x*a->x + a->y*a->y); b->x = a->y*t; b->y = -a->x*t; b->z = 0; }
	c->x = a->y*b->z - a->z*b->y;
	c->y = a->z*b->x - a->x*b->z;
	c->z = a->x*b->y - a->y*b->x;
}

//--------------------------- KFA sprite code ends ---------------------------

#if 0

void setkv6 (vx5sprite *spr)
{
	point3d r0, r1;
	long x, y, vx, vy, vz;
	kv6data *kv;
	kv6voxtype *v, *ve;

	if (spr->flags&2) return;
	kv = spr->voxnum; if (!kv) return;

	vx5.minx = ?; vx5.maxx = ?+1;
	vx5.miny = ?; vx5.maxy = ?+1;
	vx5.minz = ?; vx5.maxz = ?+1;

	v = kv->vox; //.01 is to fool rounding so they aren't all even numbers
	r0.x = spr->p.x - kv->xpiv*spr->s.x - kv->ypiv*spr->h.x - kv->zpiv*spr->f.x - .01;
	r0.y = spr->p.y - kv->xpiv*spr->s.y - kv->ypiv*spr->h.y - kv->zpiv*spr->f.y - .01;
	r0.z = spr->p.z - kv->xpiv*spr->s.z - kv->ypiv*spr->h.z - kv->zpiv*spr->f.z - .01;
	vx5.colfunc = curcolfunc;
	for(x=0;x<kv->xsiz;x++)
	{
		r1 = r0;
		for(y=0;y<kv->ysiz;y++)
		{
			for(ve=&v[kv->ylen[x*kv->ysiz+y]];v<ve;v++)
			{
				ftol(spr->f.x*v->z + r1.x,&vx);
				ftol(spr->f.y*v->z + r1.y,&vy);
				ftol(spr->f.z*v->z + r1.z,&vz);
				vx5.curcol = ((v->col&0xffffff)|0x80000000);
				setcube(vx,vy,vz,-2);
			}
			r1.x += spr->h.x; r1.y += spr->h.y; r1.z += spr->h.z;
		}
		r0.x += spr->s.x; r0.y += spr->s.y; r0.z += spr->s.z;
	}
}

#else

#ifdef _MSC_VER

	//dmulshr22 = ((a*b + c*d)>>22)
static _inline long dmulshr22 (long a, long b, long c, long d)
{
	_asm
	{
		mov eax, a
		imul b
		mov ecx, eax
		push edx
		mov eax, c
		imul d
		add eax, ecx
		pop ecx
		adc edx, ecx
		shrd eax, edx, 22
	}
}


#endif

static kv6data *gfrezkv;
static lpoint3d gfrezx, gfrezy, gfrezz, gfrezp;
static signed char gkv6colx[27] = {0,  0, 0, 0, 0, 1,-1, -1,-1,-1,-1, 0, 0, 0, 0, 1, 1, 1, 1,  1, 1, 1, 1,-1,-1,-1,-1};
static signed char gkv6coly[27] = {0,  0, 0, 1,-1, 0, 0,  0, 0,-1, 1, 1, 1,-1,-1, 0, 0,-1, 1,  1, 1,-1,-1, 1, 1,-1,-1};
static signed char gkv6colz[27] = {0,  1,-1, 0, 0, 0, 0, -1, 1, 0, 0, 1,-1, 1,-1, 1,-1, 0, 0,  1,-1, 1,-1, 1,-1, 1,-1};
long kv6colfunc (lpoint3d *p)
{
	kv6voxtype *v0, *v1, *v, *ve;
	long i, j, k, x, y, z, ox, oy, nx, ny, nz, mind, d;

	x = ((p->x*gfrezx.x + p->y*gfrezy.x + p->z*gfrezz.x + gfrezp.x)>>16);
	y = ((p->x*gfrezx.y + p->y*gfrezy.y + p->z*gfrezz.y + gfrezp.y)>>16);
	z = ((p->x*gfrezx.z + p->y*gfrezy.z + p->z*gfrezz.z + gfrezp.z)>>16);
	x = lbound0(x,gfrezkv->xsiz-1);
	y = lbound0(y,gfrezkv->ysiz-1);
	z = lbound0(z,gfrezkv->zsiz-1);

		//Process x
	v0 = gfrezkv->vox;
	if ((x<<1) < gfrezkv->xsiz) { ox = oy = j = 0; }
	else { v0 += gfrezkv->numvoxs; ox = gfrezkv->xsiz; oy = gfrezkv->ysiz; j = ox*oy; }
	v1 = v0;

	for(k=0;k<27;k++)
	{
		nx = ((long)gkv6colx[k])+x; if ((unsigned long)nx >= gfrezkv->xsiz) continue;
		ny = ((long)gkv6coly[k])+y; if ((unsigned long)ny >= gfrezkv->ysiz) continue;
		nz = ((long)gkv6colz[k])+z; if ((unsigned long)nz >= gfrezkv->zsiz) continue;

		if (nx != ox)
		{
			while (nx > ox) { v0 += gfrezkv->xlen[ox]; ox++; j += gfrezkv->ysiz; }
			while (nx < ox) { ox--; v0 -= gfrezkv->xlen[ox]; j -= gfrezkv->ysiz; }
			if ((ny<<1) < gfrezkv->ysiz) { oy = 0; v1 = v0; }
			else { oy = gfrezkv->ysiz; v1 = v0+gfrezkv->xlen[nx]; }
		}
		if (ny != oy)
		{
			while (ny > oy) { v1 += gfrezkv->ylen[j+oy]; oy++; }
			while (ny < oy) { oy--; v1 -= gfrezkv->ylen[j+oy]; }
		}

			//Process z
		for(v=v1,ve=&v1[gfrezkv->ylen[j+ny]];v<ve;v++)
			if (v->z == nz) return(v->col);
	}

		//Use brute force when all else fails.. :/
	v = gfrezkv->vox; mind = 0x7fffffff;
	for(nx=0;nx<gfrezkv->xsiz;nx++)
		for(ny=0;ny<gfrezkv->ysiz;ny++)
			for(ve=&v[gfrezkv->ylen[nx*gfrezkv->ysiz+ny]];v<ve;v++)
			{
				d = labs(x-nx)+labs(y-ny)+labs(z-v->z);
				if (d < mind) { mind = d; k = v->col; }
			}
	return(k);
}

static void kv6colfuncinit (vx5sprite *spr, float det)
{
	point3d tp, tp2;
	float f;

	gfrezkv = spr->voxnum; if (!gfrezkv) { vx5.colfunc = curcolfunc; return; }

	tp2.x = gfrezkv->xpiv + .5;
	tp2.y = gfrezkv->ypiv + .5;
	tp2.z = gfrezkv->zpiv + .5;
	tp.x = spr->p.x - spr->s.x*tp2.x - spr->h.x*tp2.y - spr->f.x*tp2.z;
	tp.y = spr->p.y - spr->s.y*tp2.x - spr->h.y*tp2.y - spr->f.y*tp2.z;
	tp.z = spr->p.z - spr->s.z*tp2.x - spr->h.z*tp2.y - spr->f.z*tp2.z;

		//spr->s.x*x + spr->h.x*y + spr->f.x*z = np.x; //Solve for x,y,z
		//spr->s.y*x + spr->h.y*y + spr->f.y*z = np.y;
		//spr->s.z*x + spr->h.z*y + spr->f.z*z = np.z;
	f = 65536.0 / det;

	tp2.x = (spr->h.y*spr->f.z - spr->h.z*spr->f.y)*f; ftol(tp2.x,&gfrezx.x);
	tp2.y = (spr->h.z*spr->f.x - spr->h.x*spr->f.z)*f; ftol(tp2.y,&gfrezy.x);
	tp2.z = (spr->h.x*spr->f.y - spr->h.y*spr->f.x)*f; ftol(tp2.z,&gfrezz.x);
	ftol(-tp.x*tp2.x - tp.y*tp2.y - tp.z*tp2.z,&gfrezp.x); gfrezp.x += 32767;

	tp2.x = (spr->f.y*spr->s.z - spr->f.z*spr->s.y)*f; ftol(tp2.x,&gfrezx.y);
	tp2.y = (spr->f.z*spr->s.x - spr->f.x*spr->s.z)*f; ftol(tp2.y,&gfrezy.y);
	tp2.z = (spr->f.x*spr->s.y - spr->f.y*spr->s.x)*f; ftol(tp2.z,&gfrezz.y);
	ftol(-tp.x*tp2.x - tp.y*tp2.y - tp.z*tp2.z,&gfrezp.y); gfrezp.y += 32767;

	tp2.x = (spr->s.y*spr->h.z - spr->s.z*spr->h.y)*f; ftol(tp2.x,&gfrezx.z);
	tp2.y = (spr->s.z*spr->h.x - spr->s.x*spr->h.z)*f; ftol(tp2.y,&gfrezy.z);
	tp2.z = (spr->s.x*spr->h.y - spr->s.y*spr->h.x)*f; ftol(tp2.z,&gfrezz.z);
	ftol(-tp.x*tp2.x - tp.y*tp2.y - tp.z*tp2.z,&gfrezp.z); gfrezp.z += 32768;
}

#define LSC3 8 //2 for testing, 8 is normal
typedef struct
{
	long xo, yo, zo, xu, yu, zu, xv, yv, zv, d, msk, pzi;
	long xmino, ymino, xmaxo, ymaxo, xusc, yusc, xvsc, yvsc;
} gfrezt;
static gfrezt gfrez[6];
typedef struct { char z[2]; long n; } slstype;
void setkv6 (vx5sprite *spr, long dacol)
{
	point3d tp, tp2; float f, det;
	long i, j, k, x, y, z, c, d, x0, y0, z0, x1, y1, z1, xi, yi, zi;
	long xo, yo, zo, xu, yu, zu, xv, yv, zv, stu, stv, tu, tv;
	long xx, yy, xmin, xmax, ymin, ymax, isrhs, ihxi, ihyi, ihzi, syshpit;
	long isx, isy, isz, ihx, ihy, ihz, ifx, ify, ifz, iox, ioy, ioz;
	long sx, sy, sx0, sy0, sz0, rx, ry, rz, pz, dcnt, dcnt2, vismask, xysiz;
	long bx0, by0, bz0, bx1, by1, bz1, *lptr, *shead, shpit, scnt, sstop;
	gfrezt *gf;
	slstype *slst;
	kv6data *kv;
	kv6voxtype *v0, *v1, *v2, *v3;
	void (*modslab)(long *, long, long);

	if (spr->flags&2) return;
	kv = spr->voxnum; if (!kv) return;

		//Calculate top-left-up corner in VXL world coordinates
	tp.x = kv->xpiv + .5;
	tp.y = kv->ypiv + .5;
	tp.z = kv->zpiv + .5;
	tp2.x = spr->p.x - spr->s.x*tp.x - spr->h.x*tp.y - spr->f.x*tp.z;
	tp2.y = spr->p.y - spr->s.y*tp.x - spr->h.y*tp.y - spr->f.y*tp.z;
	tp2.z = spr->p.z - spr->s.z*tp.x - spr->h.z*tp.y - spr->f.z*tp.z;

		//Get bounding x-y box of entire freeze area:
	bx0 = VSID; by0 = VSID; bz0 = MAXZDIM; bx1 = 0; by1 = 0; bz1 = 0;
	for(z=kv->zsiz;z>=0;z-=kv->zsiz)
		for(y=kv->ysiz;y>=0;y-=kv->ysiz)
			for(x=kv->xsiz;x>=0;x-=kv->xsiz)
			{
				ftol(spr->s.x*(float)x + spr->h.x*(float)y + spr->f.x*(float)z + tp2.x,&i);
				if (i < bx0) bx0 = i;
				if (i > bx1) bx1 = i;
				ftol(spr->s.y*(float)x + spr->h.y*(float)y + spr->f.y*(float)z + tp2.y,&i);
				if (i < by0) by0 = i;
				if (i > by1) by1 = i;
				ftol(spr->s.z*(float)x + spr->h.z*(float)y + spr->f.z*(float)z + tp2.z,&i);
				if (i < bz0) bz0 = i;
				if (i > bz1) bz1 = i;
			}
	bx0 -= 2; if (bx0 < 0) bx0 = 0;
	by0 -= 2; if (by0 < 0) by0 = 0;
	bz0 -= 2; if (bz0 < 0) bz0 = 0;
	bx1 += 2; if (bx1 > VSID) bx1 = VSID;
	by1 += 2; if (by1 > VSID) by1 = VSID;
	bz1 += 2; if (bz1 > MAXZDIM) bz1 = MAXZDIM;
	vx5.minx = bx0; vx5.maxx = bx1;
	vx5.miny = by0; vx5.maxy = by1;
	vx5.minz = bz0; vx5.maxz = bz1;

	shpit = bx1-bx0; i = (by1-by0)*shpit*sizeof(shead[0]);
		//Make sure to use array that's big enough: umost is 1MB
	shead = (long *)(((long)umost) - (by0*shpit+bx0)*sizeof(shead[0]));
	slst = (slstype *)(((long)umost)+i);
	scnt = 1; sstop = (sizeof(umost)-i)/sizeof(slstype);
	memset(umost,0,i);

	f = (float)(1<<LSC3);
	ftol(spr->s.x*f,&isx); ftol(spr->s.y*f,&isy); ftol(spr->s.z*f,&isz);
	ftol(spr->h.x*f,&ihx); ftol(spr->h.y*f,&ihy); ftol(spr->h.z*f,&ihz);
	ftol(spr->f.x*f,&ifx); ftol(spr->f.y*f,&ify); ftol(spr->f.z*f,&ifz);
	ftol(tp2.x*f,&iox);
	ftol(tp2.y*f,&ioy);
	ftol(tp2.z*f,&ioz);

		//Determine whether sprite is RHS(1) or LHS(0)
	det = (spr->h.y*spr->f.z - spr->h.z*spr->f.y)*spr->s.x +
			(spr->h.z*spr->f.x - spr->h.x*spr->f.z)*spr->s.y +
			(spr->h.x*spr->f.y - spr->h.y*spr->f.x)*spr->s.z;
	if ((*(long *)&det) > 0) isrhs = 1;
	else if ((*(long *)&det) < 0) isrhs = 0;
	else return;

	xi = (((ifx*ihy-ihx*ify)>>31)|1);
	yi = (((isx*ify-ifx*isy)>>31)|1);
	zi = (((ihx*isy-isx*ihy)>>31)|1);
	if (xi > 0) { x0 = 0; x1 = kv->xsiz; } else { x0 = kv->xsiz-1; x1 = -1; }
	if (yi > 0) { y0 = 0; y1 = kv->ysiz; } else { y0 = kv->ysiz-1; y1 = -1; }
	if (zi > 0) { z0 = 0; z1 = kv->zsiz; } else { z0 = kv->zsiz-1; z1 = -1; }

	vismask = (zi<<3)+24 + (yi<<1)+6 + (xi>>1)+2;

	dcnt = 0;
	for(j=2;j;j--)
	{
		dcnt2 = dcnt;
		vismask = ~vismask;
		for(i=1;i<64;i+=i)
		{
			if (!(vismask&i)) continue;

			if (i&0x15) { xo = yo = zo = 0; }
			else if (i == 2) { xo = isx; yo = isy; zo = isz; }
			else if (i == 8) { xo = ihx; yo = ihy; zo = ihz; }
			else             { xo = ifx; yo = ify; zo = ifz; }

				  if (i&3)  { xu = ihx; yu = ihy; zu = ihz; xv = ifx; yv = ify; zv = ifz; }
			else if (i&12) { xu = isx; yu = isy; zu = isz; xv = ifx; yv = ify; zv = ifz; }
			else           { xu = isx; yu = isy; zu = isz; xv = ihx; yv = ihy; zv = ihz; }

			if ((yu < 0) || ((!yu) && (xu < 0)))
				{ xo += xu; yo += yu; zo += zu; xu = -xu; yu = -yu; zu = -zu; }
			if ((yv < 0) || ((!yv) && (xv < 0)))
				{ xo += xv; yo += yv; zo += zv; xv = -xv; yv = -yv; zv = -zv; }
			d = xv*yu - xu*yv; if (!d) continue;
			if (d < 0)
			{
				k = xu; xu = xv; xv = k;
				k = yu; yu = yv; yv = k;
				k = zu; zu = zv; zv = k; d = -d;
			}
			xmin = ymin = xmax = ymax = 0;
			if (xu < 0) xmin += xu; else xmax += xu;
			if (yu < 0) ymin += yu; else ymax += yu;
			if (xv < 0) xmin += xv; else xmax += xv;
			if (yv < 0) ymin += yv; else ymax += yv;

			gf = &gfrez[dcnt];
			gf->xo = xo; gf->yo = yo; gf->zo = zo;
			gf->xu = xu; gf->yu = yu;
			gf->xv = xv; gf->yv = yv;
			gf->xmino = xmin; gf->ymino = ymin;
			gf->xmaxo = xmax; gf->ymaxo = ymax;
			gf->xusc = (xu<<LSC3); gf->yusc = (yu<<LSC3);
			gf->xvsc = (xv<<LSC3); gf->yvsc = (yv<<LSC3);
			gf->d = d; gf->msk = i;

			f = 1.0 / (float)d;
			ftol(((float)gf->yusc * (float)zv - (float)gf->yvsc * (float)zu) * f,&gf->pzi);
			f *= 4194304.0;
			ftol((float)zu*f,&gf->zu);
			ftol((float)zv*f,&gf->zv);

			dcnt++;
		}
	}

	ihxi = ihx*yi;
	ihyi = ihy*yi;
	ihzi = ihz*yi;

	if (xi < 0) v0 = kv->vox+kv->numvoxs; else v0 = kv->vox;
	for(x=x0;x!=x1;x+=xi)
	{
		i = (long)kv->xlen[x];
		if (xi < 0) v0 -= i;
		if (yi < 0) v1 = v0+i; else v1 = v0;
		if (xi >= 0) v0 += i;
		xysiz = x*kv->ysiz;
		sx0 = isx*x + ihx*y0 + iox;
		sy0 = isy*x + ihy*y0 + ioy;
		sz0 = isz*x + ihz*y0 + ioz;
		for(y=y0;y!=y1;y+=yi)
		{
			i = (long)kv->ylen[xysiz+y];
			if (yi < 0) v1 -= i;
			if (zi < 0) { v2 = v1+i-1; v3 = v1-1; }
					 else { v2 = v1; v3 = v1+i; }
			if (yi >= 0) v1 += i;
			while (v2 != v3)
			{
				z = v2->z; //c = v2->col;
				rx = ifx*z + sx0;
				ry = ify*z + sy0;
				rz = ifz*z + sz0;
				for(i=0;i<dcnt;i++)
				{
					gf = &gfrez[i]; if (!(v2->vis&gf->msk)) continue;
					xo = gf->xo + rx;
					yo = gf->yo + ry;
					zo = gf->zo + rz;
					xmin = ((gf->xmino + xo)>>LSC3); if (xmin < 0) xmin = 0;
					ymin = ((gf->ymino + yo)>>LSC3); if (ymin < 0) ymin = 0;
					xmax = ((gf->xmaxo + xo)>>LSC3)+1; if (xmax > VSID) xmax = VSID;
					ymax = ((gf->ymaxo + yo)>>LSC3)+1; if (ymax > VSID) ymax = VSID;
					xx = (xmin<<LSC3) - xo;
					yy = (ymin<<LSC3) - yo;
					stu = yy*gf->xu - xx*gf->yu;
					stv = yy*gf->xv - xx*gf->yv - gf->d;
					syshpit = ymin*shpit;
					for(sy=ymin;sy<ymax;sy++,syshpit+=shpit)
					{
						tu = stu; stu += gf->xusc;
						tv = stv; stv += gf->xvsc;
						sx = xmin;
						while ((tu&tv) >= 0)
						{
							sx++; if (sx >= xmax) goto freezesprcont;
							tu -= gf->yusc; tv -= gf->yvsc;
						}
						tu = ~tu; tv += gf->d;
						pz = dmulshr22(tu,gf->zv,tv,gf->zu) + zo; j = syshpit+sx;
						tu -= gf->d; tv = ~tv;
						while ((tu&tv) < 0)
						{
							if (i < dcnt2)
							{
								if (scnt >= sstop) return; //OUT OF BUFFER SPACE!
								slst[scnt].z[0] = (char)lbound0(pz>>LSC3,255);
								slst[scnt].n = shead[j]; shead[j] = scnt; scnt++;
							}
							else slst[shead[j]].z[1] = (char)lbound0(pz>>LSC3,255);
							tu += gf->yusc; tv += gf->yvsc; pz += gf->pzi; j++;
						}
freezesprcont:;
					}
				}
				v2 += zi;
			}
			sx0 += ihxi; sy0 += ihyi; sz0 += ihzi;
		}
	}

	if (dacol == -1) modslab = delslab; else modslab = insslab;

	if (vx5.colfunc == kv6colfunc) kv6colfuncinit(spr,det);

	j = by0*shpit+bx0;
	for(sy=by0;sy<by1;sy++)
		for(sx=bx0;sx<bx1;sx++,j++)
		{
			i = shead[j]; if (!i) continue;
			lptr = scum2(sx,sy);
			do
			{
				modslab(lptr,(long)slst[i].z[isrhs],(long)slst[i].z[isrhs^1]);
				i = slst[i].n;
			} while (i);
		}
	scum2finish();
	updatebbox(vx5.minx,vx5.miny,vx5.minz,vx5.maxx,vx5.maxy,vx5.maxz,dacol);
}

#endif

	//Sprite structure is already allocated
	//kv6, vox, xlen, ylen are all malloced in here!
long meltsphere (vx5sprite *spr, lpoint3d *hit, long hitrad)
{
	long i, j, x, y, z, xs, ys, zs, xe, ye, ze, sq, z0, z1;
	long oxvoxs, oyvoxs, numvoxs, cx, cy, cz, cw;
	float f, ff;
	kv6data *kv;
	kv6voxtype *voxptr;
	unsigned long *xlenptr;
	unsigned short *ylenptr;

	xs = max(hit->x-hitrad,0); xe = min(hit->x+hitrad,VSID-1);
	ys = max(hit->y-hitrad,0); ye = min(hit->y+hitrad,VSID-1);
	zs = max(hit->z-hitrad,0); ze = min(hit->z+hitrad,MAXZDIM-1);
	if ((xs > xe) || (ys > ye) || (zs > ze)) return(0);

	if (hitrad >= SETSPHMAXRAD-1) hitrad = SETSPHMAXRAD-2;

	tempfloatbuf[0] = 0.0f;
#if 0
		//Totally unoptimized
	for(i=1;i<=hitrad;i++) tempfloatbuf[i] = pow((float)i,vx5.curpow);
#else
	tempfloatbuf[1] = 1.0f;
	for(i=2;i<=hitrad;i++)
	{
		if (!factr[i][0]) tempfloatbuf[i] = exp(logint[i]*vx5.curpow);
		else tempfloatbuf[i] = tempfloatbuf[factr[i][0]]*tempfloatbuf[factr[i][1]];
	}
#endif
	*(long *)&tempfloatbuf[hitrad+1] = 0x7f7fffff; //3.4028235e38f; //Highest float

// ---------- Need to know how many voxels to allocate... SLOW!!! :( ----------
	cx = cy = cz = 0; //Centroid
	cw = 0;       //Weight (1 unit / voxel)
	numvoxs = 0;
	sq = 0; //pow(fabs(x-hit->x),vx5.curpow) + "y + "z < pow(vx5.currad,vx5.curpow)
	for(x=xs;x<=xe;x++)
	{
		ff = tempfloatbuf[hitrad]-tempfloatbuf[labs(x-hit->x)];
		for(y=ys;y<=ye;y++)
		{
			f = ff-tempfloatbuf[labs(y-hit->y)];
			if (*(long *)&f > 0) //WARNING: make sure to always write ylenptr!
			{
				while (*(long *)&tempfloatbuf[sq] <  *(long *)&f) sq++;
				while (*(long *)&tempfloatbuf[sq] >= *(long *)&f) sq--;
				z0 = max(hit->z-sq,zs); z1 = min(hit->z+sq+1,ze);
				for(z=z0;z<z1;z++)
				{
					i = getcube(x,y,z); //0:air, 1:unexposed solid, 2:vbuf col ptr
					if (i)
					{
						cx += (x-hit->x); cy += (y-hit->y); cz += (z-hit->z); cw++;
					}
					if ((i == 0) || ((i == 1) && (1))) continue; //not_on_border))) continue; //FIX THIS!!!
					numvoxs++;
				}
			}
		}
	}
	if (numvoxs <= 0) return(0); //No voxels found!
// ---------------------------------------------------------------------------

	f = 1.0 / (float)cw; //Make center of sprite the centroid
	spr->p.x = (float)hit->x + (float)cx*f;
	spr->p.y = (float)hit->y + (float)cy*f;
	spr->p.z = (float)hit->z + (float)cz*f;
	spr->s.x = 1.f; spr->h.x = 0.f; spr->f.x = 0.f;
	spr->s.y = 0.f; spr->h.y = 1.f; spr->f.y = 0.f;
	spr->s.z = 0.f; spr->h.z = 0.f; spr->f.z = 1.f;

	x = xe-xs+1; y = ye-ys+1; z = ze-zs+1;

	j = sizeof(kv6data) + numvoxs*sizeof(kv6voxtype) + x*4 + x*y*2;
	i = (long)malloc(j); if (!i) return(0); if (i&3) { free((void *)i); return(0); }
	spr->voxnum = kv = (kv6data *)i; spr->flags = 0;
	kv->leng = j;
	kv->xsiz = x;
	kv->ysiz = y;
	kv->zsiz = z;
	kv->xpiv = spr->p.x - xs;
	kv->ypiv = spr->p.y - ys;
	kv->zpiv = spr->p.z - zs;
	kv->numvoxs = numvoxs;
	kv->namoff = 0;
	kv->lowermip = 0;
	kv->vox = (kv6voxtype *)((long)spr->voxnum+sizeof(kv6data));
	kv->xlen = (unsigned long *)(((long)kv->vox)+numvoxs*sizeof(kv6voxtype));
	kv->ylen = (unsigned short *)(((long)kv->xlen) + kv->xsiz*4);

	voxptr = kv->vox; numvoxs = 0;
	xlenptr = kv->xlen; oxvoxs = 0;
	ylenptr = kv->ylen; oyvoxs = 0;

	sq = 0; //pow(fabs(x-hit->x),vx5.curpow) + "y + "z < pow(vx5.currad,vx5.curpow)
	for(x=xs;x<=xe;x++)
	{
		ff = tempfloatbuf[hitrad]-tempfloatbuf[labs(x-hit->x)];
		for(y=ys;y<=ye;y++)
		{
			f = ff-tempfloatbuf[labs(y-hit->y)];
			if (*(long *)&f > 0) //WARNING: make sure to always write ylenptr!
			{
				while (*(long *)&tempfloatbuf[sq] <  *(long *)&f) sq++;
				while (*(long *)&tempfloatbuf[sq] >= *(long *)&f) sq--;
				z0 = max(hit->z-sq,zs); z1 = min(hit->z+sq+1,ze);
				for(z=z0;z<z1;z++)
				{
					i = getcube(x,y,z); //0:air, 1:unexposed solid, 2:vbuf col ptr
					if ((i == 0) || ((i == 1) && (1))) continue; //not_on_border))) continue; //FIX THIS!!!
					voxptr[numvoxs].col = lightvox(*(long *)i);
					voxptr[numvoxs].z = z-zs;
					voxptr[numvoxs].vis = 63; //FIX THIS!!!
					voxptr[numvoxs].dir = 0; //FIX THIS!!!
					numvoxs++;
				}
			}
			*ylenptr++ = numvoxs-oyvoxs; oyvoxs = numvoxs;
		}
		*xlenptr++ = numvoxs-oxvoxs; oxvoxs = numvoxs;
	}
	return(cw);
}

	//Sprite structure is already allocated
	//kv6, vox, xlen, ylen are all malloced in here!
long meltspans (vx5sprite *spr, vspans *lst, long lstnum, lpoint3d *offs)
{
	float f;
	long i, j, x, y, z, xs, ys, zs, xe, ye, ze, z0, z1;
	long ox, oy, oxvoxs, oyvoxs, numvoxs, cx, cy, cz, cw;
	kv6data *kv;
	kv6voxtype *voxptr;
	unsigned long *xlenptr;
	unsigned short *ylenptr;

	if (lstnum <= 0) return(0);
// ---------- Need to know how many voxels to allocate... SLOW!!! :( ----------
	cx = cy = cz = 0; //Centroid
	cw = 0;       //Weight (1 unit / voxel)
	numvoxs = 0;
	xs = xe = ((long)lst[0].x)+offs->x;
	ys = ((long)lst[       0].y)+offs->y;
	ye = ((long)lst[lstnum-1].y)+offs->y;
	zs = ze = ((long)lst[0].z0)+offs->z;
	for(j=0;j<lstnum;j++)
	{
		x = ((long)lst[j].x)+offs->x;
		y = ((long)lst[j].y)+offs->y; if ((x|y)&(~(VSID-1))) continue;
			  if (x < xs) xs = x;
		else if (x > xe) xe = x;
		z0 = ((long)lst[j].z0)+offs->z;   if (z0 < 0) z0 = 0;
		z1 = ((long)lst[j].z1)+offs->z+1; if (z1 > MAXZDIM) z1 = MAXZDIM;
		if (z0 < zs) zs = z0;
		if (z1 > ze) ze = z1;
		for(z=z0;z<z1;z++) //getcube too SLOW... FIX THIS!!!
		{
			i = getcube(x,y,z); //0:air, 1:unexposed solid, 2:vbuf col ptr
			if (i) { cx += x-offs->x; cy += y-offs->y; cz += z-offs->z; cw++; }
			if (i&~1) numvoxs++;
		}
	}
	if (numvoxs <= 0) return(0); //No voxels found!
// ---------------------------------------------------------------------------

	f = 1.0 / (float)cw; //Make center of sprite the centroid
	spr->p.x = (float)offs->x + (float)cx*f;
	spr->p.y = (float)offs->y + (float)cy*f;
	spr->p.z = (float)offs->z + (float)cz*f;
	spr->x.x = 0.f; spr->y.x = 1.f; spr->z.x = 0.f;
	spr->x.y = 1.f; spr->y.y = 0.f; spr->z.y = 0.f;
	spr->x.z = 0.f; spr->y.z = 0.f; spr->z.z = 1.f;

	x = xe-xs+1; y = ye-ys+1; z = ze-zs;

	j = sizeof(kv6data) + numvoxs*sizeof(kv6voxtype) + y*4 + x*y*2;
	i = (long)malloc(j); if (!i) return(0); if (i&3) { free((void *)i); return(0); }
	spr->voxnum = kv = (kv6data *)i; spr->flags = 0;
	kv->leng = j;
	kv->xsiz = y;
	kv->ysiz = x;
	kv->zsiz = z;
	kv->xpiv = spr->p.y - ys;
	kv->ypiv = spr->p.x - xs;
	kv->zpiv = spr->p.z - zs;
	kv->numvoxs = numvoxs;
	kv->namoff = 0;
	kv->lowermip = 0;
	kv->vox = (kv6voxtype *)((long)spr->voxnum+sizeof(kv6data));
	kv->xlen = (unsigned long *)(((long)kv->vox)+numvoxs*sizeof(kv6voxtype));
	kv->ylen = (unsigned short *)(((long)kv->xlen) + kv->xsiz*4);

	voxptr = kv->vox; numvoxs = 0;
	xlenptr = kv->xlen; oxvoxs = 0;
	ylenptr = kv->ylen; oyvoxs = 0;
	ox = xs; oy = ys;
	for(j=0;j<lstnum;j++)
	{
		x = ((long)lst[j].x)+offs->x;
		y = ((long)lst[j].y)+offs->y; if ((x|y)&(~(VSID-1))) continue;
		while ((ox != x) || (oy != y))
		{
			*ylenptr++ = numvoxs-oyvoxs; oyvoxs = numvoxs; ox++;
			if (ox > xe)
			{
				*xlenptr++ = numvoxs-oxvoxs; oxvoxs = numvoxs;
				ox = xs; oy++;
			}
		}
		z0 = ((long)lst[j].z0)+offs->z;   if (z0 < 0) z0 = 0;
		z1 = ((long)lst[j].z1)+offs->z+1; if (z1 > MAXZDIM) z1 = MAXZDIM;
		for(z=z0;z<z1;z++) //getcube TOO SLOW... FIX THIS!!!
		{
			i = getcube(x,y,z); //0:air, 1:unexposed solid, 2:vbuf col ptr
			if (!(i&~1)) continue;
			voxptr[numvoxs].col = lightvox(*(long *)i);
			voxptr[numvoxs].z = z-zs;

			voxptr[numvoxs].vis = 63; //FIX THIS!!!
			//if (!isvoxelsolid(x-1,y,z)) voxptr[numvoxs].vis |= 1;
			//if (!isvoxelsolid(x+1,y,z)) voxptr[numvoxs].vis |= 2;
			//if (!isvoxelsolid(x,y-1,z)) voxptr[numvoxs].vis |= 4;
			//if (!isvoxelsolid(x,y+1,z)) voxptr[numvoxs].vis |= 8;
			//if (!isvoxelsolid(x,y,z-1)) voxptr[numvoxs].vis |= 16;
			//if (!isvoxelsolid(x,y,z+1)) voxptr[numvoxs].vis |= 32;

			voxptr[numvoxs].dir = 0; //FIX THIS!!!
			numvoxs++;
		}
	}
	while (1)
	{
		*ylenptr++ = numvoxs-oyvoxs; oyvoxs = numvoxs; ox++;
		if (ox > xe)
		{
			*xlenptr++ = numvoxs-oxvoxs; oxvoxs = numvoxs;
			ox = xs; oy++; if (oy > ye) break;
		}
	}
	return(cw);
}

static void setlighting (long x0, long y0, long z0, long x1, long y1, long z1, long lval)
{
	long i, x, y;
	char *v;

	x0 = max(x0,0); x1 = min(x1,VSID);
	y0 = max(y0,0); y1 = min(y1,VSID);
	z0 = max(z0,0); z1 = min(z1,MAXZDIM);

	lval <<= 24;

		//Set 4th byte of colors to full intensity
	for(y=y0;y<y1;y++)
		for(x=x0;x<x1;x++)
		{
			for(v=sptr[y*VSID+x];v[0];v+=v[0]*4)
				for(i=1;i<v[0];i++)
					(*(long *)&v[i<<2]) = (((*(long *)&v[i<<2])&0xffffff)|lval);
			for(i=1;i<=v[2]-v[1]+1;i++)
				(*(long *)&v[i<<2]) = (((*(long *)&v[i<<2])&0xffffff)|lval);
		}
}

	//Updates Lighting, Mip-mapping, and Floating objects list
typedef struct { long x0, y0, z0, x1, y1, z1, csgdel; } bboxtyp;
#define BBOXSIZ 256
static bboxtyp bbox[BBOXSIZ];
static long bboxnum = 0;
void updatevxl ()
{
	long i;

	for(i=bboxnum-1;i>=0;i--)
	{
		if (vx5.lightmode)
			updatelighting(bbox[i].x0,bbox[i].y0,bbox[i].z0,bbox[i].x1,bbox[i].y1,bbox[i].z1);
		if (vx5.vxlmipuse > 1)
			genmipvxl(bbox[i].x0,bbox[i].y0,bbox[i].x1,bbox[i].y1);
		if ((vx5.fallcheck) && (bbox[i].csgdel))
			checkfloatinbox(bbox[i].x0,bbox[i].y0,bbox[i].z0,bbox[i].x1,bbox[i].y1,bbox[i].z1);
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

static long lightlst[MAXLIGHTS];
static float lightsub[MAXLIGHTS];
	//Re-calculates lighting byte #4 of all voxels inside bounding box
void updatelighting (long x0, long y0, long z0, long x1, long y1, long z1)
{
	point3d tp;
	float f, g, h, fx, fy, fz;
	long i, j, x, y, z, sz0, sz1, offs, cstat, lightcnt;
	long x2, y2, x3, y3;
	char *v;

	if (!vx5.lightmode) return;
	xbsox = -17;

	x0 = max(x0-ESTNORMRAD,0); x1 = min(x1+ESTNORMRAD,VSID);
	y0 = max(y0-ESTNORMRAD,0); y1 = min(y1+ESTNORMRAD,VSID);
	z0 = max(z0-ESTNORMRAD,0); z1 = min(z1+ESTNORMRAD,MAXZDIM);

	x2 = x0; y2 = y0;
	x3 = x1; y3 = y1;
	for(y0=y2;y0<y3;y0=y1)
	{
		y1 = min(y0+64,y3);  //"justfly -" (256 lights): +1024:41sec 512:29 256:24 128:22 64:21 32:21 16:21
		for(x0=x2;x0<x3;x0=x1)
		{
			x1 = min(x0+64,x3);


			if (vx5.lightmode == 2)
			{
				lightcnt = 0; //Find which lights are close enough to affect rectangle
				for(i=vx5.numlights-1;i>=0;i--)
				{
					ftol(vx5.lightsrc[i].p.x,&x);
					ftol(vx5.lightsrc[i].p.y,&y);
					ftol(vx5.lightsrc[i].p.z,&z);
					if (x < x0) x -= x0; else if (x > x1) x -= x1; else x = 0;
					if (y < y0) y -= y0; else if (y > y1) y -= y1; else y = 0;
					if (z < z0) z -= z0; else if (z > z1) z -= z1; else z = 0;
					f = vx5.lightsrc[i].r2;
					if ((float)(x*x+y*y+z*z) < f)
					{
						lightlst[lightcnt] = i;
						lightsub[lightcnt] = 1/(sqrt(f)*f);
						lightcnt++;
					}
				}
			}

			for(y=y0;y<y1;y++)
				for(x=x0;x<x1;x++)
				{
					v = sptr[y*VSID+x]; cstat = 0;
					while (1)
					{
						if (!cstat)
						{
							sz0 = ((long)v[1]); sz1 = ((long)v[2])+1; offs = 7-(sz0<<2);
							cstat = 1;
						}
						else
						{
							sz0 = ((long)v[2])-((long)v[1])-((long)v[0])+2;
							if (!v[0]) break; v += v[0]*4;
							sz1 = ((long)v[3]); sz0 += sz1; offs = 3-(sz1<<2);
							cstat = 0;
						}
						if (z0 > sz0) sz0 = z0;
						if (z1 < sz1) sz1 = z1;
						if (vx5.lightmode < 2)
						{
							for(z=sz0;z<sz1;z++)
							{
								estnorm(x,y,z,&tp);
								ftol((tp.y*.5+tp.z)*64.f+103.5f,&i);
								v[(z<<2)+offs] = *(char *)&i;
							}
						}
						else
						{
							for(z=sz0;z<sz1;z++)
							{
								estnorm(x,y,z,&tp);
								f = (tp.y*.5+tp.z)*16+47.5;
								for(i=lightcnt-1;i>=0;i--)
								{
									j = lightlst[i];
									fx = vx5.lightsrc[j].p.x-(float)x;
									fy = vx5.lightsrc[j].p.y-(float)y;
									fz = vx5.lightsrc[j].p.z-(float)z;
									h = tp.x*fx+tp.y*fy+tp.z*fz; if (*(long *)&h >= 0) continue;
									g = fx*fx+fy*fy+fz*fz; if (g >= vx5.lightsrc[j].r2) continue;

										//g = 1.0/(g*sqrt(g))-lightsub[i]; //1.0/g;
									if (cputype&(1<<25))
									{
										_asm
										{
											movss xmm0, g        ;xmm0=g
											rcpss xmm1, xmm0     ;xmm1=1/g
											rsqrtss xmm0, xmm0   ;xmm0=1/sqrt(g)
											mulss xmm1, xmm0     ;xmm1=1/(g*sqrt(g))
											mov eax, i
											subss xmm1, lightsub[eax*4]
											movss g, xmm1
										}
									}
									else
									{
										_asm
										{
											movd mm0, g
											pfrcp mm1, mm0
											pfrsqrt mm0, mm0
											pfmul mm0, mm1
											mov eax, i
											pfsub mm0, lightsub[eax*4]
											movd g, xmm0
											femms
										}
									}
									f -= g*h*vx5.lightsrc[j].sc;
								}
								if (*(long *)&f > 0x437f0000) f = 255; //0x437f0000 is 255.0
								ftol(f,&i);
								v[(z<<2)+offs] = *(char *)&i;
							}
						}
					}
				}
		}
	}
}

//float detection & falling code begins --------------------------------------
//How to use this section of code:
//Step 1: Call checkfloatinbox after every "deleting" set* call
//Step 2: Call dofalls(); at a constant rate in movement code

	//Adds all slabs inside box (inclusive) to "float check" list
void checkfloatinbox (long x0, long y0, long z0, long x1, long y1, long z1)
{
	long x, y;
	char *ov, *v;

	if (flchkcnt >= FLCHKSIZ) return;

		//Make all off-by-1 hacks in other code unnecessary
	x0 = max(x0-1,0); x1 = min(x1+1,VSID);
	y0 = max(y0-1,0); y1 = min(y1+1,VSID);
	z0 = max(z0-1,0); z1 = min(z1+1,MAXZDIM);

		//Add local box's slabs to flchk list - checked in next dofalls()
	for(y=y0;y<y1;y++)
		for(x=x0;x<x1;x++)
		{
			v = sptr[y*VSID+x];
			while (1)
			{
				ov = v; if ((z1 <= v[1]) || (!v[0])) break;
				v += v[0]*4; if (z0 >= v[3]) continue;
				flchk[flchkcnt].x = x;
				flchk[flchkcnt].y = y;
				flchk[flchkcnt].z = ov[1];
				flchkcnt++; if (flchkcnt >= FLCHKSIZ) return;
			}
		}
}

void isnewfloatingadd (long f)
{
	long v = (((f>>(LOGHASHEAD+3))-(f>>3)) & ((1<<LOGHASHEAD)-1));
	vlst[vlstcnt].b = hhead[v]; hhead[v] = vlstcnt;
	vlst[vlstcnt].v = f; vlstcnt++;
}

long isnewfloatingot (long f)
{
	long v = hhead[((f>>(LOGHASHEAD+3))-(f>>3)) & ((1<<LOGHASHEAD)-1)];
	while (1)
	{
		if (v < 0) return(-1);
		if (vlst[v].v == f) return(v);
		v = vlst[v].b;
	}
}

	//removes a & adds b while preserving index; used only by meltfall(...)
	//Must do nothing if 'a' not in hash
void isnewfloatingchg (long a, long b)
{
	long ov, v, i, j;

	i = (((a>>(LOGHASHEAD+3))-(a>>3)) & ((1<<LOGHASHEAD)-1));
	j = (((b>>(LOGHASHEAD+3))-(b>>3)) & ((1<<LOGHASHEAD)-1));

	v = hhead[i]; ov = -1;
	while (v >= 0)
	{
		if (vlst[v].v == a)
		{
			vlst[v].v = b; if (i == j) return;
			if (ov < 0) hhead[i] = vlst[v].b; else vlst[ov].b = vlst[v].b;
			vlst[v].b = hhead[j]; hhead[j] = v;
			return;
		}
		ov = v; v = vlst[v].b;
	}
}

long isnewfloating (flstboxtype *flb)
{
	float f;
	lpoint3d p, cen;
	long i, j, nx, ny, z0, z1, fend, ovlstcnt, mass;
	char *v, *ov;

	p.x = flb->chk.x; p.y = flb->chk.y; p.z = flb->chk.z;
	v = sptr[p.y*VSID+p.x];
	while (1)
	{
		ov = v;
		if ((p.z < v[1]) || (!v[0])) return(0);
		v += v[0]*4;
		if (p.z < v[3]) break;
	}

	if (isnewfloatingot((long)ov) >= 0) return(0);
	ovlstcnt = vlstcnt;
	isnewfloatingadd((long)ov);
	if (vlstcnt >= VLSTSIZ) return(0); //EVIL HACK TO PREVENT CRASH!

		//Init: centroid, mass, bounding box
	cen = p; mass = 1;
	flb->x0 = p.x-1; flb->y0 = p.y-1; flb->z0 = p.z-1;
	flb->x1 = p.x+1; flb->y1 = p.y+1; flb->z1 = p.z+1;

	fend = 0;
	while (1)
	{
		z0 = ov[1];         if (z0 < flb->z0) flb->z0 = z0;
		z1 = ov[ov[0]*4+3]; if (z1 > flb->z1) flb->z1 = z1;

		i = z1-z0;
		cen.x += p.x*i;
		cen.y += p.y*i;
		cen.z += (((z0+z1)*i)>>1); //sum(z0 to z1-1)
		mass += i;

		for(i=0;i<8;i++) //26-connectivity
		{
			switch(i)
			{
				case 0: nx = p.x-1; ny = p.y  ; if (nx < flb->x0) flb->x0 = nx; break;
				case 1: nx = p.x  ; ny = p.y-1; if (ny < flb->y0) flb->y0 = ny; break;
				case 2: nx = p.x+1; ny = p.y  ; if (nx > flb->x1) flb->x1 = nx; break;
				case 3: nx = p.x  ; ny = p.y+1; if (ny > flb->y1) flb->y1 = ny; break;
				case 4: nx = p.x-1; ny = p.y-1; break;
				case 5: nx = p.x+1; ny = p.y-1; break;
				case 6: nx = p.x-1; ny = p.y+1; break;
				case 7: nx = p.x+1; ny = p.y+1; break;
				default: __assume(0); //tells MSVC default can't be reached
			}
			if ((unsigned long)(nx|ny) >= VSID) continue;

			v = sptr[ny*VSID+nx];
			while (1)
			{
				if (!v[0])
				{
					if (v[1] <= z1) return(0);  //This MUST be <=, (not <) !!!
					break;
				}
				ov = v; v += v[0]*4; //NOTE: this is a 'different' ov
				if ((ov[1] > z1) || (z0 > v[3])) continue; //26-connectivity
				j = isnewfloatingot((long)ov);
				if (j < 0)
				{
					isnewfloatingadd((long)ov);
					if (vlstcnt >= VLSTSIZ) return(0); //EVIL HACK TO PREVENT CRASH!
					fstk[fend].x = nx; fstk[fend].y = ny; fstk[fend].z = (long)ov;
					fend++; if (fend >= FSTKSIZ) return(0); //EVIL HACK TO PREVENT CRASH!
					continue;
				}
				if ((unsigned long)j < ovlstcnt) return(0);
			}
		}

		if (!fend)
		{
			flb->i0 = ovlstcnt;
			flb->i1 = vlstcnt;
			flb->mass = mass; f = 1.0 / (float)mass;
			flb->centroid.x = (float)cen.x*f;
			flb->centroid.y = (float)cen.y*f;
			flb->centroid.z = (float)cen.z*f;
			return(1);
		}
		fend--;
		p.x = fstk[fend].x; p.y = fstk[fend].y; ov = (char *)fstk[fend].z;
	}
}

void startfalls ()
{
	long i, z;

		//This allows clear to be MUCH faster when there isn't much falling
	if (vlstcnt < ((1<<LOGHASHEAD)>>1))
	{
		for(i=vlstcnt-1;i>=0;i--)
		{
			z = vlst[i].v;
			hhead[((z>>(LOGHASHEAD+3))-(z>>3)) & ((1<<LOGHASHEAD)-1)] = -1;
		}
	}
	else { for(z=0;z<(1<<LOGHASHEAD);z++) hhead[z] = -1; }

		//Float detection...
		//flstcnt[].i0/i1 tell which parts of vlst are floating
	vlstcnt = 0;

		//Remove any current pieces that are no longer floating
	for(i=vx5.flstnum-1;i>=0;i--)
		if (!isnewfloating(&vx5.flstcnt[i])) //Modifies flstcnt,vlst[],vlstcnt
			vx5.flstcnt[i] = vx5.flstcnt[--vx5.flstnum]; //onground, so delete flstcnt[i]

		//Add new floating pieces (while space is left on flstcnt)
	if (vx5.flstnum < FLPIECES)
		for(i=flchkcnt-1;i>=0;i--)
		{
			vx5.flstcnt[vx5.flstnum].chk = flchk[i];
			if (isnewfloating(&vx5.flstcnt[vx5.flstnum])) //Modifies flstcnt,vlst[],vlstcnt
			{
				vx5.flstcnt[vx5.flstnum].userval = -1; //New piece: let game programmer know
				vx5.flstnum++; if (vx5.flstnum >= FLPIECES) break;
			}
		}
	flchkcnt = 0;
}

	//Call 0 or 1 times (per flstcnt) between startfalls&finishfalls
void dofall (long i)
{
	long j, z;
	char *v;

		//Falling code... call this function once per piece
	vx5.flstcnt[i].chk.z++;
	for(z=vx5.flstcnt[i].i1-1;z>=vx5.flstcnt[i].i0;z--)
	{
		v = (char *)vlst[z].v; v[1]++; v[2]++;
		v = &v[v[0]*4];
		v[3]++;
		if ((v[3] == v[1]) && (vx5.flstcnt[i].i1 >= 0))
		{
			j = isnewfloatingot((long)v);
				//Make sure it's not part of the same floating object
			if ((j < vx5.flstcnt[i].i0) || (j >= vx5.flstcnt[i].i1))
				vx5.flstcnt[i].i1 = -1; //Mark flstcnt[i] for scum2 fixup
		}
	}

	if (vx5.vxlmipuse > 1)
	{
		long x0, y0, x1, y1;
		x0 = max(vx5.flstcnt[i].x0,0); x1 = min(vx5.flstcnt[i].x1+1,VSID);
		y0 = max(vx5.flstcnt[i].y0,0); y1 = min(vx5.flstcnt[i].y1+1,VSID);
		//FIX ME!!!
		//if ((x1 > x0) && (y1 > y0)) genmipvxl(x0,y0,x1,y1); //Don't replace with bbox!
	}
}

	//Sprite structure is already allocated
	//kv6, vox, xlen, ylen are all malloced in here!
long meltfall (vx5sprite *spr, long fi, long delvxl)
{
	long i, j, k, x, y, z, xs, ys, zs, xe, ye, ze;
	long oxvoxs, oyvoxs, numvoxs;
	char *v, *ov, *nv;
	kv6data *kv;
	kv6voxtype *voxptr;
	unsigned long *xlenptr;
	unsigned short *ylenptr;

	if (vx5.flstcnt[fi].i1 < 0) return(0);

	xs = max(vx5.flstcnt[fi].x0,0); xe = min(vx5.flstcnt[fi].x1,VSID-1);
	ys = max(vx5.flstcnt[fi].y0,0); ye = min(vx5.flstcnt[fi].y1,VSID-1);
	zs = max(vx5.flstcnt[fi].z0,0); ze = min(vx5.flstcnt[fi].z1,MAXZDIM-1);
	if ((xs > xe) || (ys > ye) || (zs > ze)) return(0);

		//Need to know how many voxels to allocate... SLOW :(
	numvoxs = vx5.flstcnt[fi].i0-vx5.flstcnt[fi].i1;
	for(i=vx5.flstcnt[fi].i0;i<vx5.flstcnt[fi].i1;i++)
		numvoxs += ((char *)vlst[i].v)[0];
	if (numvoxs <= 0) return(0); //No voxels found!

	spr->p = vx5.flstcnt[fi].centroid;
	spr->s.x = 1.f; spr->h.x = 0.f; spr->f.x = 0.f;
	spr->s.y = 0.f; spr->h.y = 1.f; spr->f.y = 0.f;
	spr->s.z = 0.f; spr->h.z = 0.f; spr->f.z = 1.f;

	x = xe-xs+1; y = ye-ys+1; z = ze-zs+1;

	j = sizeof(kv6data) + numvoxs*sizeof(kv6voxtype) + x*4 + x*y*2;
	i = (long)malloc(j); if (!i) return(0); if (i&3) { free((void *)i); return(0); }
	spr->voxnum = kv = (kv6data *)i; spr->flags = 0;
	kv->leng = j;
	kv->xsiz = x;
	kv->ysiz = y;
	kv->zsiz = z;
	kv->xpiv = spr->p.x - xs;
	kv->ypiv = spr->p.y - ys;
	kv->zpiv = spr->p.z - zs;
	kv->numvoxs = numvoxs;
	kv->namoff = 0;
	kv->lowermip = 0;
	kv->vox = (kv6voxtype *)((long)spr->voxnum+sizeof(kv6data));
	kv->xlen = (unsigned long *)(((long)kv->vox)+numvoxs*sizeof(kv6voxtype));
	kv->ylen = (unsigned short *)(((long)kv->xlen) + kv->xsiz*4);

	voxptr = kv->vox; numvoxs = 0;
	xlenptr = kv->xlen; oxvoxs = 0;
	ylenptr = kv->ylen; oyvoxs = 0;

	for(x=xs;x<=xe;x++)
	{
		for(y=ys;y<=ye;y++)
		{
			for(v=sptr[y*VSID+x];v[0];v=nv)
			{
				nv = v+v[0]*4;

				i = isnewfloatingot((long)v);
				if (((unsigned long)i >= vx5.flstcnt[fi].i1) || (i < vx5.flstcnt[fi].i0))
					continue;

				for(z=v[1];z<=v[2];z++)
				{
					voxptr[numvoxs].col = lightvox(*(long *)&v[((z-v[1])<<2)+4]);
					voxptr[numvoxs].z = z-zs;

					voxptr[numvoxs].vis = 0; //OPTIMIZE THIS!!!
					if (!isvoxelsolid(x-1,y,z)) voxptr[numvoxs].vis |= 1;
					if (!isvoxelsolid(x+1,y,z)) voxptr[numvoxs].vis |= 2;
					if (!isvoxelsolid(x,y-1,z)) voxptr[numvoxs].vis |= 4;
					if (!isvoxelsolid(x,y+1,z)) voxptr[numvoxs].vis |= 8;
					//if (z == v[1]) voxptr[numvoxs].vis |= 16;
					//if (z == nv[3]-1) voxptr[numvoxs].vis |= 32;
					if (!isvoxelsolid(x,y,z-1)) voxptr[numvoxs].vis |= 16;
					if (!isvoxelsolid(x,y,z+1)) voxptr[numvoxs].vis |= 32;

					voxptr[numvoxs].dir = 0; //FIX THIS!!!
					numvoxs++;
				}
				for(z=nv[3]+v[2]-v[1]-v[0]+2;z<nv[3];z++)
				{
					voxptr[numvoxs].col = lightvox(*(long *)&nv[(z-nv[3])<<2]);
					voxptr[numvoxs].z = z-zs;

					voxptr[numvoxs].vis = 0; //OPTIMIZE THIS!!!
					if (!isvoxelsolid(x-1,y,z)) voxptr[numvoxs].vis |= 1;
					if (!isvoxelsolid(x+1,y,z)) voxptr[numvoxs].vis |= 2;
					if (!isvoxelsolid(x,y-1,z)) voxptr[numvoxs].vis |= 4;
					if (!isvoxelsolid(x,y+1,z)) voxptr[numvoxs].vis |= 8;
					//if (z == v[1]) voxptr[numvoxs].vis |= 16;
					//if (z == nv[3]-1) voxptr[numvoxs].vis |= 32;
					if (!isvoxelsolid(x,y,z-1)) voxptr[numvoxs].vis |= 16;
					if (!isvoxelsolid(x,y,z+1)) voxptr[numvoxs].vis |= 32;

					voxptr[numvoxs].dir = 0; //FIX THIS!!!
					numvoxs++;
				}

#if 0
				if (delvxl) //Quick&dirty dealloc from VXL (bad for holes!)
				{
						//invalidate current vptr safely
					isnewfloatingchg((long)v,0);

					k = nv-v; //perform slng(nv) and adjust vlst at same time
					for(ov=nv;ov[0];ov+=ov[0]*4)
						isnewfloatingchg((long)ov,((long)ov)-k);

					j = (long)ov-(long)nv+(ov[2]-ov[1]+1)*4+4;

						//shift end of RLE column up
					v[0] = nv[0]; v[1] = nv[1]; v[2] = nv[2];
					for(i=4;i<j;i+=4) *(long *)&v[i] = *(long *)&nv[i];

						//remove end of RLE column from vbit
					i = ((((long)(&v[i]))-(long)vbuf)>>2); j = (k>>2)+i;
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
					nv = v;
				}
#endif
			}
			*ylenptr++ = numvoxs-oyvoxs; oyvoxs = numvoxs;
		}
		*xlenptr++ = numvoxs-oxvoxs; oxvoxs = numvoxs;
	}

	if (delvxl)
		for(x=xs;x<=xe;x++)
			for(y=ys;y<=ye;y++)
				for(v=sptr[y*VSID+x];v[0];v=nv)
				{
					nv = v+v[0]*4;

					i = isnewfloatingot((long)v);
					if (((unsigned long)i >= vx5.flstcnt[fi].i1) || (i < vx5.flstcnt[fi].i0))
						continue;

						//Quick&dirty dealloc from VXL (bad for holes!)

						//invalidate current vptr safely
					isnewfloatingchg((long)v,0);

					k = nv-v; //perform slng(nv) and adjust vlst at same time
					for(ov=nv;ov[0];ov+=ov[0]*4)
						isnewfloatingchg((long)ov,((long)ov)-k);

					j = (long)ov-(long)nv+(ov[2]-ov[1]+1)*4+4;

						//shift end of RLE column up
					v[0] = nv[0]; v[1] = nv[1]; v[2] = nv[2];
					for(i=4;i<j;i+=4) *(long *)&v[i] = *(long *)&nv[i];

						//remove end of RLE column from vbit
					i = ((((long)(&v[i]))-(long)vbuf)>>2); j = (k>>2)+i;
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
					nv = v;
				}

	vx5.flstcnt[fi].i1 = -2; //Mark flstcnt[i] invalid; no scum2 fixup

	if (vx5.vxlmipuse > 1) genmipvxl(xs,ys,xe+1,ye+1);

	return(vx5.flstcnt[fi].mass);
}

void finishfalls ()
{
	long i, x, y;

		//Scum2 box fixup: refreshes rle voxel data inside a bounding rectangle
	for(i=vx5.flstnum-1;i>=0;i--)
		if (vx5.flstcnt[i].i1 < 0)
		{
			if (vx5.flstcnt[i].i1 == -1)
			{
				for(y=vx5.flstcnt[i].y0;y<=vx5.flstcnt[i].y1;y++)
					for(x=vx5.flstcnt[i].x0;x<=vx5.flstcnt[i].x1;x++)
						scum2(x,y);
				scum2finish();
				updatebbox(vx5.flstcnt[i].x0,vx5.flstcnt[i].y0,vx5.flstcnt[i].z0,vx5.flstcnt[i].x1,vx5.flstcnt[i].y1,vx5.flstcnt[i].z1,0);
			}
			vx5.flstcnt[i] = vx5.flstcnt[--vx5.flstnum]; //onground, so delete flstcnt[i]
		}
}

//float detection & falling code ends ----------------------------------------

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

void freekv6 (kv6data *kv6)
{
	if (kv6->lowermip) freekv6(kv6->lowermip); //NOTE: dangerous - recursive!
	free((void *)kv6);
}

void uninitvoxlap ()
{
	if (sxlbuf) { free(sxlbuf); sxlbuf = 0; }

	if (vbuf) { free(vbuf); vbuf = 0; }
	if (vbit) { free(vbit); vbit = 0; }

	if (khashbuf)
	{     //Free all KV6&KFA on hash list
		long i, j;
		kfatype *kfp;
		for(i=0;i<khashpos;i+=strlen(&khashbuf[i+9])+10)
		{
			switch (khashbuf[i+8])
			{
				case 0: //KV6
					freekv6(*(kv6data **)&khashbuf[i+4]);
					break;
				case 1: //KFA
					kfp = *(kfatype **)&khashbuf[i+4];
					if (!kfp) continue;
					if (kfp->seq) free((void *)kfp->seq);
					if (kfp->frmval) free((void *)kfp->frmval);
					if (kfp->hingesort) free((void *)kfp->hingesort);
					if (kfp->hinge) free((void *)kfp->hinge);
					if (kfp->spr)
					{
						for(j=kfp->numspr-1;j>=0;j--)
							if (kfp->spr[j].voxnum)
								freekv6((kv6data *)kfp->spr[j].voxnum);
						free((void *)kfp->spr);
					}
					free((void *)kfp);
					break;
				default: __assume(0); //tells MSVC default can't be reached
			}
		}
		free(khashbuf); khashbuf = 0; khashpos = khashsiz = 0;
	}

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

		//Setsphere precalculations (factr[] tables) (Derivation in POWCALC.BAS)
		//   if (!factr[z][0]) z's prime else factr[z][0]*factr[z][1] == z
	factr[2][0] = 0; i = 1; j = 9; k = 0;
	for(z=3;z<SETSPHMAXRAD;z+=2)
	{
		if (z == j) { j += (i<<2)+12; i += 2; }
		factr[z][0] = 0; factr[k][1] = z;
		for(zz=3;zz<=i;zz=factr[zz][1])
			if (!(z%zz)) { factr[z][0] = zz; factr[z][1] = z/zz; break; }
		if (!factr[z][0]) k = z;
		factr[z+1][0] = ((z+1)>>1); factr[z+1][1] = 2;
	}
	for(z=1;z<SETSPHMAXRAD;z++) logint[z] = log((double)z);

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

		//Initialize univec normals (for KV6 lighting)
	equivecinit(255);
	//for(i=0;i<255;i++)
	//{
	//   univec[i].z = ((float)((i<<1)-254))/255.0;
	//   f = sqrt(1.0 - univec[i].z*univec[i].z);
	//   fcossin((float)i*(GOLDRAT*PI*2),&univec[i].x,&univec[i].y);
	//   univec[i].x *= f; univec[i].y *= f;
	//}
	//univec[255].x = univec[255].y = univec[255].z = 0;
	for(i=0;i<256;i++)
	{
		iunivec[i][0] = (short)(univec[i].x*4096.0);
		iunivec[i][1] = (short)(univec[i].y*4096.0);
		iunivec[i][2] = (short)(univec[i].z*4096.0);
		iunivec[i][3] = 4096;
	}
	ucossininit();

	memset(mixn,0,sizeof(mixn));

		//Initialize hash table for getkv6()
	memset(khashead,-1,sizeof(khashead));
	if (!(khashbuf = (char *)malloc(KHASHINITSIZE))) return(-1);
	khashsiz = KHASHINITSIZE;

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
	vx5.lightmode = 0;
	vx5.numlights = 0;
	vx5.kv6col = 0x808080;
	vx5.vxlmipuse = 1;
	vx5.fogcol = -1;
	vx5.fallcheck = 0;

	gmipnum = 0;

	return(0);
}

#if 0 //ndef _WIN32
	long i, j, k, l;
	char *v;

	j = k = l = 0;
	for(i=0;i<VSID*VSID;i++)
	{
		for(v=sptr[i];v[0];v+=v[0]*4) { j++; k += v[2]-v[1]+1; l += v[0]-1; }
		k += v[2]-v[1]+1; l += v[2]-v[1]+1;
	}

	printf("VOXLAP5 programmed by Ken Silverman (www.advsys.net/ken)\n");
	printf("Please DO NOT DISTRIBUTE! If this leaks, I will not be happy.\n\n");
	//printf("This copy licensed to:  \n\n");
	printf("Memory statistics upon exit: (all numbers in bytes)");
	printf("\n");
	if (screen) printf("   screen: %8ld\n",imageSize);
	printf("    radar: %8ld\n",max((((MAXXDIM*MAXYDIM*27)>>1)+7)&~7,(VSID+4)*3*SCPITCH*4+8));
	printf("  bacsptr: %8ld\n",sizeof(bacsptr));
	printf("     sptr: %8ld\n",(VSID*VSID)<<2);
	printf("     vbuf: %8ld(%8ld)\n",(j+VSID*VSID+l)<<2,VOXSIZ);
	printf("     vbit: %8ld(%8ld)\n",VOXSIZ>>5);
	printf("\n");
	printf("vbuf head: %8ld\n",(j+VSID*VSID)<<2);
	printf("vbuf cols: %8ld\n",l<<2);
	printf("     fcol: %8ld\n",k<<2);
	printf("     ccol: %8ld\n",(l-k)<<2);
	printf("\n");
	printf("%.2f bytes/column\n",(float)((j+VSID*VSID+l)<<2)/(float)(VSID*VSID));
#endif

// ----- GAME.C code begins

#define HITWALL "wav/quikdrop.wav"
#define BLOWUP "wav/blowup.wav"
#define DEBRIS "wav/debris.wav"
#define DEBRIS2 "wav/plop.wav"
#define BOUNCE "wav/bounce.wav"
#define DIVEBORD "wav/divebord.wav"
#define PICKUP "wav/pickup.wav"
#define DOOROPEN "wav/dooropen.wav"
#define SHOOTBUL "wav/shootgun.wav"
#define MACHINEGUN "wav/shoot.wav"
#define PLAYERHURT "wav/hit.wav"
#define PLAYERDIE "wav/no!.wav"
#define MONSHOOT "wav/monshoot.wav"
#define MONHURT "wav/gothit.wav"
#define MONHURT2 "wav/hurt.wav"
#define MONDIE "wav/mondie.wav"
#define TALK "wav/pop2.wav"
#define WOODRUB "wav/woodrub.wav"
static long volpercent = 100;

char cursxlnam[MAX_PATH+1];
long capturecount = 0, disablemonsts = 1;

long sxleng = 0;

	//Player position variables:
#define CLIPRAD 5
dpoint3d ipos, istr, ihei, ifor, ivel;

	//Tile variables: (info telling where status bars & menus are stored)
typedef struct { long f, p, x, y; } tiletype;
#define FONTXDIM 9
#define FONTYDIM 12
long showtarget = 1;

#define MAXSPRITES 1024 //NOTE: this shouldn't be static!
#ifdef __cplusplus
struct spritetype : vx5sprite //Note: C++!
{
	point3d v, r;  //other attributes (not used by voxlap engine)
	long owner, tim, tag;
};
#else
typedef struct
{
	point3d p; long flags;
	static union { point3d s, x; }; static union { kv6data *voxnum; kfatype *kfaptr; };
	static union { point3d h, y; }; long kfatim;
	static union { point3d f, z; }; long okfatim;
//----------------------------------------------------
	point3d v, r;  //other attributes (not used by voxlap engine)
	long owner, tim, tag;
} spritetype;
#endif

long numsprites, spr2goaltim = 0;
//long sortorder[MAXSPRITES];

long lockanginc = 0;

	//Mouse button state global variables:
long obstatus = 0, bstatus = 0;

	//Timer global variables:
double odtotclk, dtotclk;
float fsynctics;
long totclk;

void findrandomspot (long *x, long *y, long *z)
{
	long cnt;

	cnt = 64;
	do
	{
		(*x) = (rand()&(VSID-1));
		(*y) = (rand()&(VSID-1));
		for((*z)=255;(*z)>=0;(*z)--)
			if (!isvoxelsolid(*x,*y,*z)) break;
		cnt--;
	} while (((*z) < 0) && (cnt > 0));
}

	//Heightmap for bottom
unsigned char bothei[32*32], botheimin;
long botcol[32*32];
long botcolfunc (lpoint3d *p) { return(botcol[(p->y&(32-1))*32+(p->x&(32-1))]); }
void botinit ()
{
	float f, g;
	long i, j, x, y, c;
	float fsc[4+14] = {.4,.45,.55,.45,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
	long col[4+14] = {0x80453010,0x80403510,0x80403015,0x80423112,
							0x80403013,0x80433010,0x80403313,0x80403310,0x80433013,
							0x80403310,0x80403310,0x80433013,0x80433010,0x80403310,
							0x80433010,0x80403013,0x80433310,0x80403013,
							};
	long xpos[4+14] = {0,16, 3,20,  8,24, 0, 7,25,19,12,26,10, 3,29,11,18,25};
	long ypos[4+14] = {0, 5,16,22,  1, 3, 8, 8, 9,13,15,16,12,24,24,28,29,29};
	long rad[4+14] = {68,51,78,67, 20,21,25,19,22,26,18,24,25,23,25,17,23,25};

		//Initialize solid heightmap at bottom of map
	botheimin = 255;
	for(y=0;y<32;y++)
		for(x=0;x<32;x++)
		{
			f = 254.0; c = 0x80403010;
			for(j=0;j<4+14;j++)
			{
				i = (((xpos[j]-x)&31)-16)*(((xpos[j]-x)&31)-16)+
					 (((ypos[j]-y)&31)-16)*(((ypos[j]-y)&31)-16);
				if (i >= rad[j]) continue;
				g = 255.0-sqrt((float)(rad[j]-i))*fsc[j];
				if (g < f) { f = g; c = col[j]; }
			}
			i = y*32+x; botcol[i] = (((rand()<<15)+rand())&0x30303)+c;
			bothei[i] = (unsigned char)f;
			if (bothei[i] < botheimin) botheimin = bothei[i];
		}
}

long initmap ()
{
	kv6data *tempkv6;
	lpoint3d lp;
	float f, g;
	long i, j, k;
	char tempnam[MAX_PATH], *vxlnam, *skynam, *kv6nam, *userst;

		//Initialize solid height&color map at bottom of world
	botinit();

	i = strlen(cursxlnam); numsprites = 0;
	if ((i >= 4) && ((cursxlnam[i-3] == 'S') || (cursxlnam[i-3] == 's')))
	{
		if (loadsxl(cursxlnam,&vxlnam,&skynam,&userst))
		{  //.SXL valid so load sprite info out of .SXL

			if (!loadvxl(vxlnam,&ipos,&istr,&ihei,&ifor))
				printf("failed to load '%s'\n", vxlnam);

		} else
			printf("failed to load '%s'\n", cursxlnam);
	}
	else if ((i >= 4) && ((cursxlnam[i-3] == 'V') || (cursxlnam[i-3] == 'v')))
	{
		if (!loadvxl(cursxlnam,&ipos,&istr,&ihei,&ifor))
			printf("failed to load '%s'\n", cursxlnam);
	}
	else
		printf("failed to load anything\n");

	vx5.lightmode = 1;
	vx5.vxlmipuse = 9; vx5.mipscandist = 192;
	vx5.fallcheck = 1;
	updatevxl();

	vx5.maxscandist = (long)(VSID*1.42);

	//vx5.fogcol = 0x6f6f7f; vx5.maxscandist = 768; //TEMP HACK!!!

	ivel.x = ivel.y = ivel.z = 0;

	return(0);
}

long initapp (long argc, char **argv)
{
	long i, j, k, z, argfilindex, cpuoption = -1;

	prognam = "\"Ken-VOXLAP\" test game";
	xres = 640; yres = 480; colbits = 32; fullscreen = 0; argfilindex = -1;
	for(i=argc-1;i>0;i--)
	{
		if ((argv[i][0] != '/') && (argv[i][0] != '-')) { argfilindex = i; continue; }
		if (!stricmp(&argv[i][1],"win")) { fullscreen = 0; continue; }
		if (!stricmp(&argv[i][1],"3dn")) { cpuoption = 0; continue; }
		if (!stricmp(&argv[i][1],"sse")) { cpuoption = 1; continue; }
		if (!stricmp(&argv[i][1],"sse2")) { cpuoption = 2; continue; }
		//if (!stricmp(&argv[i][1],"?")) { showinfo(); return(-1); }
		if ((argv[i][1] >= '0') && (argv[i][1] <= '9'))
		{
			k = 0; z = 0;
			for(j=1;;j++)
			{
				if ((argv[i][j] >= '0') && (argv[i][j] <= '9'))
					{ k = (k*10+argv[i][j]-48); continue; }
				switch (z)
				{
					case 0: xres = k; break;
					case 1: yres = k; break;
					//case 2: colbits = k; break;
				}
				if (!argv[i][j]) break;
				z++; if (z > 2) break;
				k = 0;
			}
		}
	}
	if (xres > MAXXDIM) xres = MAXXDIM;
	if (yres > MAXYDIM) yres = MAXYDIM;

	if (initvoxlap() < 0) return(-1);

	setsideshades(0,4,1,3,2,2);

		//AthlonXP 2000+: SSE:26.76ms, 3DN:27.34ms, SSE2:28.93ms
	extern long cputype;
	switch(cpuoption)
	{
		case 0: cputype &= ~((1<<25)|(1<<26)); cputype |= ((1<<30)|(1<<31)); break;
		case 1: cputype |= (1<<25); cputype &= ~(1<<26); cputype &= ~((1<<30)|(1<<31)); break;
		case 2: cputype |= ((1<<25)|(1<<26)); cputype &= ~((1<<30)|(1<<31)); break;
		default:;
	}

	if (argfilindex >= 0) strcpy(cursxlnam,argv[argfilindex]);
						  else strcpy(cursxlnam,"vxl/default.sxl");
	if (initmap() < 0) return(-1);

		//Init klock
	readklock(&dtotclk);
	totclk = (long)(dtotclk*1000.0);

	fsynctics = 1.f;

	return(0);
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
			do
			{
				findrandomspot(&lp.x,&lp.y,&lp.z);
				ipos.x = lp.x; ipos.y = lp.y; ipos.z = lp.z-CLIPRAD;
			} while (findmaxcr(ipos.x,ipos.y,ipos.z,CLIPRAD) < CLIPRAD*.9);
			ivel.x = ivel.y = ivel.z = 0;
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
