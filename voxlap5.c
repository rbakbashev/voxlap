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

#define MAXXDIM 1024
#define MAXYDIM 768
#define VSID 1024   //Maximum .VXL dimensions in both x & y direction
#define MAXZDIM 256 //Maximum .VXL dimensions in z direction (height)

#pragma pack(push,1)

typedef struct { long x, y, z; } lpoint3d;
typedef struct { float x, y, z; } point3d;
typedef struct { float x, y, z, z2; } point4d;
typedef struct { double x, y, z; } dpoint3d;

#pragma pack(pop)

#define MAXFRM 1024 //MUST be even number for alignment!

	//Voxlap5 shared global variables:
struct
{
		//Opticast variables:
	long anginc, sideshademode, mipscandist, maxscandist, vxlmipuse;
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

	//VXL MISC functions:
extern void updatebbox (long, long, long, long, long, long, long);
extern void updatevxl ();
extern void genmipvxl (long, long, long, long);

	//ZIP functions:
extern int kzopen (const char *);
extern int kzread (void *, int);
extern int kzfilelength ();
extern int kztell ();
extern void kzclose ();

#include "sysmain.h"


//VOXLAP engine by Ken Silverman (http://advsys.net/ken)

#define PREC (256*4096)
#define CMPPREC (256*4096)
#define FPREC (256*4096)
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

static char nullst = 0; //nullst always NULL string

#pragma pack(push,1)
	//Rendering variables:
typedef struct { long col, dist; } castdat;
typedef struct { castdat *i0, *i1; long z0, z1, cx0, cy0, cx1, cy1; } cftype;
typedef struct { unsigned short x, y; } uspoint2d;
typedef struct { long x, y; } lpoint2d;
typedef struct { float x, y; } point2d;
#pragma pack(pop)

cftype cf[256];

	//Screen related variables:
static long xres, yres, bytesperline, frameplace, xres4;
long ylookup[MAXYDIM+1];

static lpoint3d glipos;
static point3d gipos, gistr, gihei, gifor;
static point3d gixs, giys, gizs, giadd;
static float gihx, gihy, gihz, gposxfrac[2], gposyfrac[2], grd;
static long gposz, giforzsgn, gstartz0, gstartz1, gixyi[2];
static char *gstartv;

	//Norm flash variables
#define GSIZ 512  //NOTE: GSIZ should be 1<<x, and must be <= 65536
static long bbuf[GSIZ][GSIZ>>5], p2c[32], p2m[32];      //bbuf: 2.0K

	//Opticast global variables:
	//radar: 320x200 requires  419560*2 bytes (area * 6.56*2)
	//radar: 400x300 requires  751836*2 bytes (area * 6.27*2)
	//radar: 640x480 requires 1917568*2 bytes (area * 6.24*2)
#define SCPITCH 256
long *radar = 0, *radarmem = 0;
static long *zbuffermem = 0, zbuffersiz = 0;
static castdat *angstart[MAXXDIM*4], *gscanptr;
#define CMPRECIPSIZ MAXXDIM+32
static float cmprecip[CMPRECIPSIZ], wx0, wy0, wx1, wy1;
static long iwx0, iwy0, iwx1, iwy1;
static point3d gcorn[4];
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

long zbufoff;
#ifdef __cplusplus
}
#endif
#define gi0 (((long *)&gi)[0])
#define gi1 (((long *)&gi)[1])

#ifdef _MSC_VER

#pragma warning(disable:4799) //I know how to use EMMS

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

void gline (long leng, float x0, float y0, float x1, float y1)
{
	unsigned __int64 q;
	float f, f1, f2, vd0, vd1, vz0, vx1, vy1, vz1;
	long j;
	cftype *c;

	long gx, ogx, gy, ixy, col, dax, day;
	cftype *c2, *ce;
	char *v;

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

#ifdef _MSC_VER

static __declspec(align(16)) point4d opti4[5];
static __declspec(align(16)) void* opti4asm = opti4;

void (*hrend)(long,long,long,long,long,long);
void (*vrend)(long,long,long,long,long);

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
	for(i=0;i<256+4;i++) gylookup[i] = (i*PREC-gposz);
	gmaxscandist = min(max(vx5.maxscandist,1),2047)*PREC;

	hrend = hrendzsse;
	vrend = vrendzsse;
	nskypic = skyoff = 0;

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

	gmipnum = 1;
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
}

//------------------------- SXL parsing code begins --------------------------

static char *sxlbuf = 0;
static long sxlparspos, sxlparslen;

long loadsxl (const char *sxlnam, char **vxlnam, char **skynam, char **globst)
{
	long j, k, m, n;

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

	if ((b != ylookup[1]) || (x != xres) || (y != yres))
	{
		bytesperline = b; xres = x; yres = y; xres4 = (xres<<2);
		ylookup[0] = 0; for(i=0;i<yres;i++) ylookup[i+1] = ylookup[i]+bytesperline;
		//gihx = gihz = (float)xres*.5f; gihy = (float)yres*.5f; //BAD!!!
		if ((ylookup[yres]+256 > zbuffersiz) || (!zbuffermem))  //Increase Z buffer size if too small
		{
			if (zbuffermem) { free(zbuffermem); zbuffermem = 0; }
			zbuffersiz = ylookup[yres]+256;
			if (!(zbuffermem = (long *)malloc(zbuffersiz))) evilquit("voxsetframebuffer: allocation too big");
		}
	}
		//zbuffer aligns its memory to the same pixel boundaries as the screen!
		//WARNING: Pentium 4's L2 cache has severe slowdowns when 65536-64 <= (zbufoff&65535) < 64
	zbufoff = (((((long)zbuffermem)-frameplace-128)+255)&~255)+128;
	uurend = &uurendmem[((frameplace&4)^(((long)uurendmem)&4))>>2];
}

void uninitvoxlap ()
{
	if (sxlbuf) { free(sxlbuf); sxlbuf = 0; }

	if (vbuf) { free(vbuf); vbuf = 0; }
	if (vbit) { free(vbit); vbit = 0; }

	if (skylng) { free((void *)skylng); skylng = 0; }
	if (skypic) { free((void *)skypic); skypic = skyoff = 0; }

	if (zbuffermem) { free(zbuffermem); zbuffermem = 0; }
	if (radarmem) { free(radarmem); radarmem = 0; radar = 0; }
}

long initvoxlap ()
{
	__int64 q;
	long i, j, k, z, zz;
	float f, ff;

	  //WARNING: xres&yres are local to VOXLAP5.C so don't rely on them here!
	if (!(radarmem = (long *)malloc(max((((MAXXDIM*MAXYDIM*27)>>1)+7)&~7,(VSID+4)*3*SCPITCH*4+8))))
		return(-1);
	radar = (long *)((((long)radarmem)+7)&~7);

		//Lookup table to save 1 divide for gline()
	for(i=1;i<CMPRECIPSIZ;i++) cmprecip[i] = CMPPREC/(float)i;

	for(z=0;z<32;z++) { p2c[z] = (1<<z); p2m[z] = p2c[z]-1; }

	memset(mixn,0,sizeof(mixn));

	vx5.anginc = 1; //Higher=faster (1:full,2:half)
	vx5.sideshademode = 0; setsideshades(0,0,0,0,0,0);
	vx5.mipscandist = 128;
	vx5.maxscandist = 256; //must be <= 2047
	vx5.vxlmipuse = 1;

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
	dpoint3d dp;
	float f, fmousx, fmousy;
	long i, j, k, l;

	if (!startdirectdraw(&i,&j,&k,&l)) goto skipalldraw;

	voxsetframebuffer(i,j,k,l);
	setcamera(&ipos,&istr,&ihei,&ifor,xres*.5,yres*.5,xres*.5);
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

	ivel.x = ivel.y = ivel.z = 0;

	f = fsynctics*60.0;
	if (keystatus[0x1e]) { ivel.x -= istr.x*f; ivel.y -= istr.y*f; ivel.z -= istr.z*f; } // A
	if (keystatus[0x20]) { ivel.x += istr.x*f; ivel.y += istr.y*f; ivel.z += istr.z*f; } // D
	if (keystatus[0x11]) { ivel.x += ifor.x*f; ivel.y += ifor.y*f; ivel.z += ifor.z*f; } // W
	if (keystatus[0x1f]) { ivel.x -= ifor.x*f; ivel.y -= ifor.y*f; ivel.z -= ifor.z*f; } // S
	if (keystatus[0x12]) { ivel.x -= ihei.x*f; ivel.y -= ihei.y*f; ivel.z -= ihei.z*f; } // E
	if (keystatus[0x10]) { ivel.x += ihei.x*f; ivel.y += ihei.y*f; ivel.z += ihei.z*f; } // Q

	f = fsynctics*60.0;
	ipos.x += ivel.x*f;
	ipos.y += ivel.y*f;
	ipos.z += ivel.z*f;

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
