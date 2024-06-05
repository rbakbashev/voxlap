//VOXLAP engine by Ken Silverman (http://advsys.net/ken)

#include <fcntl.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "sysmain.h"

#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#else
#include <stdarg.h>
#define MAX_PATH 260
#endif

#ifdef __GNUC__
  #include <stdint.h>
  typedef  int16_t i16;
  typedef uint16_t u16;
  typedef  int32_t i32;
  typedef uint32_t u32;
  typedef  int64_t i64;
  typedef uint64_t u64;
#else
  typedef              short i16;
  typedef     unsigned short u16;
  typedef               long i32;
  typedef      unsigned long u32;
  typedef          long long i64;
  typedef unsigned long long u64;
#endif

#if !defined(_WIN32) && !defined(__DOS__)
#include <unistd.h>
#include <dirent.h>
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
typedef struct { double x, y, z; } dpoint3d;

#pragma pack(pop)

	//Voxlap5 shared global variables:
struct
{
		//Opticast variables:
	long anginc, maxscandist;
} vx5;

	//File related functions:
static long loadvxl (const char *, dpoint3d *, dpoint3d *, dpoint3d *, dpoint3d *);

	//Screen related functions:
static void voxsetframebuffer (long, long, long, long);
static void setcamera (dpoint3d *, dpoint3d *, dpoint3d *, dpoint3d *, float, float, float);
static void opticast ();

	//Physics helper functions:
static void dorthorotate (double, double, double, dpoint3d *, dpoint3d *, dpoint3d *);

	//ZIP functions:
static int kzopen (const char *);
static int kzread (void *, int);
static int kzfilelength ();
static int kztell ();
static void kzclose ();

#define PREC (256*4096)
#define CMPPREC (256*4096)

#define VOXSIZ VSID*VSID*128
static char *sptr[(VSID*VSID*4)/3];
static long *vbuf = 0;

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

#pragma pack(push,1)
	//Rendering variables:
typedef struct { long col, dist; } castdat;
typedef struct { castdat *i0, *i1; long z0, z1, cx0, cy0, cx1, cy1; } cftype;
#pragma pack(pop)

static cftype cf[256];

	//Screen related variables:
static long bytesperline, frameplace, xres4;
static long ylookup[MAXYDIM+1];

static lpoint3d glipos;
static point3d gipos, gistr, gihei, gifor;
static point3d gixs, giys, gizs, giadd;
static float gihx, gihy, gihz, gposxfrac[2], gposyfrac[2], grd;
static long gposz, giforzsgn, gstartz0, gstartz1, gixyi[2];
static char *gstartv;

	//Norm flash variables
static long p2c[32], p2m[32];      //bbuf: 2.0K

	//Opticast global variables:
	//radar: 320x200 requires  419560*2 bytes (area * 6.56*2)
	//radar: 400x300 requires  751836*2 bytes (area * 6.27*2)
	//radar: 640x480 requires 1917568*2 bytes (area * 6.24*2)
#define SCPITCH 256
static long *radar = 0, *radarmem = 0;
static long *zbuffermem = 0, zbuffersiz = 0;
static castdat *angstart[MAXXDIM*4], *gscanptr;
#define CMPRECIPSIZ MAXXDIM+32
static float cmprecip[CMPRECIPSIZ], wx0, wy0, wx1, wy1;
static long iwx0, iwy0, iwx1, iwy1;
static point3d gcorn[4];
static long lastx[max(MAXYDIM,VSID)], uurendmem[MAXXDIM*2+8], *uurend;

static long gylookup[512+36]; //256+4+128+4+64+4+...
static long gpz[2], gdz[2], gxmax, gixy[2], gpixy;
static long gmaxscandist;

static long zbufoff;
static long gi0;
static long gi1;

static _inline void dcossin (double a, double *c, double *s)
{
	*c = cos(a);
	*s = sin(a);
}

static _inline void ftol (float f, long *a)
{
	*a = (long)(f + 0.5f);
}

static _inline i32 mulshr16 (i32 a, i32 d)
{
	return ((i64)a * (i64)d) >> 16;
}

static _inline i64 mul64 (i32 a, i32 d)
{
	return (i64)a * (i64)d;
}

static _inline long shldiv16 (long a, long b)
{
	i64 div = ((i64)(a >> 16) << 32) + (i64)(a << 16);

	return div / (i64)b;
}

static _inline long isshldiv16safe (long a, long b)
{
	if (a >= 0)
		a = -a;

	if (b >= 0)
		b = -b;

	return (b - (a >> 14)) >> 31;
}

static _inline long dmulrethigh (long b, long c, long a, long d)
{
	return (u64)(c * (i64)b - d * (i64)a) >> 32;
}

static _inline void clearbuf (void *d, long c, long a)
{
	memset(d, a, c << 2);
}

	//if (a < 0) return(0); else if (a > b) return(b); else return(a);
static _inline long lbound0 (long a, long b) //b MUST be >= 0
{
	if ((unsigned long)a <= b) return(a);
	return((~(a>>31))&b);
}

static _inline long signbit(float f)
{
	u32 l;
	memcpy(&l, &f, 4);
	return (l >> 31);
}

static _inline long signbiti(float f)
{
	i32 l;
	memcpy(&l, &f, 4);
	return (l >> 31);
}

static void gline (long leng, float x0, float y0, float x1, float y1)
{
	i64 q;
	float f, f1, f2, vd0, vd1, vz0, vx1, vy1, vz1;
	long j;
	cftype *c;

	long gx, ogx = 0, gy, ixy, col, dax, day;
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
	if (vd0 < 0) vd0 = 0; //vd0 MUST NOT be negative: bad for asm
	vd1 = f;
	ftol(fabs(f1)*PREC,&gdz[0]);
	ftol(fabs(f2)*PREC,&gdz[1]);

	gixy[0] = (signbiti(vx1)<<3)+4; //=sgn(vx1)*4
	gixy[1] = gixyi[signbit(vy1)]; //=sgn(vy1)*4*VSID
	if (gdz[0] <= 0) { ftol(gposxfrac[signbit(vx1)]*fabs(f1)*PREC,&gpz[0]); if (gpz[0] <= 0) gpz[0] = 0x7fffffff; gdz[0] = 0x7fffffff-gpz[0]; } //Hack for divide overflow
	else ftol(gposxfrac[signbit(vx1)]*(float)gdz[0],&gpz[0]);
	if (gdz[1] <= 0) { ftol(gposyfrac[signbit(vy1)]*fabs(f2)*PREC,&gpz[1]); if (gpz[1] <= 0) gpz[1] = 0x7fffffff; gdz[1] = 0x7fffffff-gpz[1]; } //Hack for divide overflow
	else ftol(gposyfrac[signbit(vy1)]*(float)gdz[1],&gpz[1]);

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
	q = mul64(gdz[0],j); q += (u64)gpz[0];
	if (q < (u64)gxmax)
	{
		gxmax = (long)q;
	}
	if (gixy[1] < 0) j = glipos.y; else j = VSID-1-glipos.y;
	q = mul64(gdz[1],j); q += (u64)gpz[1];
	if (q < (u64)gxmax)
	{
		gxmax = (long)q;
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

static void hline (float x0, float y0, float x1, float y1, long *ix0, long *ix1)
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

static void vline (float x0, float y0, float x1, float y1, long *iy0, long *iy1)
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

static void hrendz (long sx, long sy, long p1, long plc, long incr, long j)
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

static void vrendz (long sx, long sy, long p1, long iplc, long iinc)
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

static void setcamera (dpoint3d *ipo, dpoint3d *ist, dpoint3d *ihe, dpoint3d *ifo,
                       float dahx, float dahy, float dahz)
{
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

static void opticast ()
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

	ftol(cx*65536,&cx16);
	ftol(cy*65536,&cy16);

	ftol((x1-x0)/vx5.anginc,&j);
	if ((fy < 0) && (j > 0)) //(cx,cy),(x0,wy0),(x1,wy0)
	{
		ff = (x1-x0) / (float)j; grd = 1.0f / (wy0-cy);
		gscanptr = (castdat *)radar;
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
				if (p0 < p1) hrendz(p0,sy,p1,u,ui,i);
			}
		}
	}

	ftol((y2-y1)/vx5.anginc,&j);
	if ((gx > 0) && (j > 0)) //(cx,cy),(wx1,y1),(wx1,y2)
	{
		ff = (y2-y1) / (float)j; grd = 1.0f / (wx1-cx);
		gscanptr = (castdat *)radar;
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
				  { for(sy=p0;sy<p1;sy++) vrendz(lastx[sy],sy,xres,lastx[sy],1); }
			else { for(sy=p0;sy<p1;sy++) vrendz(lastx[sy],sy,xres,-lastx[sy],-1); }
		}
	}

	ftol((x2-x3)/vx5.anginc,&j);
	if ((gy > 0) && (j > 0)) //(cx,cy),(x2,wy1),(x3,wy1)
	{
		ff = (x2-x3) / (float)j; grd = 1.0f / (wy1-cy);
		gscanptr = (castdat *)radar;
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
				if (p0 < p1) hrendz(p0,sy,p1,u,ui,i);
			}
		}
	}

	ftol((y3-y0)/vx5.anginc,&j);
	if ((fx < 0) && (j > 0)) //(cx,cy),(wx0,y3),(wx0,y0)
	{
		ff = (y3-y0) / (float)j; grd = 1.0f / (wx0-cx);
		gscanptr = (castdat *)radar;
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
			for(sy=p0;sy<p1;sy++) vrendz(0,sy,lastx[sy]+1,0,giforzsgn);
		}
	}
}

static long loadvxl (const char *lodfilnam, dpoint3d *ipo, dpoint3d *ist, dpoint3d *ihe, dpoint3d *ifo)
{
	long i, fsiz;
	char *v;

	if (!vbuf) { vbuf = (long *)malloc((VOXSIZ>>2)<<2); if (!vbuf) evilquit("vbuf malloc failed"); }

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

	return(1);
}

static void dorthorotate (double ox, double oy, double oz, dpoint3d *ist, dpoint3d *ihe, dpoint3d *ifo)
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

//----------------------------------------------------------------------------

static void voxsetframebuffer (long p, long b, long x, long y)
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

// ----- GAME.C code begins

	//Player position variables:
static dpoint3d ipos, istr, ihei, ifor;

	//Timer global variables:
static double curtime;
static float dt;

static long initmap ()
{
	const char *vxlnam = "vxl/untitled.vxl";

	if (!loadvxl(vxlnam, &ipos, &istr, &ihei, &ifor)) {
		printf("failed to load '%s'\n", vxlnam);
		return -1;
	}

	vx5.maxscandist = (long)(VSID*1.42);

	return 0;
}

long initapp (long argc, char **argv)
{
	if (sizeof(u16) != 2 && sizeof(i16) != 2 && \
			sizeof(u32) != 4 && sizeof(i32) != 4 && \
			sizeof(u64) != 8 && sizeof(i64) != 8) {
		puts("int sizes incorrect");
		exit(1);
	}

	prognam = "\"Ken-VOXLAP\" test game";
	xres = 640; yres = 480; colbits = 32; fullscreen = 0;

	  //WARNING: xres&yres are local to VOXLAP5.C so don't rely on them here!
	if (!(radarmem = (long *)malloc(max((((MAXXDIM*MAXYDIM*27)>>1)+7)&~7,(VSID+4)*3*SCPITCH*4+8))))
		return(-1);
	radar = (long *)((((long)radarmem)+7)&~7);

		//Lookup table to save 1 divide for gline()
	for(long i=1;i<CMPRECIPSIZ;i++) cmprecip[i] = CMPPREC/(float)i;

	for(long z=0;z<32;z++) { p2c[z] = (1<<z); p2m[z] = p2c[z]-1; }

	vx5.anginc = 1; //Higher=faster (1:full,2:half)
	vx5.maxscandist = 256; //must be <= 2047

	if (initmap() < 0) return(-1);

		//Init klock
	readklock(&curtime);

	dt = 1.f;

	return 0;
}

void uninitapp ()
{
	if (vbuf) { free(vbuf); vbuf = 0; }

	if (zbuffermem) { free(zbuffermem); zbuffermem = 0; }
	if (radarmem) { free(radarmem); radarmem = 0; radar = 0; }
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
	readmouse(&fmousx,&fmousy,0,0);

	double oldtime = curtime;
	readklock(&curtime);
	dt = (float)(curtime - oldtime);

		//Rotate player's view
	dp.x = istr.z*.1; dp.y = fmousy*.008; dp.z = fmousx*.008;
	dorthorotate(dp.x,dp.y,dp.z,&istr,&ihei,&ifor);

	dpoint3d ivel = { 0, 0, 0 };

	f = dt*60.0;
	if (keystatus[0x1e]) { ivel.x -= istr.x*f; ivel.y -= istr.y*f; ivel.z -= istr.z*f; } // A
	if (keystatus[0x20]) { ivel.x += istr.x*f; ivel.y += istr.y*f; ivel.z += istr.z*f; } // D
	if (keystatus[0x11]) { ivel.x += ifor.x*f; ivel.y += ifor.y*f; ivel.z += ifor.z*f; } // W
	if (keystatus[0x1f]) { ivel.x -= ifor.x*f; ivel.y -= ifor.y*f; ivel.z -= ifor.z*f; } // S
	if (keystatus[0x12]) { ivel.x -= ihei.x*f; ivel.y -= ihei.y*f; ivel.z -= ihei.z*f; } // E
	if (keystatus[0x10]) { ivel.x += ihei.x*f; ivel.y += ihei.y*f; ivel.z += ihei.z*f; } // Q

	f = dt*60.0;
	ipos.x += ivel.x*f;
	ipos.y += ivel.y*f;
	ipos.z += ivel.z*f;
}

/// ------- KPLIB code begins

typedef struct
{
	FILE *fil;   //0:no file open, !=0:open file (either stand-alone or zip)
	int leng;    //Uncompressed file size (bytes)
	int pos;     //Current uncompressed relative file position (0<=pos<=leng)
	int i;       //For stand-alone/ZIP comptyp#0, this is like "uncomptell"
	             //For ZIP comptyp#8&btype==0 "<64K store", this saves i state
} kzfilestate;
static kzfilestate kzfs;

static int kzopen (const char *filnam)
{
	kzfs.fil = fopen(filnam,"rb");

	if (!kzfs.fil)
	{
		printf("file '%s' not found\n", filnam);
		exit(1);
	}

	kzfs.leng = filelength(_fileno(kzfs.fil));
	kzfs.pos = 0;
	kzfs.i = 0;

	return (int)kzfs.fil;
}

	//returns number of bytes copied
static int kzread (void *buffer, int leng)
{
	int i;

	if ((!kzfs.fil) || (leng <= 0)) return(0);

	if (kzfs.pos != kzfs.i) //Seek only when position changes
		fseek(kzfs.fil,kzfs.pos,SEEK_SET);
	i = min(kzfs.leng-kzfs.pos,leng);
	fread(buffer,i,1,kzfs.fil);
	kzfs.i += i; //kzfs.i is a local copy of ftell(kzfs.fil);

	i = kzfs.pos;
	kzfs.pos += leng; if (kzfs.pos > kzfs.leng) kzfs.pos = kzfs.leng;
	return(kzfs.pos-i);
}

static int kzfilelength ()
{
	if (!kzfs.fil) return(0);
	return(kzfs.leng);
}

static int kztell ()
{
	if (!kzfs.fil) return(-1);
	return(kzfs.pos);
}

static void kzclose ()
{
	if (kzfs.fil) { fclose(kzfs.fil); kzfs.fil = 0; }
}
