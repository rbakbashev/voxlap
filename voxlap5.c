// VOXLAP engine by Ken Silverman (http://advsys.net/ken)

#include <fcntl.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#include "sysmain.h"

#if !defined(max)
#define max(a, b) (((a) > (b)) ? (a) : (b))
#endif
#if !defined(min)
#define min(a, b) (((a) < (b)) ? (a) : (b))
#endif

#define MAXXDIM 1024
#define MAXYDIM 768
#define VSID 1024 // Maximum .VXL dimensions in both x & y direction
#define MAXZDIM 256 // Maximum .VXL dimensions in z direction (height)

#define PREC (256 * 4096)

#define VOXSIZ VSID * VSID * 128

#define SCPITCH 256

typedef uintptr_t usize;

typedef struct { int32_t x, y, z; } lpoint3d;
typedef struct { float x, y, z; } point3d;

typedef struct { uint32_t col; float dist; } castdat;
typedef struct { castdat *i0, *i1; int32_t z0, z1, cx0, cy0, cx1, cy1; } cftype;

// Player position variables:
static point3d ipos, istr, ihei, ifor;

// Timer global variables:
static double curtime;
static float dt = 1.f;

static float optistrx, optistry, optiheix, optiheiy, optiaddx, optiaddy;

// Opticast variables:
static int32_t anginc, maxscandist;

static uint8_t* voxbuf;
static uint8_t* slabptr[(VSID * VSID * 4) / 3];

//                     +--------+--------+--------+--------+
//      voxbuf format: |   0:   |   1:   |   2:   |   3:   |
//+--------------------+--------+--------+--------+--------+
//|      First header: | nextptr|   z1   |   z1c  |  dummy |
//|           Color 1: |    b   |    g   |    r   | intens |
//|           Color 2: |    b   |    g   |    r   | intens |
//|             ...    |    b   |    g   |    r   | intens |
//|           Color n: |    b   |    g   |    r   | intens |
//| Additional header: | nextptr|   z1   |   z1c  |   z0   |
//+--------------------+--------+--------+--------+--------+
//  nextptr: add this # <<2 to index to get to next header (0 if none)
//       z1: z floor (top of floor color list)
//      z1c: z bottom of floor color list MINUS 1! - needed to calculate
//             slab size with slng() and used as a separator for fcol/ccol
//       z0: z ceiling (bottom of ceiling color list)

// Rendering variables:
static cftype cf[256];

// Screen related variables:
static uint32_t* pixels;

static lpoint3d iposl;
static float halfxres, halfyres, halfzres, grd;
static int32_t gstartz0, gstartz1;
static uint8_t* gstartv;

// Opticast global variables:
// radar: 320x200 requires  419560*2 bytes (area * 6.56*2)
// radar: 400x300 requires  751836*2 bytes (area * 6.27*2)
// radar: 640x480 requires 1917568*2 bytes (area * 6.24*2)
static int32_t *radar = 0, *radarmem = 0;
static float* zbuffermem;
static size_t zbuffersiz;
static castdat *angstart[MAXXDIM * 4], *gscanptr;
static float wx0, wy0, wx1, wy1;
static int32_t iwx0, iwy0, iwx1, iwy1;
static point3d gcorn[4];
static int32_t lastx[max(MAXYDIM, VSID)], uurend[MAXXDIM * 2 + 8];

static int32_t gpz[2], gdz[2], gxmax, gixy[2];
static uintptr_t gpixy;
static int32_t gmaxscandist;

static int32_t gi0;
static int32_t gi1;

static inline void dcossin(double a, double* c, double* s)
{
	*c = cos(a);
	*s = sin(a);
}

static inline void ftol(float f, int32_t* a)
{
	*a = (int32_t)(f + 0.5f);
}

static inline int32_t mulshr16(int32_t a, int32_t d)
{
	return ((int64_t)a * (int64_t)d) >> 16;
}

static inline int64_t mul64(int32_t a, int32_t d)
{
	return (int64_t)a * (int64_t)d;
}

static inline int32_t shldiv16(int32_t a, int32_t b)
{
	int64_t div = ((int64_t)(a >> 16) << 32) + (int64_t)(a << 16);

	return div / (int64_t)b;
}

static inline int32_t isshldiv16safe(int32_t a, int32_t b)
{
	if (a >= 0)
		a = -a;

	if (b >= 0)
		b = -b;

	return (b - (a >> 14)) >> 31;
}

static inline int32_t dmulrethigh(int32_t b, int32_t c, int32_t a, int32_t d)
{
	return (uint64_t)(c * (int64_t)b - d * (int64_t)a) >> 32;
}

// if (a < 0) return(0); else if (a > b) return(b); else return(a);
static inline int32_t lbound0(int32_t a, int32_t b) // b MUST be >= 0
{
	if ((uint32_t)a <= b)
		return (a);
	return ((~(a >> 31)) & b);
}

static inline int32_t signbit(float f)
{
	uint32_t l;
	memcpy(&l, &f, 4);
	return (l >> 31);
}

static inline int32_t signbiti(float f)
{
	int32_t l;
	memcpy(&l, &f, 4);
	return (l >> 31);
}

static int drawfwall(uint8_t *v, cftype *c, int32_t ogx)
{
	int32_t col;

	if (v[1] != c->z1) {
		if (v[1] > c->z1)
			c->z1 = v[1];
		else {
			do {
				c->z1--;
				col = *(int32_t*)&v[(c->z1 - v[1]) * 4 + 4];
				while (dmulrethigh(c->z1 * PREC - (int32_t)(ipos.z * PREC), c->cx1, c->cy1, ogx) < 0) {
					c->i1->col = col;
					c->i1--;
					if (c->i0 > c->i1)
						return 1;
					c->cx1 -= gi0;
					c->cy1 -= gi1;
				}
			} while (v[1] != c->z1);
		}
	}

	return 0;
}

static int drawcwall(uint8_t *v, cftype *c, int32_t ogx)
{
	int32_t col;

	if (v[3] != c->z0) {
		if (v[3] < c->z0)
			c->z0 = v[3];
		else {
			do {
				c->z0++;
				col = *(int32_t*)&v[(c->z0 - v[3]) * 4 - 4];
				while (dmulrethigh(c->z0 * PREC - (int32_t)(ipos.z * PREC), c->cx0, c->cy0, ogx) >= 0) {
					c->i0->col = col;
					c->i0++;
					if (c->i0 > c->i1)
						return 1;
					c->cx0 += gi0;
					c->cy0 += gi1;
				}
			} while (v[3] != c->z0);
		}
	}

	return 0;
}

static int drawceil(uint8_t *v, cftype *c, int32_t gx)
{
	while (dmulrethigh(c->z0 * PREC - (int32_t)(ipos.z * PREC), c->cx0, c->cy0, gx) >= 0) {
		c->i0->col = (*(int32_t*)&v[-4]);
		c->i0++;
		if (c->i0 > c->i1)
			return 1;
		c->cx0 += gi0;
		c->cy0 += gi1;
	}

	return 0;
}

static int drawflor(uint8_t *v, cftype *c, int32_t gx)
{
	while (dmulrethigh(c->z1 * PREC - (int32_t)(ipos.z * PREC), c->cx1, c->cy1, gx) < 0) {
		c->i1->col = *(int32_t*)&v[4];
		c->i1--;
		if (c->i0 > c->i1)
			return 1;
		c->cx1 -= gi0;
		c->cy1 -= gi1;
	}

	return 0;
}

static int afterdel(uint8_t **v, cftype **c, cftype *ce, uintptr_t *ixy, int32_t *j, int32_t *ogx, int32_t *gx)
{
	(*c)--;
	if ((*c) < &cf[128]) {
		*ixy += gixy[(*j)];
		gpz[(*j)] += gdz[(*j)];
		*j = (((uint32_t)(gpz[1] - gpz[0])) >> 31);
		*ogx = *gx;
		*gx = gpz[(*j)];

		if (*gx > gxmax)
			return 1;
		*v = (uint8_t*)*(uintptr_t*)*ixy;
		(*c) = ce;
	}

	return 0;
}

static int find_highest_intersecting_slab(uint8_t **v, cftype *c, int32_t ogx)
{
	while (1) {
		if (!(*v)[0])
			return 1;
		if (dmulrethigh(((*v)[2] + 1) * PREC - (int32_t)(ipos.z * PREC), c->cx0, c->cy0, ogx) >= 0)
			break;
		*v += (*v)[0] * 4;
	}

	return 0;
}

static int split_cf(uint8_t *v, cftype **c, int32_t ogx, cftype **ce)
{
	castdat *col;
	int32_t gy, dax, day;

	// If next slab ALSO intersects, split cf!
	gy = (v[v[0] * 4 + 3]) * PREC - (int32_t)(ipos.z * PREC);
	if (dmulrethigh(gy, (*c)->cx1, (*c)->cy1, ogx) < 0) {
		col = (*c)->i1;
		dax = (*c)->cx1;
		day = (*c)->cy1;
		while (dmulrethigh((v[2] + 1) * PREC - (int32_t)(ipos.z * PREC), dax, day, ogx) < 0) {
			col -= 1;
			dax -= gi0;
			day -= gi1;
		}
		(*ce)++;
		if ((*ce) >= &cf[192])
			return 1; // Give it max=64 entries like ASM
		for (cftype *c2 = (*ce); c2 > (*c); c2--)
			c2[0] = c2[-1];
		(*c)[1].i1 = col;
		(*c)->i0 = col + 1;
		(*c)[1].cx1 = dax;
		(*c)->cx0 = dax + gi0;
		(*c)[1].cy1 = day;
		(*c)->cy0 = day + gi1;
		(*c)[1].z1 = (*c)->z0 = v[v[0] * 4 + 3];
		(*c)++;
	}

	return 0;
}

static void clearcol(cftype *ce)
{
	for (cftype *c = ce; c >= &cf[128]; c--)
		while (c->i0 <= c->i1) {
			c->i0->col = 0;
			c->i0++;
		}
}

static int deletez(cftype **ce, cftype *c)
{
	(*ce)--;
	if ((*ce) < &cf[128])
		return 1;
	for (cftype *c2 = c; c2 <= (*ce); c2++)
		c2[0] = c2[1];

	return 0;
}

static void gline(int32_t leng, float x0, float y0, float x1, float y1)
{
	int64_t q;
	float f, f1, f2, vd0, vd1, vz0, vx1, vy1, vz1;
	int32_t j;
	cftype* c;

	int32_t gx, ogx = 0;
	uintptr_t ixy;
	cftype *ce;
	uint8_t* v;

	vd0 = x0 * istr.x + y0 * ihei.x + gcorn[0].x;
	vd1 = x0 * istr.y + y0 * ihei.y + gcorn[0].y;
	vz0 = x0 * istr.z + y0 * ihei.z + gcorn[0].z;
	vx1 = x1 * istr.x + y1 * ihei.x + gcorn[0].x;
	vy1 = x1 * istr.y + y1 * ihei.y + gcorn[0].y;
	vz1 = x1 * istr.z + y1 * ihei.z + gcorn[0].z;

	f = sqrt(vx1 * vx1 + vy1 * vy1);
	f1 = f / vx1;
	f2 = f / vy1;
	if (fabs(vx1) > fabs(vy1))
		vd0 = vd0 * f1;
	else
		vd0 = vd1 * f2;
	if (vd0 < 0)
		vd0 = 0; // vd0 MUST NOT be negative: bad for asm
	vd1 = f;
	ftol(fabs(f1) * PREC, &gdz[0]);
	ftol(fabs(f2) * PREC, &gdz[1]);

	gixy[0] = (signbiti(vx1) << 3) + 4; //=sgn(vx1)*4
	gixy[1] = vy1 < 0 ? -(VSID << 2) : VSID << 2; //=sgn(vy1)*4*VSID

	float posxfrac = vx1 < 0 ? (ipos.x - (float)iposl.x) : 1 - (ipos.x - (float)iposl.x);
	float posyfrac = vy1 < 0 ? (ipos.y - (float)iposl.y) : 1 - (ipos.y - (float)iposl.y);

	if (gdz[0] <= 0) {
		ftol(posxfrac * fabs(f1) * PREC, &gpz[0]);
		if (gpz[0] <= 0)
			gpz[0] = 0x7fffffff;
		gdz[0] = 0x7fffffff - gpz[0];
	} // Hack for divide overflow
	else
		ftol(posxfrac * (float)gdz[0], &gpz[0]);

	if (gdz[1] <= 0) {
		ftol(posyfrac * fabs(f2) * PREC, &gpz[1]);
		if (gpz[1] <= 0)
			gpz[1] = 0x7fffffff;
		gdz[1] = 0x7fffffff - gpz[1];
	} // Hack for divide overflow
	else
		ftol(posyfrac * (float)gdz[1], &gpz[1]);

	c = &cf[128];
	c->i0 = gscanptr;
	c->i1 = &gscanptr[leng];
	c->z0 = gstartz0;
	c->z1 = gstartz1;
	if (ifor.z < 0) {
		ftol((vd1 - vd0) * PREC / leng, &gi0);
		ftol(vd0 * PREC, &c->cx0);
		ftol((vz1 - vz0) * PREC / leng, &gi1);
		ftol(vz0 * PREC, &c->cy0);
	} else {
		ftol((vd0 - vd1) * PREC / leng, &gi0);
		ftol(vd1 * PREC, &c->cx0);
		ftol((vz0 - vz1) * PREC / leng, &gi1);
		ftol(vz1 * PREC, &c->cy0);
	}
	c->cx1 = leng * gi0 + c->cx0;
	c->cy1 = leng * gi1 + c->cy0;

	gxmax = gmaxscandist;

	// Clip borders safely (MUST use integers!) - don't wrap around
	if (gixy[0] < 0)
		j = iposl.x;
	else
		j = VSID - 1 - iposl.x;
	q = mul64(gdz[0], j);
	q += (uint64_t)gpz[0];
	if (q < (uint64_t)gxmax)
		gxmax = (int32_t)q;
	if (gixy[1] < 0)
		j = iposl.y;
	else
		j = VSID - 1 - iposl.y;
	q = mul64(gdz[1], j);
	q += (uint64_t)gpz[1];
	if (q < (uint64_t)gxmax)
		gxmax = (int32_t)q;

	//------------------------------------------------------------------------
	ce = c;
	v = gstartv;
	j = (((uint32_t)(gpz[1] - gpz[0])) >> 31);
	gx = gpz[j];
	ixy = gpixy;

	// 0: none, 1: floor, 2: ceil
	int drawmode = 0;

	if (v == (uint8_t*)*(uintptr_t*)gpixy)
		drawmode = 1;
	else
		drawmode = 2;

	while (1) {
		do {
			if (drawmode == 0) {
				if (drawfwall(v, c, ogx)) {
					if (deletez(&ce, c))
						return;
					break;
				}

				if (v == (uint8_t*)*(uintptr_t*)ixy)
					drawmode = 1;
			}

			if (drawmode == 0) {
				if (drawcwall(v, c, ogx)) {
					if (deletez(&ce, c))
						return;
					break;
				}

				drawmode = 2;
			}

			if (drawmode == 2) {
				if (drawceil(v, c, gx)) {
					if (deletez(&ce, c))
						return;
					break;
				}
				drawmode = 1;
			}

			if (drawmode == 1)
				if (drawflor(v, c, gx)) {
					if (deletez(&ce, c))
						return;
					break;
				}
		} while (0);

		drawmode = 0;

		if (afterdel(&v, &c, ce, &ixy, &j, &ogx, &gx))
			break;

		if (find_highest_intersecting_slab(&v, c, ogx))
			continue;

		if (split_cf(v, &c, ogx, &ce))
			return;
	}

	clearcol(ce);
	return;
}

static void hline(float x0, float y0, float x1, float y1, int32_t* ix0, int32_t* ix1)
{
	float dyx;

	dyx = (y1 - y0) * grd; // grd = 1/(x1-x0)

	if (y0 < wy0)
		ftol((wy0 - y0) / dyx + x0, ix0);
	else if (y0 > wy1)
		ftol((wy1 - y0) / dyx + x0, ix0);
	else
		ftol(x0, ix0);

	if (y1 < wy0)
		ftol((wy0 - y0) / dyx + x0, ix1);
	else if (y1 > wy1)
		ftol((wy1 - y0) / dyx + x0, ix1);
	else
		ftol(x1, ix1);

	if ((*ix0) < iwx0)
		(*ix0) = iwx0;
	if ((*ix0) > iwx1)
		(*ix0) = iwx1; //(*ix1) = min(max(*ix1,wx0),wx1);

	gline(
			labs((*ix1) - (*ix0)),
			(float)(*ix0),
			((*ix0) - x1) * dyx + y1,
			(float)(*ix1),
			((*ix1) - x1) * dyx + y1
		);
}

static void vline(float x0, float y0, float x1, float y1, int32_t* iy0, int32_t* iy1)
{
	float dxy;

	dxy = (x1 - x0) * grd; // grd = 1/(y1-y0)

	if (x0 < wx0)
		ftol((wx0 - x0) / dxy + y0, iy0);
	else if (x0 > wx1)
		ftol((wx1 - x0) / dxy + y0, iy0);
	else
		ftol(y0, iy0);

	if (x1 < wx0)
		ftol((wx0 - x0) / dxy + y0, iy1);
	else if (x1 > wx1)
		ftol((wx1 - x0) / dxy + y0, iy1);
	else
		ftol(y1, iy1);

	if ((*iy0) < iwy0)
		(*iy0) = iwy0;
	if ((*iy0) > iwy1)
		(*iy0) = iwy1;

	gline(
			labs((*iy1) - (*iy0)),
			((*iy0) - y1) * dxy + x1,
			(float)(*iy0),
			((*iy1) - y1) * dxy + x1,
			(float)(*iy1)
		);
}

static void hrendz(int32_t sx, int32_t sy, int32_t x_end, int32_t plc, int32_t incr, int32_t j)
{
	float* zb    = zbuffermem + sy * xres + sx;
	uint32_t* p0 = pixels + sy * xres + sx;
	uint32_t* p1 = pixels + sy * xres + x_end;
	float dirx   = optistrx * (float)sx + optiheix * (float)sy + optiaddx;
	float diry   = optistry * (float)sx + optiheiy * (float)sy + optiaddy;

	do {
		*p0 = angstart[plc >> 16][j].col;
		*zb = angstart[plc >> 16][j].dist / sqrt(dirx * dirx + diry * diry);
		dirx += optistrx;
		diry += optistry;
		plc += incr;
		p0++;
	} while (p0 != p1);
}

static void vrendz(int32_t sx, int32_t sy, int32_t x_end, int32_t iplc, int32_t iinc)
{
	float* zb    = zbuffermem + sy * xres + sx;
	uint32_t* p0 = pixels + sy * xres + sx;
	uint32_t* p1 = pixels + sy * xres + x_end;
	float dirx   = optistrx * (float)sx + optiheix * (float)sy + optiaddx;
	float diry   = optistry * (float)sx + optiheiy * (float)sy + optiaddy;

	while (p0 < p1) {
		*p0 = angstart[uurend[sx] >> 16][iplc].col;
		*zb = angstart[uurend[sx] >> 16][iplc].dist / sqrt(dirx * dirx + diry * diry);
		dirx += optistrx;
		diry += optistry;
		uurend[sx] += uurend[sx + MAXXDIM];
		p0++;
		iplc += iinc;
		sx++;
	}
}

static void setcamera(point3d* ipos, point3d* istr, point3d* ihei, point3d* ifor,
	float dahx, float dahy, float dahz)
{
	halfxres = dahx;
	halfyres = dahy;
	halfzres = dahz;

	gcorn[0].x = -halfxres * istr->x - halfyres * ihei->x + halfzres * ifor->x;
	gcorn[0].y = -halfxres * istr->y - halfyres * ihei->y + halfzres * ifor->y;
	gcorn[0].z = -halfxres * istr->z - halfyres * ihei->z + halfzres * ifor->z;
	gcorn[1].x =  xres * istr->x + gcorn[0].x;
	gcorn[1].y =  xres * istr->y + gcorn[0].y;
	gcorn[1].z =  xres * istr->z + gcorn[0].z;
	gcorn[2].x =  yres * ihei->x + gcorn[1].x;
	gcorn[2].y =  yres * ihei->y + gcorn[1].y;
	gcorn[2].z =  yres * ihei->z + gcorn[1].z;
	gcorn[3].x =  yres * ihei->x + gcorn[0].x;
	gcorn[3].y =  yres * ihei->y + gcorn[0].y;
	gcorn[3].z =  yres * ihei->z + gcorn[0].z;
}

static void casty1(float x0, float x1, float fy, float cx, float cy, float cx16, float cy16)
{
	int32_t j, i, p0, p1, kadd, sy, kmul, ui, u, u1;
	float ff, f;

	ftol((x1 - x0) / anginc, &j);
	if ((fy < 0) && (j > 0)) //(cx,cy),(x0,wy0),(x1,wy0)
	{
		ff = (x1 - x0) / (float)j;
		grd = 1.0f / (wy0 - cy);
		gscanptr = (castdat*)radar;
		for (i = 0, f = x0 + ff * .5f; i < j; f += ff, i++) {
			vline(cx, cy, f, wy0, &p0, &p1);
			if (ifor.z < 0)
				angstart[i] = gscanptr + p0;
			else
				angstart[i] = gscanptr - p1;
			gscanptr += labs(p1 - p0) + 1;
		}

		j <<= 16;
		f = (float)j / ((x1 - x0) * grd);
		ftol((cx - x0) * grd * f, &kadd);
		ftol(cx - .5f, &p1);
		p0 = lbound0(p1 + 1, xres);
		p1 = lbound0(p1, xres);
		ftol(cy - 0.50005f, &sy);
		if (sy >= yres)
			sy = yres - 1;
		ff = (fabs((float)p1 - cx) + 1) * f / 2147483647.0 + cy; // Anti-crash hack
		while ((ff < sy) && (sy >= 0))
			sy--;
		if (sy >= 0) {
			ftol(f, &kmul);
			for (; sy >= 0; sy--)
				if (isshldiv16safe(kmul, (sy << 16) - cy16))
					break; // Anti-crash hack
			if (ifor.z < 0)
				i = -sy;
			else
				i = sy;
			for (; sy >= 0; sy--, i -= (ifor.z < 0 ? -1 : 1)) {
				ui = shldiv16(kmul, (sy << 16) - cy16);
				u = mulshr16((p0 << 16) - cx16, ui) + kadd;
				while ((p0 > 0) && (u >= ui)) {
					u -= ui;
					p0--;
				}
				u1 = (p1 - p0) * ui + u;
				while ((p1 < xres) && (u1 < j)) {
					u1 += ui;
					p1++;
				}
				if (p0 < p1)
					hrendz(p0, sy, p1, u, ui, i);
			}
		}
	}
}

static void castx1(float y1, float y2, float gx, float cx, float cy, float cx16, float cy16)
{
	int32_t j, i, p0, p1, kadd, sx, kmul, ui, u;
	float ff, f;

	ftol((y2 - y1) / anginc, &j);
	if ((gx > 0) && (j > 0)) //(cx,cy),(wx1,y1),(wx1,y2)
	{
		ff = (y2 - y1) / (float)j;
		grd = 1.0f / (wx1 - cx);
		gscanptr = (castdat*)radar;
		for (i = 0, f = y1 + ff * .5f; i < j; f += ff, i++) {
			hline(cx, cy, wx1, f, &p0, &p1);
			if (ifor.z < 0)
				angstart[i] = gscanptr - p0;
			else
				angstart[i] = gscanptr + p1;
			gscanptr += labs(p1 - p0) + 1;
		}

		j <<= 16;
		f = (float)j / ((y2 - y1) * grd);
		ftol((cy - y1) * grd * f, &kadd);
		ftol(cy - .5f, &p1);
		p0 = lbound0(p1 + 1, yres);
		p1 = lbound0(p1, yres);
		ftol(cx + 0.50005f, &sx);
		if (sx < 0)
			sx = 0;
		ff = (fabs((float)p1 - cy) + 1) * f / 2147483647.0 + cx; // Anti-crash hack
		while ((ff > sx) && (sx < xres))
			sx++;
		if (sx < xres) {
			ftol(f, &kmul);
			for (; sx < xres; sx++)
				if (isshldiv16safe(kmul, (sx << 16) - cx16))
					break; // Anti-crash hack
			for (; sx < xres; sx++) {
				ui = shldiv16(kmul, (sx << 16) - cx16);
				u = mulshr16((p0 << 16) - cy16, ui) + kadd;
				while ((p0 > 0) && (u >= ui)) {
					u -= ui;
					lastx[--p0] = sx;
				}
				uurend[sx] = u;
				uurend[sx + MAXXDIM] = ui;
				u += (p1 - p0) * ui;
				while ((p1 < yres) && (u < j)) {
					u += ui;
					lastx[p1++] = sx;
				}
			}
			if (ifor.z < 0)
				for (int32_t sy = p0; sy < p1; sy++)
					vrendz(lastx[sy], sy, xres, lastx[sy], 1);
			else
				for (int32_t sy = p0; sy < p1; sy++)
					vrendz(lastx[sy], sy, xres, -lastx[sy], -1);
		}
	}
}

static void casty2(float x2, float x3, float gy, float cx, float cy, float cx16, float cy16)
{
	int32_t j, i, p0, p1, kadd, sy, kmul, ui, u, u1;
	float ff, f;

	ftol((x2 - x3) / anginc, &j);
	if ((gy > 0) && (j > 0)) //(cx,cy),(x2,wy1),(x3,wy1)
	{
		ff = (x2 - x3) / (float)j;
		grd = 1.0f / (wy1 - cy);
		gscanptr = (castdat*)radar;
		for (i = 0, f = x3 + ff * .5f; i < j; f += ff, i++) {
			vline(cx, cy, f, wy1, &p0, &p1);
			if (ifor.z < 0)
				angstart[i] = gscanptr - p0;
			else
				angstart[i] = gscanptr + p1;
			gscanptr += labs(p1 - p0) + 1;
		}

		j <<= 16;
		f = (float)j / ((x2 - x3) * grd);
		ftol((cx - x3) * grd * f, &kadd);
		ftol(cx - .5f, &p1);
		p0 = lbound0(p1 + 1, xres);
		p1 = lbound0(p1, xres);
		ftol(cy + 0.50005f, &sy);
		if (sy < 0)
			sy = 0;
		ff = (fabs((float)p1 - cx) + 1) * f / 2147483647.0 + cy; // Anti-crash hack
		while ((ff > sy) && (sy < yres))
			sy++;
		if (sy < yres) {
			ftol(f, &kmul);
			for (; sy < yres; sy++)
				if (isshldiv16safe(kmul, (sy << 16) - cy16))
					break; // Anti-crash hack
			if (ifor.z < 0)
				i = sy;
			else
				i = -sy;
			for (; sy < yres; sy++, i -= (ifor.z < 0 ? -1 : 1)) {
				ui = shldiv16(kmul, (sy << 16) - cy16);
				u = mulshr16((p0 << 16) - cx16, ui) + kadd;
				while ((p0 > 0) && (u >= ui)) {
					u -= ui;
					p0--;
				}
				u1 = (p1 - p0) * ui + u;
				while ((p1 < xres) && (u1 < j)) {
					u1 += ui;
					p1++;
				}
				if (p0 < p1)
					hrendz(p0, sy, p1, u, ui, i);
			}
		}
	}
}

static void castx2(float y0, float y3, float fx, float cx, float cy, float cx16, float cy16)
{
	int32_t j, i, p0, p1, kadd, sx, kmul, ui, u;
	float ff, f;

	ftol((y3 - y0) / anginc, &j);
	if ((fx < 0) && (j > 0)) //(cx,cy),(wx0,y3),(wx0,y0)
	{
		ff = (y3 - y0) / (float)j;
		grd = 1.0f / (wx0 - cx);
		gscanptr = (castdat*)radar;
		for (i = 0, f = y0 + ff * .5f; i < j; f += ff, i++) {
			hline(cx, cy, wx0, f, &p0, &p1);
			if (ifor.z < 0)
				angstart[i] = gscanptr + p0;
			else
				angstart[i] = gscanptr - p1;
			gscanptr += labs(p1 - p0) + 1;
		}

		j <<= 16;
		f = (float)j / ((y3 - y0) * grd);
		ftol((cy - y0) * grd * f, &kadd);
		ftol(cy - .5f, &p1);
		p0 = lbound0(p1 + 1, yres);
		p1 = lbound0(p1, yres);
		ftol(cx - 0.50005f, &sx);
		if (sx >= xres)
			sx = xres - 1;
		ff = (fabs((float)p1 - cy) + 1) * f / 2147483647.0 + cx; // Anti-crash hack
		while ((ff < sx) && (sx >= 0))
			sx--;
		if (sx >= 0) {
			ftol(f, &kmul);
			for (; sx >= 0; sx--)
				if (isshldiv16safe(kmul, (sx << 16) - cx16))
					break; // Anti-crash hack
			for (; sx >= 0; sx--) {
				ui = shldiv16(kmul, (sx << 16) - cx16);
				u = mulshr16((p0 << 16) - cy16, ui) + kadd;
				while ((p0 > 0) && (u >= ui)) {
					u -= ui;
					lastx[--p0] = sx;
				}
				uurend[sx] = u;
				uurend[sx + MAXXDIM] = ui;
				u += (p1 - p0) * ui;
				while ((p1 < yres) && (u < j)) {
					u += ui;
					lastx[p1++] = sx;
				}
			}
			for (int32_t sy = p0; sy < p1; sy++)
				vrendz(0, sy, lastx[sy] + 1, 0, (ifor.z < 0 ? -1 : 1));
		}
	}
}

static void opticast()
{
	float f, cx, cy, fx, fy, gx, gy, x0, y0, x1, y1, x2, y2, x3, y3;
	int32_t cx16, cy16;

	iposl.x = (int32_t)ipos.x;
	iposl.y = (int32_t)ipos.y;
	iposl.z = (int32_t)ipos.z;

	gpixy = (uintptr_t)&slabptr[iposl.y * VSID + iposl.x];

	gmaxscandist = min(max(maxscandist, 1), 2047) * PREC;

	gstartv = (uint8_t*)*(uintptr_t*)gpixy;
	if (iposl.z >= gstartv[1]) {
		do {
			if (!gstartv[0])
				return;
			gstartv += gstartv[0] * 4;
		} while (iposl.z >= gstartv[1]);
		if (iposl.z < gstartv[3])
			return;
		gstartz0 = gstartv[3];
	} else
		gstartz0 = 0;
	gstartz1 = gstartv[1];

	if (ifor.z == 0)
		f = 32000;
	else
		f = halfzres / ifor.z;
	f = min(max(f, -32000), 32000);
	cx = istr.z * f + halfxres;
	cy = ihei.z * f + halfyres;

	wx0 = (float)(-(anginc));
	wx1 = (float)(xres - 1 + (anginc));
	wy0 = (float)(-(anginc));
	wy1 = (float)(yres - 1 + (anginc));
	ftol(wx0, &iwx0);
	ftol(wx1, &iwx1);
	ftol(wy0, &iwy0);
	ftol(wy1, &iwy1);

	fx = wx0 - cx;
	fy = wy0 - cy;
	gx = wx1 - cx;
	gy = wy1 - cy;
	x0 = x3 = wx0;
	y0 = y1 = wy0;
	x1 = x2 = wx1;
	y2 = y3 = wy1;
	if (fy < 0) {
		if (fx < 0) {
			f = sqrt(fx * fy);
			x0 = cx - f;
			y0 = cy - f;
		}
		if (gx > 0) {
			f = sqrt(-gx * fy);
			x1 = cx + f;
			y1 = cy - f;
		}
	}
	if (gy > 0) {
		if (gx > 0) {
			f = sqrt(gx * gy);
			x2 = cx + f;
			y2 = cy + f;
		}
		if (fx < 0) {
			f = sqrt(-fx * gy);
			x3 = cx - f;
			y3 = cy + f;
		}
	}
	if (x0 > x1) {
		if (fx < 0)
			y0 = fx / gx * fy + cy;
		else
			y1 = gx / fx * fy + cy;
	}
	if (y1 > y2) {
		if (fy < 0)
			x1 = fy / gy * gx + cx;
		else
			x2 = gy / fy * gx + cx;
	}
	if (x2 < x3) {
		if (fx < 0)
			y3 = fx / gx * gy + cy;
		else
			y2 = gx / fx * gy + cy;
	}
	if (y3 < y0) {
		if (fy < 0)
			x0 = fy / gy * fx + cx;
		else
			x3 = gy / fy * fx + cx;
	}
	// This makes precision errors cause pixels to overwrite rather than omit
	x0 -= .01;
	x1 += .01;
	y1 -= .01;
	y2 += .01;
	x3 -= .01;
	x2 += .01;
	y0 -= .01;
	y3 += .01;

	f = (float)PREC / halfzres;
	optistrx = istr.x * f;
	optiheix = ihei.x * f;
	optiaddx = gcorn[0].x * f;
	optistry = istr.y * f;
	optiheiy = ihei.y * f;
	optiaddy = gcorn[0].y * f;

	ftol(cx * 65536, &cx16);
	ftol(cy * 65536, &cy16);

	casty1(x0, x1, fy, cx, cy, cx16, cy16);
	castx1(y1, y2, gx, cx, cy, cx16, cy16);
	casty2(x2, x3, gy, cx, cy, cx16, cy16);
	castx2(y0, y3, fx, cx, cy, cx16, cy16);
}

static inline int filelength(int h)
{
	struct stat st;
	if (fstat(h, &st) < 0)
		return -1;
	return st.st_size;
}

static int32_t loadvxl(const char* lodfilnam, point3d* ipos, point3d* istr, point3d* ihei, point3d* ifor)
{
	FILE* fil;
	int32_t i, fsiz;
	usize vbyte;
	double posd[3], strd[3], heid[3], ford[3];

	voxbuf = malloc(VOXSIZ);
	if (!voxbuf)
		evilquit("voxbuf malloc failed");

	fil = fopen(lodfilnam, "rb");
	if (!fil) {
		printf("file '%s' not found\n", lodfilnam);
		exit(1);
	}

	fsiz = filelength(fileno(fil));
	int pos = 0;

	pos += 4;
	fread(&i, 4, 1, fil);
	if (i != 0x09072000)
		return (0);

	pos += 4;
	fread(&i, 4, 1, fil);
	if (i != VSID)
		return (0);

	pos += 4;
	fread(&i, 4, 1, fil);
	if (i != VSID)
		return (0);

	fread(posd, 24, 1, fil);
	fread(strd, 24, 1, fil);
	fread(heid, 24, 1, fil);
	fread(ford, 24, 1, fil);
	pos += 24 * 4;

	ipos->x = posd[0];
	ipos->y = posd[1];
	ipos->z = posd[2];
	istr->x = strd[0];
	istr->y = strd[1];
	istr->z = strd[2];
	ihei->x = heid[0];
	ihei->y = heid[1];
	ihei->z = heid[2];
	ifor->x = ford[0];
	ifor->y = ford[1];
	ifor->z = ford[2];

	fread(voxbuf, fsiz - pos, 1, fil);

	fclose(fil);

	vbyte = 0;

	memset(slabptr, 0, sizeof(slabptr));

	for (i = 0; i < VSID * VSID; i++) {
		slabptr[i] = voxbuf + vbyte;

		while (voxbuf[vbyte] != 0) // nextptr
			vbyte += (usize)voxbuf[vbyte] << 2;

		usize z_bot_floor_col_list = (usize)(voxbuf[vbyte + 2]) + 1;
		usize z_top_floor_col_list = (usize)(voxbuf[vbyte + 1]); // floor

		vbyte += (z_bot_floor_col_list - z_top_floor_col_list + 1) << 2;
	}

	return (1);
}

static void dorthorotate(double ox, double oy, double oz, point3d* istr, point3d* ihei, point3d* ifor)
{
	double f, t, dx, dy, dz, rr[9];

	dcossin(ox, &ox, &dx);
	dcossin(oy, &oy, &dy);
	dcossin(oz, &oz, &dz);
	f = ox * oz;
	t = dx * dz;
	rr[0] = t * dy + f;
	rr[7] = -f * dy - t;
	f = ox * dz;
	t = dx * oz;
	rr[1] = -f * dy + t;
	rr[6] = t * dy - f;
	rr[2] = dz * oy;
	rr[3] = -dx * oy;
	rr[4] = ox * oy;
	rr[8] = oz * oy;
	rr[5] = dy;
	ox = istr->x;
	oy = ihei->x;
	oz = ifor->x;
	istr->x = ox * rr[0] + oy * rr[3] + oz * rr[6];
	ihei->x = ox * rr[1] + oy * rr[4] + oz * rr[7];
	ifor->x = ox * rr[2] + oy * rr[5] + oz * rr[8];
	ox = istr->y;
	oy = ihei->y;
	oz = ifor->y;
	istr->y = ox * rr[0] + oy * rr[3] + oz * rr[6];
	ihei->y = ox * rr[1] + oy * rr[4] + oz * rr[7];
	ifor->y = ox * rr[2] + oy * rr[5] + oz * rr[8];
	ox = istr->z;
	oy = ihei->z;
	oz = ifor->z;
	istr->z = ox * rr[0] + oy * rr[3] + oz * rr[6];
	ihei->z = ox * rr[1] + oy * rr[4] + oz * rr[7];
	ifor->z = ox * rr[2] + oy * rr[5] + oz * rr[8];
}

//----------------------------------------------------------------------------

static void voxsetframebuffer(uint32_t* _pixels, int pitch, int x, int y)
{
	pixels = _pixels;

	// This sucks, but it crashes without it
	x = min(x, MAXXDIM);
	y = min(y, MAXYDIM);

	if (x != xres || y != yres) {
		xres = x;
		yres = y;
	}

	// Increase Z buffer size if too small
	size_t newzbufsiz = pitch * yres;
	if (newzbufsiz > zbuffersiz || zbuffermem == 0) {
		if (zbuffermem)
			free(zbuffermem);

		zbuffersiz = newzbufsiz;

		if (!(zbuffermem = malloc(zbuffersiz)))
			evilquit("voxsetframebuffer: allocation too big");
	}
}

static int initmap()
{
	const char* vxlnam = "vxl/untitled.vxl";

	if (!loadvxl(vxlnam, &ipos, &istr, &ihei, &ifor)) {
		printf("failed to load '%s'\n", vxlnam);
		return -1;
	}

	maxscandist = (int32_t)(VSID * 1.42);

	return 0;
}

long initapp(long argc, char** argv)
{
	prognam = "\"Ken-VOXLAP\" test game";
	xres = 640;
	yres = 480;
	colbits = 32;
	fullscreen = 0;

	size_t radarmemsz = max(
			(((MAXXDIM * MAXYDIM * 27) >> 1) + 7) & ~7,
			(VSID + 4) * 3 * SCPITCH * 4 + 8
		);
	if (!(radarmem = malloc(radarmemsz)))
		return (-1);
	radar = (int32_t*)((((uintptr_t)radarmem) + 7) & ~7);

	anginc = 1; // Higher=faster (1:full,2:half)
	maxscandist = 256; // must be <= 2047

	if (initmap() < 0)
		return (-1);

	// Init klock
	readklock(&curtime);

	return 0;
}

void uninitapp()
{
	free(voxbuf);
	free(zbuffermem);
	free(radarmem);
}

void doframe()
{
	point3d dp;
	float f, fmousx, fmousy;
	int pitch, screen_w, screen_h;

	if (startdirectdraw(&pixels, &pitch, &screen_w, &screen_h)) {
		voxsetframebuffer(pixels, pitch, screen_w, screen_h);
		setcamera(&ipos, &istr, &ihei, &ifor, xres * .5, yres * .5, xres * .5);
		opticast();

		stopdirectdraw();
		nextpage();
	}

	// Read keyboard, mouse, and timer
	readkeyboard();
	readmouse(&fmousx, &fmousy, 0, 0);

	double oldtime = curtime;
	readklock(&curtime);
	dt = (float)(curtime - oldtime);

	// Rotate player's view
	dp.x = istr.z * .1;
	dp.y = fmousy * .008;
	dp.z = fmousx * .008;
	dorthorotate(dp.x, dp.y, dp.z, &istr, &ihei, &ifor);

	point3d ivel = { 0, 0, 0 };

	f = dt * 60.0;
	if (keystatus[0x1e]) {
		ivel.x -= istr.x * f;
		ivel.y -= istr.y * f;
		ivel.z -= istr.z * f;
	} // A
	if (keystatus[0x20]) {
		ivel.x += istr.x * f;
		ivel.y += istr.y * f;
		ivel.z += istr.z * f;
	} // D
	if (keystatus[0x11]) {
		ivel.x += ifor.x * f;
		ivel.y += ifor.y * f;
		ivel.z += ifor.z * f;
	} // W
	if (keystatus[0x1f]) {
		ivel.x -= ifor.x * f;
		ivel.y -= ifor.y * f;
		ivel.z -= ifor.z * f;
	} // S
	if (keystatus[0x12]) {
		ivel.x -= ihei.x * f;
		ivel.y -= ihei.y * f;
		ivel.z -= ihei.z * f;
	} // E
	if (keystatus[0x10]) {
		ivel.x += ihei.x * f;
		ivel.y += ihei.y * f;
		ivel.z += ihei.z * f;
	} // Q

	f = dt * 60.0;
	ipos.x += ivel.x * f;
	ipos.y += ivel.y * f;
	ipos.z += ivel.z * f;
}
