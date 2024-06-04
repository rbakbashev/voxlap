#if defined(_MSC_VER)
  #define SDL_MAIN_HANDLED
  #include "SDL2\SDL.h"
#else
  #include <SDL2/SDL.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#undef _WIN32
#include "sysmain.h"

#define evalmacro(x) evalmacrox(x)
#define evalmacrox(x) #x

char *prognam = "WinMain App (by KS/TD)";
long progresiz = 1;

extern long initapp (long argc, char **argv);
extern void uninitapp ();
extern void doframe ();

static long quitprogram = 0, quitparam;
void breath();

long ActiveApp = 1, alwaysactive = 0;
long xres = 640, yres = 480, colbits = 8, fullscreen = 1, maxpages = 8;

void initklock ()
{
}

void readklock (double *tim)
{
	*tim = (double)(SDL_GetTicks64() / 1000.);
}

//SDL Video VARIABLES & CODE--------------------------------------------------

PALETTEENTRY pal[256];

static SDL_Window *screen = NULL;
static SDL_Renderer *renderer = NULL;
static SDL_Texture *texture = NULL;
static Uint32 *pixels = NULL;
static int pitch = 0;
static int surflocked = 0;

#define MAXVALIDMODES 256
static validmodetype validmodelist[MAXVALIDMODES];
static long validmodecnt = 0;
validmodetype curvidmodeinfo;


static void gvmladd(long x, long y, char c)
{
	if (validmodecnt == MAXVALIDMODES) return;
	memset(&validmodelist[validmodecnt], 0, sizeof(validmodetype));
	validmodelist[validmodecnt].x = x;
	validmodelist[validmodecnt].y = y;
	validmodelist[validmodecnt].c = c;
	validmodecnt++;
}

long getvalidmodelist (validmodetype **davalidmodelist)
{
	return(validmodecnt);
}

void updatepalette (long start, long danum)
{
}

void stopdirectdraw()
{
	if (!surflocked) return;
	SDL_UnlockTexture(texture);
	surflocked = 0;
}

long startdirectdraw(long *vidplc, long *dabpl, long *daxres, long *dayres)
{
	if (surflocked) stopdirectdraw();

	if (SDL_LockTexture(texture, 0, (void**)&pixels, &pitch) < 0) return 0;
	
	*vidplc = (long)pixels; *dabpl = pitch;
	*daxres = xres; *dayres = yres; surflocked = 1;
	return 1;
}

void nextpage()
{
	if (!screen) return;
	if (surflocked) stopdirectdraw();

	SDL_RenderClear(renderer);
	SDL_RenderCopy(renderer, texture, NULL, NULL);
	SDL_RenderPresent(renderer);
}

long clearscreen(long fillcolor)
{
	SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
	SDL_RenderClear(renderer);
	return 1;
}

long initdirectdraw(long daxres, long dayres, long dacolbits)
{
	Uint32 winbits;

	extern long mouse_acquire;
	if (mouse_acquire) {
		SDL_SetRelativeMouseMode(SDL_FALSE);
	}

	xres = daxres; yres = dayres; colbits = dacolbits;

	winbits = 0;
	if (fullscreen) winbits |= SDL_WINDOW_FULLSCREEN;
	if (progresiz) winbits |= SDL_WINDOW_RESIZABLE;

	screen = SDL_CreateWindow(
			prognam,
			SDL_WINDOWPOS_UNDEFINED,
			SDL_WINDOWPOS_UNDEFINED,
			daxres, dayres,
			winbits);
	if (!screen) {
		puts("video init failed!");
		return 0;
	}

	renderer = SDL_CreateRenderer(screen, -1, 0);
	if (!renderer) {
		puts("failed to create renderer");
		return 0;
	}

	texture = SDL_CreateTexture(renderer,
			SDL_PIXELFORMAT_RGB888,
			SDL_TEXTUREACCESS_STREAMING,
			daxres, dayres);
	if (!texture) {
		puts("failed to create texture");
		return 0;
	}

	memset(&curvidmodeinfo, 0, sizeof(validmodetype));

	if (colbits == 8) updatepalette(0,256);

	if (mouse_acquire) SDL_SetRelativeMouseMode(SDL_TRUE);

	return 1;
}

long changeres(long daxres, long dayres, long dacolbits, long dafullscreen)
{
	fullscreen = dafullscreen;
	return initdirectdraw(daxres, dayres, dacolbits);
}

//SDL Input VARIABLES & CODE--------------------------------------------------
char keystatus[256];
#define KEYBUFSIZ 256
static long keybuf[KEYBUFSIZ], keybufr = 0, keybufw = 0, keybufw2 = 0;

char ext_keystatus[256]; // +TD
char ext_mbstatus[8] = {0}; // +TD extended mouse button status
long mouse_acquire = 1, kbd_acquire = 1;
static void (*setmousein)(long, long) = NULL;
static long mousex = 0, mousey = 0, gbstatus = 0;
long mousmoth=0;
float mousper=0.0;

long readkeyboard ()
{
	// all gets done in sdlmsgloop now
	return 0;
}

void readmouse (float *fmousx, float *fmousy, float *fmousz, long *bstatus)
{
	if (!mouse_acquire) { *fmousx = *fmousy = 0; *bstatus = 0; return; }
	
	*fmousx = (float)mousex;
	*fmousy = (float)mousey;
	if (fmousz)
		*fmousz = 0.f;
	*bstatus = gbstatus;

	mousex = mousey = 0;
}

//Quitting routines ----------------------------------------------------------

void quitloop ()
{
	SDL_Event ev;
	ev.type = SDL_QUIT;
	SDL_PushEvent(&ev);
}

void setvolume (long) {}
void playsound (const char *, long, float, void *, long) {}
void setears3d (float, float, float, float, float, float, float, float, float) {}
void playsoundupdate (void *, void *) {}
void kensoundclose() {}

void evilquit (const char *str) //Evil because this function makes awful assumptions!!!
{
	kensoundclose();
	SDL_DestroyTexture(texture);
	SDL_DestroyRenderer(renderer);
	SDL_DestroyWindow(screen);
	SDL_Quit();
	if (str) fprintf(stderr,"fatal error! %s\n",str);
	uninitapp();
	exit(0);
}

//GENERAL CODE-------------------------------------------------------------------

void setacquire (long mouse, long kbd)
{
	// SDL doesn't let the mouse and keyboard be grabbed independantly
	if (mouse || kbd) {
		SDL_SetRelativeMouseMode(SDL_TRUE);
		mouse_acquire = 1;
	} else {
		SDL_SetRelativeMouseMode(SDL_FALSE);
		mouse_acquire = 0;
	}
}

void setmouseout (void (*in)(long,long), long x, long y)
{
	if (fullscreen) return;
	setmousein = in;
	setacquire(0, kbd_acquire);

	SDL_WarpMouseInWindow(screen, x, y);
}

static unsigned char keytranslation[] = {
	0, 0, 0, 0, 0, 0, 0, 0, 14, 15, 0, 0, 0, 28, 0, 0, 0, 0, 0, 89, 0, 0, 0, 
	0, 0, 0, 0, 1, 0, 0, 0, 0, 57, 2, 40, 4, 5, 6, 8, 40, 10, 11, 9, 13, 51,
	12, 52, 53, 11, 2, 3, 4, 5, 6, 7, 8, 9, 10, 39, 39, 51, 13, 52, 53, 3, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	26, 43, 27, 7, 12, 41, 30, 48, 46, 32, 18, 33, 34, 35, 23, 36, 37, 38, 50,
	49, 24, 25, 16, 19, 31, 20, 22, 47, 17, 45, 21, 44, 0, 0, 0, 0, 211, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 82, 79, 80, 81, 75, 76, 77, 71, 72, 73, 83, 181, 55, 74, 78, 156, 0,
	200, 208, 205, 203, 210, 199, 207, 201, 209, 59, 60, 61, 62, 63, 64, 65,
	66, 67, 68, 87, 88, 0, 0, 0, 0, 0, 0, 69, 58, 70, 54, 42, 157, 29, 184, 56,
	0, 0, 219, 220, 0, 0, 0, /*-2*/0, 84, 183, 221, 0, 0, 0
};


static void sdlmsgloop(void)
{
	SDL_Event ev;
	unsigned char sc;
	long shkeystatus;

	while (SDL_PollEvent(&ev)) {
		switch (ev.type) {
			case SDL_KEYDOWN:
				if (ev.key.keysym.sym < 324) {
					sc = keytranslation[ev.key.keysym.sym];
					keystatus[ sc ] = 1;   // FIXME: verify this is kosher
					ext_keystatus[ sc ] = 1|2;
				}
				switch (ev.key.keysym.sym) {
					case SDLK_ESCAPE: quitprogram = 1; break;
					case SDLK_LSHIFT: shkeystatus |= (1<<16); break;
					case SDLK_RSHIFT: shkeystatus |= (1<<17); break;
					case SDLK_LCTRL:  shkeystatus |= (1<<18); break;
					case SDLK_RCTRL:  shkeystatus |= (1<<19); break;
					case SDLK_LALT:   shkeystatus |= (1<<20); break;
					case SDLK_RALT:   shkeystatus |= (1<<21); break;
					default:
						{
							long i = ((keybufw2+1)&(KEYBUFSIZ-1));
							keybuf[keybufw2&(KEYBUFSIZ-1)] = ((long)sc << 8)+shkeystatus;   // FIXME: verify
							if (i != keybufr) keybufw2 = i; //prevent fifo overlap
						}
						break;
				}
				break;
			case SDL_KEYUP:
				if (ev.key.keysym.sym < 324) {
					sc = keytranslation[ev.key.keysym.sym];
					keystatus[ sc ] = 0;   // FIXME: verify this is kosher
					ext_keystatus[ sc ] &= ~1; // preserve bit 2 only
				}
				switch (ev.key.keysym.sym) {
					case SDLK_LSHIFT: shkeystatus &= ~(1<<16); break;
					case SDLK_RSHIFT: shkeystatus &= ~(1<<17); break;
					case SDLK_LCTRL:  shkeystatus &= ~(1<<18); break;
					case SDLK_RCTRL:  shkeystatus &= ~(1<<19); break;
					case SDLK_LALT:   shkeystatus &= ~(1<<20); break;
					case SDLK_RALT:   shkeystatus &= ~(1<<21); break;
					default:
						break;
				}
				break;
			case SDL_MOUSEMOTION:
				if (ActiveApp) {
					mousex += ev.motion.xrel;
					mousey += ev.motion.yrel;
				}
				break;

			case SDL_MOUSEBUTTONDOWN:
			case SDL_MOUSEBUTTONUP:
				{
					int i=-1;
					switch (ev.button.button) {
						case SDL_BUTTON_LEFT:   i = 0; break;
						case SDL_BUTTON_RIGHT:  i = 1; break;
						case SDL_BUTTON_MIDDLE: i = 2; break;
					}
					if (i<0) break;
					
					if (ev.type == SDL_MOUSEBUTTONDOWN) {
						ext_mbstatus[i] = 1|2;
						gbstatus |= 1<<i;
					} else {
						ext_mbstatus[i] &= 2;
						gbstatus &= ~(1<<i);
					}
				}
				break;
			case SDL_QUIT:
				kensoundclose();
				SDL_DestroyTexture(texture);
				SDL_DestroyRenderer(renderer);
				SDL_DestroyWindow(screen);
				SDL_Quit();
				uninitapp();
				quitprogram = 1;
				break;
		}
	}
}

long keyread ()
{
	long i;

	if (keybufr == keybufw) return(0);
	i = keybuf[keybufr]; keybufr = ((keybufr+1)&(KEYBUFSIZ-1));
	return(i);
}

void breath ()
{
	sdlmsgloop();
}

	//Call like this: arg2filename(argv[1],".ksm",curfilename);
	//Make sure curfilename length is lon(wParam&255)g enough!
#if !defined(_MSC_VER)
static void strlwr(char *s)
{
		for ( ; *s; s++)
				if ((*s) >= 'A' || (*s) <= 'Z')
						*s &= ~0x20;
}
#endif
void arg2filename (const char *oarg, const char *ext, char *narg)
{
	long i;

		//Chop off quotes at beginning and end of filename
	for(i=strlen(oarg);i>0;i--)
		if (oarg[i-1] == '\\') break;
	if ((!i) && (oarg[0] == '\"')) i = 1;
	strcpy(narg,&oarg[i]);
	if (narg[0] == 0) return;
	if (narg[strlen(narg)-1] == '\"') narg[strlen(narg)-1] = 0;
	strlwr(narg);
	if (!strstr(narg,ext)) strcat(narg,ext);
}

	//Precision: bits 8-9:, Rounding: bits 10-11:
	//00 = 24-bit    (0)   00 = nearest/even (0)
	//01 = reserved  (1)   01 = -inf         (4)
	//10 = 53-bit    (2)   10 = inf          (8)
	//11 = 64-bit    (3)   11 = 0            (c)
long fpuasm[2];

#if defined(_MSC_VER)
static _inline void fpuinit (long a)
{
	_asm
	{
		mov eax, a
		fninit
		fstcw fpuasm
		and byte ptr fpuasm[1], 0f0h
		or byte ptr fpuasm[1], al
		fldcw fpuasm
	}
}
#else
void fpuinit (long a)
{
	__asm__ __volatile__ (
		"fninit; fstcw fpuasm; andb $240, fpuasm(,1); orb %%al, fpuasm(,1); fldcw fpuasm;"
		: : "a" (a) : "memory","cc");
}
#endif

int main(int argc, char *argv[])
{
	Uint32 sdlinitflags;
	int i;
	
	sdlinitflags = SDL_INIT_TIMER;
	for(i=0;i<256;i++) keystatus[i] = 0;
	for(i=0;i<256;i++) ext_keystatus[i] = 0;
	sdlinitflags |= SDL_INIT_VIDEO;

	if (SDL_Init(sdlinitflags) < 0) { fputs("Failure initialising SDL.",stderr); return -1; }
	atexit(SDL_Quit);   // In case we exit() somewhere
	
	initklock();

	if (initapp(argc, argv) < 0) return -1;

	if (!initdirectdraw(xres,yres,colbits))
	{
		SDL_Quit();
		return 0;
	}

	breath();
	while (!quitprogram)
	{
		SDL_Delay(10);

		if (alwaysactive || ActiveApp) { fpuinit(0x2); doframe(); }
		else SDL_WaitEvent(NULL);
		breath();
	}

	return quitparam;
}
