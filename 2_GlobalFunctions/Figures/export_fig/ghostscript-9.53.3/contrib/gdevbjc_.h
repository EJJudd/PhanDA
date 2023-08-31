/* BJC-210/240/250/265/1000 Bubble Jet Printer driver for GhostScript
 * header file
 *
 * Copyright 2000, 2001, 2002 Gergely Sz�sz (Gergely Sza'sz)
 * mailto://szaszg@hu.inter.net http://bjc250gs.sourceforge.net
 *
 *   This program may be distributed and/or modified under the terms of
 *   the GNU General Public License as published by the Free Software
 *   Foundation (the "GPL"); either version 2 of the GPL, or (at your option)
 *   any later version.
 *
 *   When distributed under the terms of the GPL, this program is distributed
 *   in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
 *   even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 *   PURPOSE.  See the GPL for more details.
 *
 *   If this program is being distributed under the terms of the GPL, you
 *   should have received a copy of the GPL along with this program, normally
 *   in a plain ASCII text file named COPYING; if not, write to the:
 *
 *   Free Software Foundation, Inc.,
 *   51 Franklin Street, Fifth Floor
 *   Boston, MA 02110-1301
 *   USA
 *
 */

/* Copyright (C) 1989, 2000 Aladdin Enterprises.  All rights reserved.

   This program may also be distributed as part of AFPL Ghostscript, under the
   terms of the Aladdin Free Public License (the "License").

   AFPL Ghostscript is distributed with NO WARRANTY OF ANY KIND.  No
   author or distributor accepts any responsibility for the consequences of
   using it, or for whether it serves any particular purpose or works at
   all, unless he or she says so in writing.  Refer to the License for full
   details.

   Every copy of AFPL Ghostscript must include a copy of the License,
   normally in a plain ASCII text file named PUBLIC.  The License grants you
   the right to copy, modify and redistribute AFPL Ghostscript, but only
   under certain conditions described in the License.  Among other things,
   the License requires that the copyright notice and this notice be
   preserved on all copies.
*/

/* BJC printer drivers definitions and utility interfaces */

#ifndef gdevbjcUNDERSCORE_INCLUDED
#  define gdevbjcUNDERSCORE_INCLUDED

/* Define the default X and Y resolution. */
#define X_DPI 360
#define Y_DPI 360

#define CMYK_K 'K'
#define CMYK_C 'C'
#define CMYK_M 'M'
#define CMYK_Y 'Y'

#define INK_K 8
#define INK_C 1
#define INK_M 2
#define INK_Y 4

typedef struct {
    gs_param_string p_string;
    int             p_value;
} stringParamDescription;

/* Definitions of device parameters */
struct skip_s {
    bool skipC;
    bool skipM;
    bool skipY;
    bool skipK;
};

typedef struct skip_s skip_t;

typedef struct {
    unsigned char l;
    unsigned char c;
} media_t;

struct gx_device_bjc_printer_s {
    gx_device_common;
    gx_prn_device_common;
    int printerType;         /* BJC-250, BJC-250ex, BJC1000, ... */
    int feeder;              /* Needed paper feeder */
    int mediaType;           /* paper type */
    int quality;
    int ink;
    int rnd;
    bool inverse;
    bool compress;
    bool smooth;
    bool limit;               /* paper color correction limit check */
    bool compose;             /* color K decomposition */
    float gamma;
    float redGamma;
    float greenGamma;
    float blueGamma;
    struct {
        int red;
        int green;
        int blue;
    } paperColor;                      /* paper color for color correction */
    /* Global variables from gdevbjca moved here */
    int bjc_j; /* =0 */
    int bjc_k; /* =31 */
    int bjc_treshold[1024];
    bool FloydSteinbergDirectionForward;
    int *FloydSteinbergErrorsC;
    int *FloydSteinbergErrorsM;
    int *FloydSteinbergErrorsY;
    int *FloydSteinbergErrorsK;
    int *FloydSteinbergErrorsG;
    int FloydSteinbergC;
    int FloydSteinbergM;
    int FloydSteinbergY;
    int FloydSteinbergK;
    int FloydSteinbergG;
    int bjc_gamma_tableC[256];
    int bjc_gamma_tableM[256];
    int bjc_gamma_tableY[256];
};

typedef struct gx_device_bjc_printer_s gx_device_bjc_printer;

#define bjc_device_margins_body(dtype, procs, dname, w10, h10, xdpi, ydpi, lo, to, lm, bm, rm, tm, ncomp, depth, mg, mc, dg, dc, print_page)\
        std_device_full_body_type(dtype, &procs, dname, &st_device_printer,\
          (int)((long)(w10) * (xdpi) / 10),\
          (int)((long)(h10) * (ydpi) / 10),\
          xdpi, ydpi,\
          ncomp, depth, mg, mc, dg, dc,\
          -(lo) * (xdpi), -(to) * (ydpi),\
          (lm) * 72.0f, (bm) * 72.0f,\
          (rm) * 72.0f, (tm) * 72.0f\
        ),\
        prn_device_body_rest_(print_page)

#define bjc_device_margins(procs, dname, w10, h10, xdpi, ydpi, lo, to, lm, bm, rm, tm, ncomp, depth, mg, mc, dg, dc, print_page, def_ink)\
{ bjc_device_margins_body(gx_device_bjc_printer, procs, dname,\
    w10, h10, xdpi, ydpi, lo, to, lm, bm, rm, tm, \
    ncomp, depth, mg, mc, dg, dc, print_page),\
    (int) 0,                /* def, printer: BJC250 */ \
    (int) 0x10,             /* Auto Feeder */          \
    (int) 0,                /*  Plain paper */         \
    (int) 0,                /*  Normal quality */      \
    (int) def_ink,          /* Black monochrome ink */ \
    (int) 15,               /* Randomness is 15% */    \
    (bool) false,           /* Normal printing */      \
    (bool) true,            /* Use compression */      \
    (bool) false,           /* no smoothing */         \
    (bool) false,           /* no limit check */       \
    (bool) true,            /* yes we run decomp */    \
    (float) 1.0,            /* initiaal gamma */       \
    (float) 1.0,            /* initiaal gamma */       \
    (float) 1.0,            /* initiaal gamma */       \
    (float) 1.0,            /* initiaal gamma */       \
    {  (int) 255,           /* White  paper  */        \
       (int) 255,           /*               */        \
       (int) 255 }          /*               */        \
}

#define bjc_device(procs, dname, w10, h10, xdpi, ydpi, lm, bm, rm, tm, ncomp, depth, mg, mc, dg, dc, print_page, def_ink)\
  bjc_device_margins(procs, dname, w10, h10, xdpi, ydpi,\
    lm, tm, lm, bm, rm, tm, ncomp, depth, mg, mc, dg, dc, print_page, def_ink)

#define bjc_cmyk_param_procs(v_prn_open, v_prn_output_page, v_prn_close, \
                 p_map_color_rgb, p_map_cmyk_color, \
                 v_prn_get_params, v_prn_put_params)\
   {v_prn_open, NULL, NULL, v_prn_output_page, v_prn_close,\
    NULL, p_map_color_rgb, NULL, NULL, NULL, NULL, NULL, NULL,\
    v_prn_get_params, v_prn_put_params,\
    p_map_cmyk_color, NULL, NULL, NULL, gx_page_device_get_page_device}

/* There are the definitions of commands for the Canon BJC printers. */

#define BJC_BJL_CLEANING     0x00  /*   "@Cleaning=1\n" */
#define BJC_BJL_POWER_OFF    0x01  /*   "@PowerOff\n" */
#define BJC_BJL_ROCLEANING   0x02  /*   "@RollerCleaning\n" */
#define BJC_BJL_NOZZLECHK    0x03  /*   "@TestPrint=NozzleCheck\n" */
#define BJC_BJL_TESTA        0x04  /*   "@TestPrint=A\012" */
#define BJC_BJL_TESTB        0x05  /*   "@TestPrint=B\012" */
#define BJC_BJL_TESTC        0x06  /*   "@TestPrint=A\012" */
#define BJC_BJL_DEMO         0x07  /*   "@TestPrint=DemoPrint\012" */

#define BJC_BJL_ON_ENABLE    0x20  /*   "AutoPowerOn=Enable\n" */
#define BJC_BJL_ON_DISABLE   0x21  /*   "AutoPowerOn=Disable\n" */
#define BJC_BJL_OFF_1MIN     0x22  /*   "AutoPowerOff=1\n" */
#define BJC_BJL_OFF_10MIN    0x23  /*   "AutoPowerOff=10\n" */
#define BJC_BJL_OFF_30MIN    0x24  /*   "AutoPowerOff=30\n" */
#define BJC_BJL_OFF_60MIN    0x25  /*   "AutoPowerOff=60\n" */
#define BJC_BJL_OFF_DISABLE  0x26  /*   "AutoPowerOff=Disable\n" */

typedef struct {
    const char *string;
    int         numeric;
    int         length; } BJL_command;

/* Put a BJL command to stream */

void bjc_put_bjl_command(gp_file * file, int bjl_command);

void bjc_put_LF(gp_file *file);
void bjc_put_FF(gp_file *file);
void bjc_put_CR(gp_file *file);
void bjc_put_initialize(gp_file *file);
void bjc_put_set_initial(gp_file *file);
void bjc_put_set_compression(gp_file *file, char compression);
void bjc_put_print_method_short(gp_file *file, char color);
void bjc_put_print_method(gp_file *file, char color, char media, char quality, char density);
void bjc_put_raster_resolution(gp_file *file, int x_resolution, int y_resolution);
void bjc_put_raster_skip(gp_file *file, int skip);
void bjc_put_page_margins(gp_file *file, int length, int lm, int rm, int top);
void bjc_put_media_supply(gp_file *file, char supply, char type);
void bjc_put_cmyk_image(gp_file *file, char component, const byte *data, int count);
void bjc_put_move_lines(gp_file *file, int lines);
void bjc_put_move_lines_unit(gp_file *file, int unit);
void bjc_put_extended_margins(gp_file *file, int length, int lm, int rm, int top);
void bjc_put_image_format(gp_file *file, char depth, char format, char ink);
void bjc_put_page_id(gp_file *file, int id);
void bjc_put_continue_image(gp_file *file, const byte *data, int count);
void bjc_put_indexed_image(gp_file *file, int dot_rows, int dot_cols, int layers);

bool bjc_invert_bytes(byte *row, uint raster, bool inverse, byte lastmask);
bool bjc_invert_cmyk_bytes(byte *rowC,byte *rowM, byte *rowY, byte *rowK, uint raster,
                           bool inverse, byte lastmask, skip_t *skip);
uint bjc_compress(const byte *row, uint raster, byte *compressed);

int  FloydSteinbergInitG(gx_device_printer * pdev);
void FloydSteinbergDitheringG(gx_device_bjc_printer *dev,
                              byte *row, byte *dithered, uint width, uint raster, bool limit_extr);
void FloydSteinbergCloseG(gx_device_printer *pdev);

int  FloydSteinbergForwardCut(int error, int *Errors, int i, byte *dithered, byte bitmask);
int  FloydSteinbergBckwardCut(int error, int *Errors, int i, byte *dithered, byte bitmask);
int  FloydSteinbergInitC(gx_device_printer * pdev);
void FloydSteinbergDitheringC(gx_device_bjc_printer *dev,
                              byte *row, byte *dithered, uint width, uint raster,
                              bool limit_extr, bool composeK);
void FloydSteinbergCloseC(gx_device_printer *pdev);

void bjc_build_gamma_table(gx_device_bjc_printer *dev, float gamma, char color);
void bjc_rgb_to_cmy (byte r, byte g, byte b, int *c, int *m, int *y);
void bjc_rgb_to_gray(byte r, byte g, byte b, int *k);
uint bjc_rand(gx_device_bjc_printer *dev);
void bjc_init_tresh(gx_device_bjc_printer *dev, int rnd);
#endif				/* gdevbjcUNDERSCORE_INCLUDED */
