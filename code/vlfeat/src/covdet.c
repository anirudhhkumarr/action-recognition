/** @internal
 ** @file     covdet.c
 ** @author   Andrea Vedaldi
 ** @author   Karel Lenc
 ** @brief    Driver for covariant feature detectors
 **/

/*
Copyright (C) 2007-12 Andrea Vedaldi and Brian Fulkerson.
All rights reserved.

This file is part of the VLFeat library and is made available under
the terms of the BSD license (see the COPYING file).
*/

#define VL_COVDET_DRIVER_VERSION 0.1

#include "generic-driver.h"
#include "scalespace-export.h"

#include <vl/generic.h>
#include <vl/stringop.h>
#include <vl/pgm.h>
#include <vl/covdet.h>
#include <vl/scalespace.h>
#include <vl/getopt_long.h>

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

/* ----------------------------------------------------------------- */
/* help message */
char const help_message [] =
  "Usage: %s [options] files ...\n"
  "\n"
  "Options include:\n"
  "Driver options:\n"
  " --verbose -v      Be verbose\n"
  " --help -h         Print this help message\n"
  " --output -o       Specify output file\n"
  " --frames          Specify frames file\n"
  " --descriptors     Specify descriptors file\n"
  " --meta            Specify meta file\n"
  " --gss             Output Gaussian scale space files\n"
  " --read-frames     Specify a file from which to read frames\n"
  " --ifrms-type      Specify the type of frames in the input file\n"
  "\n"
  "Covariant frames detector options:\n"
  " --method -m       Method of frames localisation [log | hessian]\n"
  " --affine-adapt -a Calculate affine covariant frames \n"
  " --orientation -r  Calculate rotation covariant frames \n"
  " --octaves -O      Number of octaves, -1 for max. (-1)\n"
  " --levels -S       Number of levels per octave (3)\n"
  " --first-octave    Index of the first octave (-1)\n"
  " --edge-thresh     Specify the edge threshold (0 for DoG, 28.5 for Hessian)\n"
  " --peak-thresh     Specift the peak threshold (10)\n"
  " --magnif          Specify the magnification factor for SIFT descriptor\n"
  "                   (3 for discs, 3*sqrt(3) for ellipses)"
  " --norm-thresh     Specify the normalisation threshold for SIFT descriptor (-inf)\n"
  " --calc-inv-sm     Export shape matrix as its inverse (Oxford format)"
  " --no-descriptors  Do not calculate descriptors\n"
  "\n"
  "Disc frames (without affine shape) options:"
  " --win-size        Gaussian window size for the sift descriptor (2)\n"
  "\n"
  "Ellipse frames options:"
  " --aff-win-size    Size of the window for affine shape est (19)\n"
  " --aff-max-iter    Max. num. of iteration for aff. shape est (16)\n"
  " --aff-conv-thr    Convergence threshold for aff. shape est. (0.05)\n"
  " --sift-patch-size Size of the patch used for orient. comp. and SIFT desc (41)\n"
  "\n" ;

/* ----------------------------------------------------------------- */
/* long options codes */
enum {
  opt_meta = 1000,
  opt_frames,
  opt_descriptors,
  opt_gss,
  opt_first_octave,
  opt_edge_thresh,
  opt_peak_thresh,
  opt_magnif,
  opt_aff_win_size,
  opt_aff_max_iter,
  opt_aff_conv_thr,
  opt_patch_size,
  opt_norm_thresh,
  opt_window_size,
  opt_calc_inv_sm,
  opt_no_descr,
  opt_read_frames,
  opt_ifrms_type
} ;

/* short options */
char const opts [] = "vhO:S:o:m:ar" ;

/* long options */
struct option const longopts [] = {
  { "verbose",         no_argument,            0,          'v'              },
  { "help",            no_argument,            0,          'h'              },
  { "octaves",         required_argument,      0,          'O'              },
  { "levels",          required_argument,      0,          'S'              },
  { "output",          optional_argument,      0,          'o'              },
  { "method",          required_argument,      0,          'm'              },
  { "affine-adapt",    no_argument,            0,          'a'              },
  { "orientation",     no_argument,            0,          'r'              },
  { "meta",            optional_argument,      0,          opt_meta         },
  { "frames",          optional_argument,      0,          opt_frames       },
  { "descriptors",     optional_argument,      0,          opt_descriptors  },
  { "gss",             optional_argument,      0,          opt_gss          },
  { "first-octave",    required_argument,      0,          opt_first_octave },
  { "edge-thresh",     required_argument,      0,          opt_edge_thresh  },
  { "peak-thresh",     required_argument,      0,          opt_peak_thresh  },
  { "magnif",          required_argument,      0,          opt_magnif       },
  { "norm-thresh",     required_argument,      0,          opt_norm_thresh  },
  { "win-size",        required_argument,      0,          opt_window_size  },
  { "aff-win-size",    required_argument,      0,          opt_aff_win_size },
  { "aff-max-iter",    required_argument,      0,          opt_aff_max_iter },
  { "aff-conv-thr",    required_argument,      0,          opt_aff_conv_thr },
  { "sift-patch-size", required_argument,      0,          opt_patch_size   },
  { "calc-inv-sm",     no_argument,            0,          opt_calc_inv_sm  },
  { "no-descriptors",  no_argument,            0,          opt_no_descr     },
  { "read-frames",     required_argument,      0,          opt_read_frames  },
  { "ifrms-type",      required_argument,      0,          opt_ifrms_type   },
  { 0,                 0,                      0,          0                }
} ;

/* ----------------------------------------------------------------- */
/** @internal
 ** @brief Keypoint ordering
 ** Sort callback for ordering keypoints in a list.
 **/
int
korder (void const* a, void const* b) {
  double x = ((double*) a) [2] - ((double*) b) [2] ;
  if (x < 0) return -1 ;
  if (x > 0) return +1 ;
  return 0 ;
}

/** @internal
 ** @brief Save a frame into file
 **
 ** @param self Meta file structure
 ** @param frame Pointer to the frame
 ** @param frm_type Frame type
 **/
VL_INLINE void const *
vl_file_meta_put_frame(VlFileMeta * self, void const* frame,
                       VlFrameType frameType)
{
  switch (frameType) {
  case VL_FRAMETYPE_DISC: {
    VlFrameDisc const *frm = (VlFrameDisc const *)frame;
    vl_file_meta_put_double (self, frm -> x) ;
    vl_file_meta_put_double (self, frm -> y) ;
    vl_file_meta_put_double (self, frm -> sigma) ;
    frame = (void const *)(frm + 1);
  } break;
  case VL_FRAMETYPE_ORIENTED_DISC:{
    VlFrameOrientedDisc const *frm = (VlFrameOrientedDisc const *)frame;
    vl_file_meta_put_double (self, frm -> x) ;
    vl_file_meta_put_double (self, frm -> y)  ;
    vl_file_meta_put_double (self, frm -> sigma) ;
    vl_file_meta_put_double (self, frm->angle);
    frame = (void const *)(frm + 1);
  } break;
  case VL_FRAMETYPE_ELLIPSE:{
    VlFrameEllipse const *frm = (VlFrameEllipse const *)frame;
    vl_file_meta_put_double (self, frm -> x) ;
    vl_file_meta_put_double (self, frm -> y) ;
    vl_file_meta_put_double (self, frm -> e11) ;
    vl_file_meta_put_double (self, frm -> e12) ;
    vl_file_meta_put_double (self, frm -> e22) ;
    frame = (void const *)(frm + 1);
  } break;
  case VL_FRAMETYPE_ORIENTED_ELLIPSE:{
    VlFrameOrientedEllipse const *frm =
        (VlFrameOrientedEllipse const *)frame;
    vl_file_meta_put_double (self, frm -> x) ;
    vl_file_meta_put_double (self, frm -> y) ;
    vl_file_meta_put_double (self, frm -> a11) ;
    vl_file_meta_put_double (self, frm -> a12) ;
    vl_file_meta_put_double (self, frm -> a21) ;
    vl_file_meta_put_double (self, frm -> a22) ;
    frame = (void const *)(frm + 1);
  } break;
  default:
    VL_ASSERT(0,"Invalid frame type.");
    break;
  }

  return frame;
}


int vl_file_meta_read_double_array(VlFileMeta* file, double* array, int length)
{
  int i, err;

  for (i = 0; i < length; ++i){
    err = vl_file_meta_get_double (file, array++);
    if (err) {
      return err;
    }
  }
  return err;
}

/* ---------------------------------------------------------------- */
/** @brief SIFT driver entry point
 **/
int
main(int argc, char **argv)
{
  /* algorithm parameters */
  double   edge_thresh  = -1 ;
  double   peak_thresh  = -1 ;
  int      O = -1, S = 3, o_min = -1 ;
  double   magnif        = -1 ;
  double   norm_thresh    = -1;
  double   window_size      = -1;
  int      aff_win_size  = 19 ;
  int      aff_max_iter  = -1 ;
  double   aff_conv_thr  = -1 ;
  int      patch_size    = 41 ;

  VlFrameType frm_type = VL_FRAMETYPE_DISC;
  VlFrameType inputFrameType = VL_FRAMETYPE_DISC;
  VlCovDetMethod method = VL_COVDET_METHOD_DOG;

  vl_bool err = VL_ERR_OK ;
  char err_msg [1024] ;
  int n ;
  int exit_code = 0 ;
  int verbose = 0 ;
  vl_bool force_output = VL_FALSE ;
  vl_bool calc_descs = VL_TRUE ;
  vl_bool calc_inv_sm = VL_FALSE ;
  vl_bool calc_orient = VL_FALSE;
  vl_bool calc_affine = VL_FALSE;
  VlEnumerator * pair ;

  VlFileMeta out  = {1, "%.sift",  VL_PROT_ASCII, "", 0} ;
  VlFileMeta frm  = {0, "%.frame", VL_PROT_ASCII, "", 0} ;
  VlFileMeta dsc  = {0, "%.descr", VL_PROT_ASCII, "", 0} ;
  VlFileMeta met  = {0, "%.meta",  VL_PROT_ASCII, "", 0} ;
  VlFileMeta gsf  = {0, "%.pgm",   VL_PROT_ASCII, "", 0} ;
  VlFileMeta ifr  = {0, "%.frame", VL_PROT_ASCII, "", 0} ;

#define ERRF(msg, arg) {                                        \
    err = VL_ERR_BAD_ARG ;                                      \
    snprintf(err_msg, sizeof(err_msg), msg, arg) ;              \
    break ;                                                     \
  }

#define ERR(msg) {                                              \
    err = VL_ERR_BAD_ARG ;                                      \
    snprintf(err_msg, sizeof(err_msg), msg) ;                   \
    break ;                                                     \
}

  /* -----------------------------------------------------------------
   *                                                     Parse options
   * -------------------------------------------------------------- */

  while (!err) {
    int ch = getopt_long(argc, argv, opts, longopts, 0) ;

    /* If there are no files passed as input, print the help and settings */
    if (ch == -1 && argc - optind == 0)
      ch = 'h';

    /* end of option list? */
    if (ch == -1) break;

    switch (ch) {

    case '?' :
      /* unkown option ............................................ */
      ERRF("Invalid option '%s'.", argv [optind - 1]) ;
      break ;

    case ':' :
      /* missing argument ......................................... */
      ERRF("Missing mandatory argument for option '%s'.",
          argv [optind - 1]) ;
      break ;

    case 'h' :
      /* --help ................................................... */
      printf (help_message, argv [0]) ;
      printf ("SIFT         filespec: `%s'\n", out.pattern) ;
      printf ("Frames       filespec: `%s'\n", frm.pattern) ;
      printf ("Descriptors  filespec: `%s'\n", dsc.pattern) ;
      printf ("Meta         filespec: `%s'\n", met.pattern) ;
      printf ("GSS          filespec: '%s'\n", gsf.pattern) ;
      printf ("Read frames  filespec: '%s'\n", ifr.pattern) ;
      printf ("Version: driver %s; libvl %s\n",
              VL_XSTRINGIFY(VL_SIFT_DRIVER_VERSION),
              vl_get_version_string()) ;
      exit (0) ;
      break ;

    case 'v' :
      /* --verbose ................................................ */
      ++ verbose ;
      break ;

    case 'o' :
      /* --output  ................................................ */
      err = vl_file_meta_parse (&out, optarg) ;
      if (err)
        ERRF("The arguments of '%s' is invalid.", argv [optind - 1]) ;
      force_output = 1 ;
      break ;

    case opt_frames :
      /* --frames  ................................................ */
      err = vl_file_meta_parse (&frm, optarg) ;
      if (err)
        ERRF("The arguments of '%s' is invalid.", argv [optind - 1]) ;
      break ;

    case opt_descriptors :
      /* --descriptor ............................................. */
      err = vl_file_meta_parse (&dsc, optarg) ;
      if (err)
        ERRF("The arguments of '%s' is invalid.", argv [optind - 1]) ;
      break;

    case opt_meta :
      /* --meta ................................................... */
      err = vl_file_meta_parse (&met, optarg) ;
      if (err)
        ERRF("The arguments of '%s' is invalid.", argv [optind - 1]) ;

      if (met.protocol != VL_PROT_ASCII)
        ERR("meta file supports only ASCII protocol") ;
      break ;

    case opt_read_frames :
      /* --read_frames ............................................ */
      err = vl_file_meta_parse (&ifr, optarg) ;
      if (err)
        ERRF("The arguments of '%s' is invalid.", argv [optind - 1]) ;
      break ;

    case opt_ifrms_type:
      /* --method ................................................ */
      pair = vl_enumeration_get_casei(vlFrameTypes, optarg);
      if (pair == NULL) {
        ERRF("Invalid method of frame localisation '%s'.", optarg) ;
        }
      inputFrameType = pair->value ;
      break;
    case opt_gss :
      /* --gss .................................................... */
      err = vl_file_meta_parse (&gsf, optarg) ;
      if (err)
        ERRF("The arguments of '%s' is invalid.", argv [optind - 1]) ;
      break ;


    case 'O' :
      /* --octaves ............................................... */
      n = sscanf (optarg, "%d", &O) ;
      if (n == 0 || O < 0)
        ERRF("The argument of '%s' must be a non-negative integer.",
            argv [optind - 1]) ;
      break ;

    case 'S' :
      /* --levels ............................................... */
      n = sscanf (optarg, "%d", &S) ;
      if (n == 0 || S < 0)
        ERRF("The argument of '%s' must be a non-negative integer.",
            argv [optind - 1]) ;
      break ;

    case 'm':
      /* --method ................................................ */
      pair = vl_enumeration_get_casei(vlCovdetMethods, optarg);
      if (pair == NULL) {
        ERRF("Invalid method of frame localisation '%s'.", optarg) ;
        }
      method = pair->value ;

    case 'a' :
      /* --affine-adapt............................................ */
      calc_affine = !calc_affine;
      break ;

    case 'r' :
      /* --orientation............................................ */
      calc_orient = !calc_orient;
      break ;

    case opt_first_octave :
      /* --first-octave ......................................... */
      n = sscanf (optarg, "%d", &o_min) ;
      if (n == 0)
        ERRF("The argument of '%s' must be an integer.",
            argv [optind - 1]) ;
      break ;

    case opt_edge_thresh :
      /* --edge-thresh ........................................... */
      n = sscanf (optarg, "%lf", &edge_thresh) ;
      if (n == 0 || edge_thresh < 1)
        ERRF("The argument of '%s' must be not smaller than 1.",
            argv [optind - 1]) ;
      break ;

    case opt_peak_thresh :
      /* --edge-thresh ........................................... */
      n = sscanf (optarg, "%lf", &peak_thresh) ;
      if (n == 0 || peak_thresh < 0)
        ERRF("The argument of '%s' must be a non-negative float.",
            argv [optind - 1]) ;
      break ;

    case opt_magnif :
      /* --magnif  .............................................. */
      n = sscanf (optarg, "%lf", &magnif) ;
      if (n == 0 || magnif < 1)
        ERRF("The argument of '%s' must be a non-negative float.",
            argv [optind - 1]) ;
      break ;

    case opt_window_size :
      /* --win-size  ............................................ */
      n = sscanf (optarg, "%lf", &window_size) ;
      if (n == 0 || window_size < 1)
        ERRF("The argument of '%s' must be a non-negative float.",
            argv [optind - 1]) ;
      break ;


    case opt_norm_thresh :
      /* --win-size  ............................................ */
      n = sscanf (optarg, "%lf", &norm_thresh) ;
      if (n == 0 || norm_thresh < 1)
        ERRF("The argument of '%s' must be a non-negative float.",
            argv [optind - 1]) ;
      break ;

    case opt_aff_win_size :
      n = sscanf (optarg, "%d", &aff_win_size) ;
      if (n == 0 || aff_win_size < 1) {
        ERRF("The argument of '%s' must be a positive integer.",
             argv [optind - 1]) ;
      }
      break ;

    case opt_aff_max_iter :
      n = sscanf (optarg, "%d", &aff_max_iter) ;
      if (n == 0 || aff_max_iter < 1) {
        ERRF("The argument of '%s' must be a positive integer.",
             argv [optind - 1]) ;
      }
      break ;

    case opt_aff_conv_thr :
      n = sscanf (optarg, "%lf", &aff_conv_thr) ;
      if (n == 0 || aff_conv_thr < 0) {
        ERRF("The argument of '%s' must be a non-negative real.",
             argv [optind - 1]) ;
      }
      break ;

    case opt_patch_size :
      n = sscanf (optarg, "%d", &patch_size) ;
      if (n == 0 || patch_size < 1) {
       ERRF("The argument of '%s' must be a positive integer.",
            argv [optind - 1]) ;
      }
      break ;

    case opt_calc_inv_sm :
      /* --calc-inv-sm .......................................... */
      calc_inv_sm = !calc_inv_sm;
      break ;

    case opt_no_descr :
      /* --no-descriptors ....................................... */
      calc_descs = !calc_descs ;
      break ;

    case 0 :
    default :
      /* should not get here ...................................... */
      fprintf(stderr, "Invalid argument.\n") ;
      abort() ;
    }
  }

  /* check for parsing errors */
  if (err) {
    fprintf(stderr, "%s: error: %s (%d)\n",
            argv [0],
            err_msg, err) ;
    exit (1) ;
  }

  /* parse other arguments (filenames) */
  argc -= optind ;
  argv += optind ;

  /*
     if --output is not specified, specifying --frames or --descriptors
     prevent the aggregate outout file to be produced.
  */
  if (! force_output && (frm.active || dsc.active)) {
    out.active = 0 ;
  }

  frm_type = vl_frame_get_type(calc_affine, calc_orient);

  if (verbose > 1) {
#define PRNFO(name,fm)                                                  \
    printf("covdet: " name) ;                                           \
    printf("%3s ",  (fm).active ? "yes" : "no") ;                       \
    printf("%-6s ", vl_string_protocol_name ((fm).protocol)) ;          \
    printf("%-10s\n", (fm).pattern) ;

    PRNFO("write aggregate . ", out) ;
    PRNFO("write frames .... ", frm) ;
    PRNFO("write descriptors ", dsc) ;
    PRNFO("write meta ...... ", met) ;
    PRNFO("write GSS ....... ", gsf) ;
    PRNFO("read  frames .... ", ifr) ;

  }

  /* ------------------------------------------------------------------
   *                                         Process one image per time
   * --------------------------------------------------------------- */

  while (argc--) {

    char             basename [1024] ;
    char const      *name = *argv++ ;

    FILE            *in    = 0 ;
    vl_uint8        *data  = 0 ;
    float           *fdata = 0 ;
    VlPgmImage       pim ;

    VlCovDet        *det = 0;
    vl_size          q ;

    float           *inputFrames = 0 ;
    int              numInputFrames = 0, inputFramesStorageSize = 0 ;

    /* ...............................................................
     *                                                 Determine files
     * ............................................................ */

    /* get basenmae from filename */
    q = vl_string_basename (basename, sizeof(basename), name, 1) ;

    err = (q >= sizeof(basename)) ;

    if (err) {
      snprintf(err_msg, sizeof(err_msg),
               "Basename of '%s' is too long", name);
      err = VL_ERR_OVERFLOW ;
      goto done ;
    }

    if (verbose) {
      printf ("covdet: <== '%s'\n", name) ;
    }

    if (verbose > 1) {
      printf ("covdet: basename is '%s'\n", basename) ;
    }

    /* open input file */
    in = fopen (name, "rb") ;
    if (!in) {
      err = VL_ERR_IO ;
      snprintf(err_msg, sizeof(err_msg),
               "Could not open '%s' for reading.", name) ;
      goto done ;
    }

    /* ...............................................................
     *                                                       Read data
     * ............................................................ */

    /* read PGM header */
    err = vl_pgm_extract_head (in, &pim) ;

    if (err) {
      switch (vl_get_last_error()) {
      case  VL_ERR_PGM_IO :
        snprintf(err_msg, sizeof(err_msg),
                 "Cannot read from '%s'.", name) ;
        err = VL_ERR_IO ;
        break ;

      case VL_ERR_PGM_INV_HEAD :
        snprintf(err_msg, sizeof(err_msg),
                 "'%s' contains a malformed PGM header.", name) ;
        err = VL_ERR_IO ;
        goto done ;
      }
    }

    if (verbose)
      printf ("covdet: image is %d by %d pixels\n",
              (signed)pim. width,
              (signed)pim. height) ;

    /* allocate buffer */
    data  = malloc(vl_pgm_get_npixels (&pim) *
                   vl_pgm_get_bpp       (&pim) * sizeof (vl_uint8)   ) ;
    fdata = malloc(vl_pgm_get_npixels (&pim) *
                   vl_pgm_get_bpp       (&pim) * sizeof (float)) ;

    if (!data || !fdata) {
      err = VL_ERR_ALLOC ;
      snprintf(err_msg, sizeof(err_msg),
               "Could not allocate enough memory.") ;
      goto done ;
    }

    /* read PGM body */
    err  = vl_pgm_extract_data (in, &pim, data) ;

    if (err) {
      snprintf(err_msg, sizeof(err_msg), "PGM body malformed.") ;
      err = VL_ERR_IO ;
      goto done ;
    }

    /* convert data type */
    for (q = 0 ; q < (unsigned) (pim.width * pim.height) ; ++q) {
      fdata [q] = data [q] ;
    }

    /* ...............................................................
     *                                     Optionally source frames
     * ............................................................ */

#define WERR(name,op)                                           \
    if (err == VL_ERR_OVERFLOW) {                               \
      snprintf(err_msg, sizeof(err_msg),                        \
               "Output file name too long.") ;                  \
      goto done ;                                               \
    } else if (err) {                                           \
      snprintf(err_msg, sizeof(err_msg),                        \
               "Could not open '%s' for " #op, name) ;          \
      goto done ;                                               \
    }

#define QERR(err,name)\
  if (err == VL_ERR_EOF) {                                       \
    /* End of file, no more points */                           \
    hasNewFrames = VL_FALSE;                                    \
    break;                                                      \
  } else if(err) {                                              \
    snprintf(err_msg, sizeof(err_msg),                          \
             "Corrupted file '%s'.", name) ;                    \
    goto done ;                                                 \
  } else { numInputFrames++; }                                  \

    if (ifr.active) {
      vl_size inputFrameSize = vl_frame_size(inputFrameType);
      vl_bool hasNewFrames = VL_TRUE;
      double buffer[6];
      inputFramesStorageSize = 0;
      numInputFrames = 0;

      /* open file */
      err = vl_file_meta_open (&ifr, basename, "rb") ;
      WERR(ifr.name, reading) ;

      while (hasNewFrames) {

        if (inputFramesStorageSize <= numInputFrames) {
          inputFramesStorageSize += 1000;
          inputFrames = vl_realloc(inputFrames,
                                   inputFramesStorageSize * inputFrameSize);
        }

        switch (inputFrameType) {
        case VL_FRAMETYPE_DISC: {
          VlFrameDisc *frame;
          frame = (VlFrameDisc*)inputFrames + numInputFrames;
          err = vl_file_meta_read_double_array(&ifr, buffer, 3);
          QERR(err,ifr.name);
          frame->x = buffer[0];
          frame->y = buffer[1];
          frame->sigma = buffer[2];
        } break;
        case VL_FRAMETYPE_ORIENTED_DISC: {
          VlFrameOrientedDisc *frm;
          frm = (VlFrameOrientedDisc*)inputFrames + numInputFrames;
          err = vl_file_meta_read_double_array(&ifr, buffer, 4);
          QERR(err,ifr.name);
          frm->x = buffer[0];
          frm->y = buffer[1];
          frm->sigma = buffer[2];
          frm->angle = buffer[3];
        } break;
        case VL_FRAMETYPE_ELLIPSE: {
          VlFrameEllipse *frm;
          frm = (VlFrameEllipse*)inputFrames + numInputFrames;
          err = vl_file_meta_read_double_array(&ifr, buffer, 5);
          QERR(err,ifr.name);
          frm->x = buffer[0];
          frm->y = buffer[1];
          frm->e11 = buffer[2];
          frm->e12 = buffer[3];
          frm->e22 = buffer[4];
        } break;
        case VL_FRAMETYPE_ORIENTED_ELLIPSE: {
          VlFrameOrientedEllipse *frm;
          frm = (VlFrameOrientedEllipse*)inputFrames + numInputFrames;
          err = vl_file_meta_read_double_array(&ifr, buffer, 6);
          QERR(err,ifr.name);
          frm->x = buffer[0];
          frm->y = buffer[1];
          frm->a11 = buffer[2];
          frm->a12 = buffer[3];
          frm->a21 = buffer[4];
          frm->a22 = buffer[5];
        } break;
        default:
          break;
        }
      }

      if (verbose) {
        printf ("covdet: read %d frames from '%s'\n", numInputFrames, ifr.name) ;
      }

      /* close file */
      vl_file_meta_close (&ifr) ;
    }
    /* ...............................................................
     *                                               Open output files
     * ............................................................ */

    err = vl_file_meta_open (&out, basename, "wb") ; WERR(out.name, writing) ;
    err = vl_file_meta_open (&dsc, basename, "wb") ; WERR(dsc.name, writing) ;
    err = vl_file_meta_open (&frm, basename, "wb") ; WERR(frm.name, writing) ;
    err = vl_file_meta_open (&met, basename, "wb") ; WERR(met.name, writing) ;

    if (verbose > 1) {
      if (out.active) printf("covdet: writing all ....... to . '%s'\n", out.name);
      if (frm.active) printf("covdet: writing frames .... to . '%s'\n", frm.name);
      if (dsc.active) printf("covdet: writing descriptors to . '%s'\n", dsc.name);
      if (met.active) printf("covdet: writign meta ...... to . '%s'\n", met.name);
    }

    /* ...............................................................
     *                                                   Make detector
     * ............................................................ */

    {
      VlAffineShapeEstimator         *aff_det = 0 ;
      VlAffinePatchNormalizer        *aff_norm = 0 ;
      VlScaleSpace        *gss;
      void const          *frames = 0 ;
      float const         *descrs = 0 ;
      int                  nframes = numInputFrames, i, err ;
      vl_size              desc_sz = 0;

      /* create a filter to process the image */
      if (calc_affine){
        det = vl_covdet_new_ellipse_detector (pim.width, pim.height,
                                                method, O, o_min, S,
                                                aff_win_size, patch_size,
                                                calc_orient, calc_descs) ;
      } else {
        det = vl_covdet_new_disc_detector (pim.width, pim.height,
                                             method, O, o_min, S,
                                             calc_orient, calc_descs) ;
      }
      gss = vl_covdet_get_gss(det);
      if (calc_affine) {
        aff_det = vl_covdet_get_affdet(det);
        aff_norm = vl_covdet_get_affnorm(det);
      }
      if (calc_descs){
        desc_sz = vl_covdet_get_descriptor_size(det);
      }

      if (peak_thresh >= 0)  vl_covdet_set_peak_thresh (det, peak_thresh) ;
      if (edge_thresh >= 0)  vl_covdet_set_edge_thresh (det, edge_thresh) ;
      if (calc_descs) {
        if (norm_thresh >= 0)  vl_covdet_set_norm_thresh (det, norm_thresh);
        if (magnif      >= 0)  vl_covdet_set_magnif      (det, magnif) ;
        if (!calc_affine) {
          /* This parameter used only for ellipse detection  */
          if (window_size >= 0)  vl_covdet_set_window_size (det, window_size);
        }
      }
      if (calc_affine) {
        if (aff_max_iter >= 0) vl_affineshapeestimator_set_max_iter    (aff_det, aff_max_iter);
        if (aff_conv_thr >= 0) vl_affineshapeestimator_set_conv_thresh (aff_det, aff_conv_thr);
        if (calc_inv_sm  >  0) vl_covdet_set_calc_inverse_sm (det, calc_inv_sm);
      }

      if (verbose) {
        printf("covdet: filter settings:\n") ;
        printf("covdet:   method                = %s\n",
               vlCovdetMethods[method].name) ;
        printf("covdet:   output frame type     = %s\n",
               vlFrameTypes[frm_type].name) ;
        printf("covdet:   octaves      (O)      = %d\n",
               (int)vl_scalespace_get_octaves_num   (gss)) ;
        printf("covdet:   levels       (S)      = %d\n",
               (int)vl_scalespace_get_levels_num    (gss)) ;
        printf("covdet:   first octave (o_min)  = %d\n",
               (int)vl_scalespace_get_octave_min    (gss)) ;
        printf("covdet:   edge thresh           = %g\n",
               vl_covdet_get_edge_thresh      (det)) ;
        printf("covdet:   peak thresh           = %g\n",
               vl_covdet_get_peak_thresh      (det)) ;
        if (calc_descs) {
          printf("covdet:   norm thresh           = %g\n",
                    vl_covdet_get_norm_thresh     (det)) ;
          printf("covdet:   magnif                = %g\n",
                    vl_covdet_get_magnif          (det)) ;
          if (!calc_affine) {
            printf("covdet:   window size           = %g\n",
                      vl_covdet_get_window_size     (det)) ;
          }
        }
        if (calc_affine) {
          printf("covdet:   aff_win_size          = %d\n",
                    (int)vl_affineshapeestimator_get_window_size  (aff_det)) ;
          printf("covdet:   max_iter              = %d\n",
                    (int)vl_affineshapeestimator_get_max_iter     (aff_det)) ;
          printf("covdet:   conv_thr              = %g\n",
                    vl_affineshapeestimator_get_conv_thresh  (aff_det)) ;
          printf("covdet:   patch_size            = %d\n",
                    (int)vl_affinepatchnormalizer_get_patch_size (aff_norm)) ;
          printf("covdet:   calc_inverse_sm       = %s\n",
                 vl_covdet_get_calc_inverse_sm(det) ? "yes" : "no") ;
        }

        printf("covdet:   calc_decriptors       = %s\n",
               calc_descs ? "yes" : "no") ;

        printf((nframes > 0) ?
                  "covdet: will source frames? yes (%d read)\n" :
                  "covdet: will source frames? no\n", nframes) ;
      }

      /* ...............................................................
       *                                                         Do job
       * ............................................................ */
      {
        void const *current_frame, *next_frame = 0;

        /* run detector ............................................. */
        if (nframes == 0) {
          vl_covdet_detect(det, fdata);

          frames  = vl_covdet_get_frames_storage (det) ;
          nframes = vl_covdet_get_frames_num (det) ;
          descrs = vl_covdet_get_descriptors (det);

          if (verbose > 1) {
            printf ("covdet: detected %d '%s' frames\n", nframes,
                    vlFrameNames[frm_type]) ;
          }
        } else {
          /* TODO solve frames import */
          vl_covdet_convert_frames(det, fdata, inputFrames, numInputFrames,
                                   inputFrameType);

          frames  = vl_covdet_get_frames_storage (det) ;
          nframes = vl_covdet_get_frames_num (det) ;
          descrs = vl_covdet_get_descriptors (det);
        }

        /* optionally save GSS */
        if (gsf.active) {
          err = save_gss (gss, &gsf, basename, verbose, VL_FALSE) ;
          if (err) {
            snprintf (err_msg, sizeof(err_msg),
                      "Could not write GSS to PGM file.") ;
            goto done ;
          }
        }

        current_frame = frames;

        /* For each keypoint ........................................ */
        for (i = 0; i < nframes ; ++i) {

          if (out.active) {
            unsigned int l ;
            next_frame = vl_file_meta_put_frame(&out, current_frame, frm_type);

            if (calc_descs) {
              for (l = 0 ; l < desc_sz ; ++l) {
                vl_file_meta_put_uint8 (&out, (vl_uint8) (512.0 * descrs [l])) ;
              }
            }
            if (out.protocol == VL_PROT_ASCII) fprintf(out.file, "\n") ;
          }


          if (frm.active) {
            next_frame = vl_file_meta_put_frame(&frm, current_frame, frm_type);
            if (frm.protocol == VL_PROT_ASCII) fprintf(frm.file, "\n") ;
          }

          if (calc_descs && dsc.active) {
            unsigned int l ;
            for (l = 0 ; l < 128 ; ++l) {
              double x = 512.0 * descrs[l] ;
              x = (x < 255.0) ? x : 255.0 ;
              vl_file_meta_put_uint8 (&dsc, (vl_uint8) (x)) ;
            }
            if (dsc.protocol == VL_PROT_ASCII) fprintf(dsc.file, "\n") ;
          }

          if (calc_descs){
            descrs += desc_sz;
          }

          current_frame = next_frame;
        } /* next keypoint */
      }
    }

    /* ...............................................................
     *                                                       Finish up
     * ............................................................ */

    if (met.active) {
      fprintf(met.file, "<covdet\n") ;
      fprintf(met.file, "  input       = '%s'\n", name) ;
      if (dsc.active) {
        fprintf(met.file, "  descriptors = '%s'\n", dsc.name) ;
      }
      if (frm.active) {
        fprintf(met.file,"  frames      = '%s'\n", frm.name) ;
      }
      fprintf(met.file, ">\n") ;
    }

  done :
    /* release input keys buffer */
    if (inputFrames) {
      free (inputFrames) ;
      inputFramesStorageSize = numInputFrames = 0 ;
      inputFrames = 0;
    }

    /* release filter */
    if (det) {
      vl_covdet_delete (det) ;
      det = 0 ;
    }

    /* release image data */
    if (fdata) {
      free (fdata) ;
      fdata = 0 ;
    }

    /* release image data */
    if (data) {
      free (data) ;
      data = 0 ;
    }

    /* close files */
    if (in) {
      fclose (in) ;
      in = 0 ;
    }

    vl_file_meta_close (&out) ;
    vl_file_meta_close (&frm) ;
    vl_file_meta_close (&dsc) ;
    vl_file_meta_close (&met) ;
    vl_file_meta_close (&gsf) ;
    vl_file_meta_close (&ifr) ;

    /* if bad print error message */
    if (err) {
      fprintf
        (stderr,
         "covdet: err: %s (%d)\n",
         err_msg,
         err) ;
      exit_code = 1 ;
    }
  }

  /* quit */
  return exit_code ;
}
