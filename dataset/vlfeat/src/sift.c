/** @internal @file sift.c
 ** @author Andrea Vedaldi
 ** @author Karel Lenc
 ** @brief Scale Invariant Feature Transform (SIFT) - Driver
 **/

/*
Copyright (C) 2007-12 Andrea Vedaldi and Brian Fulkerson.
All rights reserved.

This file is part of the VLFeat library and is made available under
the terms of the BSD license (see the COPYING file).
*/

#define VL_SIFT_DRIVER_VERSION 0.2

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
  " --verbose -v      Be verbose\n"
  " --help -h         Print this help message\n"
  " --output -o       Specify output file\n"
  " --frames          Specify frames file\n"
  " --descriptors     Specify descriptors file\n"
  " --meta            Specify meta file\n"
  " --gss             Specify Gaussian scale space files\n"
  " --octaves -O      Number of octaves\n"
  " --levels -S       Number of levels per octave\n"
  " --first-octave    Index of the first octave\n"
  " --edge-thresh     Specify the edge threshold\n"
  " --peak-thresh     Specift the peak threshold\n"
  " --magnif          Specify the magnification factor\n"
  " --no-descriptors  Do not calculate descriptors\n"
  " --read-frames     Specify a file from which to read frames\n"
  " --orientations    Force the computation of the orientations\n"
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
  opt_peakThreshold,
  opt_magnif,
  opt_no_descr,
  opt_read_frames,
  opt_orientations
} ;

/* short options */
char const opts [] = "vhO:S:o:" ;

/* long options */
struct option const longopts [] = {
  { "verbose",         no_argument,            0, 'v'              },
  { "help",            no_argument,            0, 'h'              },
  { "octaves",         required_argument,      0, 'O'              },
  { "levels",          required_argument,      0, 'S'              },
  { "output",          optional_argument,      0, 'o'              },
  { "meta",            optional_argument,      0, opt_meta         },
  { "frames",          optional_argument,      0, opt_frames       },
  { "descriptors",     optional_argument,      0, opt_descriptors  },
  { "gss",             optional_argument,      0, opt_gss          },
  { "first-octave",    required_argument,      0, opt_first_octave },
  { "edge-thresh",     required_argument,      0, opt_edge_thresh  },
  { "peak-thresh",     required_argument,      0, opt_peakThreshold  },
  { "magnif",          required_argument,      0, opt_magnif       },
  { "no-descriptors",  no_argument,            0, opt_no_descr     },
  { "read-frames",     required_argument,      0, opt_read_frames  },
  { "orientations",    no_argument,            0, opt_orientations },
  { 0,                 0,                      0, 0                }
} ;


/* ----------------------------------------------------------------- */
/** @brief Keypoint ordering
 ** @internal
 **/
int
korder (void const* a, void const* b) {
  double x = ((double*) a) [2] - ((double*) b) [2] ;
  if (x < 0) return -1 ;
  if (x > 0) return +1 ;
  return 0 ;
}

/* ---------------------------------------------------------------- */
/** @brief SIFT driver entry point
 **/
int
main(int argc, char **argv)
{
  double edgeThreshold = -1 ;
  double peakThreshold = -1 ;
  double magnif = -1 ;
  int numOctaves = -1, numLevels = 3, firstOctave = -1 ;

  int  verbose = 0 ;
  vl_bool err = VL_ERR_OK ;
  char err_msg [1024] ;
  int n ;
  int exit_code = 0 ;
  vl_bool force_output = VL_FALSE ;
  vl_bool force_orientations = VL_FALSE ;
  vl_bool calc_desc = VL_TRUE ;

  VlFileMeta out  = {1, "%.sift",  VL_PROT_ASCII, "", 0} ;
  VlFileMeta frm  = {0, "%.frame", VL_PROT_ASCII, "", 0} ;
  VlFileMeta dsc  = {0, "%.descr", VL_PROT_ASCII, "", 0} ;
  VlFileMeta met  = {0, "%.meta",  VL_PROT_ASCII, "", 0} ;
  VlFileMeta gss  = {0, "%.pgm",   VL_PROT_ASCII, "", 0} ;
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
    if (ch == -1 && argc - optind == 0) {
      ch = 'h';
    }

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
      printf ("GSS          filespec: '%s'\n", gss.pattern) ;
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

    case opt_gss :
      /* --gss .................................................... */
      err = vl_file_meta_parse (&gss, optarg) ;
      if (err)
        ERRF("The arguments of '%s' is invalid.", argv [optind - 1]) ;
      break ;

    case 'O' :
      /* --octaves ............................................... */
      n = sscanf (optarg, "%d", &numOctaves) ;
      if (n == 0 || numOctaves < 0)
        ERRF("The argument of '%s' must be a non-negative integer.",
            argv [optind - 1]) ;
      break ;

    case 'S' :
      /* --levels ............................................... */
      n = sscanf (optarg, "%d", &numLevels) ;
      if (n == 0 || numLevels < 0)
        ERRF("The argument of '%s' must be a non-negative integer.",
            argv [optind - 1]) ;
      break ;

    case opt_first_octave :
      /* --first-octave ......................................... */
      n = sscanf (optarg, "%d", &firstOctave) ;
      if (n == 0)
        ERRF("The argument of '%s' must be an integer.",
            argv [optind - 1]) ;
      break ;

    case opt_edge_thresh :
      /* --edge-thresh ........................................... */
      n = sscanf (optarg, "%lf", &edgeThreshold) ;
      if (n == 0 || edgeThreshold < 1)
        ERRF("The argument of '%s' must be not smaller than 1.",
            argv [optind - 1]) ;
      break ;

    case opt_peakThreshold :
      /* --edge-thresh ........................................... */
      n = sscanf (optarg, "%lf", &peakThreshold) ;
      if (n == 0 || peakThreshold < 0)
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

    case opt_no_descr :
      /* --no-descriptors ....................................... */
      calc_desc = 0 ;
      break ;

    case opt_orientations :
      /* --orientations ......................................... */
      force_orientations = 1 ;
      break ;

    case 0 :
    default :
      /* should not get here ...................................... */
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

  if (verbose > 1) {
#define PRNFO(name,fm)                                                  \
    printf("sift: " name) ;                                             \
    printf("%3s ",  (fm).active ? "yes" : "no") ;                       \
    printf("%-6s ", vl_string_protocol_name ((fm).protocol)) ;          \
    printf("%-10s\n", (fm).pattern) ;

    PRNFO("write aggregate . ", out) ;
    PRNFO("write frames .... ", frm) ;
    PRNFO("write descriptors ", dsc) ;
    PRNFO("write meta ...... ", met) ;
    PRNFO("write GSS ....... ", gss) ;
    PRNFO("read  frames .... ", ifr) ;

    if (force_orientations)
      printf("sift: will compute orientations\n") ;
  }

  /* ------------------------------------------------------------------
   *                                         Process one image per time
   * --------------------------------------------------------------- */

  while (argc--) {

    char             basename [1024] ;
    char const      *name = *argv++ ;

    FILE            *in    = 0 ;
    vl_uint8        *data  = 0 ;
    float     *fdata = 0 ;
    VlPgmImage       pim ;

    VlCovDet      *det = 0;
    vl_size          q ;
    int              i ;

    double           *ikeys = 0 ;
    int              nikeys = 0, ikeys_size = 0 ;

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
      printf ("sift: <== '%s'\n", name) ;
    }

    if (verbose > 1) {
      printf ("sift: basename is '%s'\n", basename) ;
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
      printf ("sift: image is %" VL_FMT_SIZE " by %" VL_FMT_SIZE " pixels\n",
              pim. width,
              pim. height) ;

    /* allocate buffer */
    data  = malloc(vl_pgm_get_npixels (&pim) *
                   vl_pgm_get_bpp (&pim) * sizeof(vl_uint8)   ) ;
    fdata = malloc(vl_pgm_get_npixels (&pim) *
                   vl_pgm_get_bpp (&pim) * sizeof(float)) ;

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
     *                                     Optionally source keypoints
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

    if (ifr.active) {

      /* open file */
      err = vl_file_meta_open (&ifr, basename, "rb") ;
      WERR(ifr.name, reading) ;

#define QERR                                                            \
      if (err ) {                                                       \
        snprintf (err_msg, sizeof(err_msg),                             \
                  "'%s' malformed", ifr.name) ;                         \
        err = VL_ERR_IO ;                                               \
        goto done ;                                                     \
      }

      while (1) {
        double x, y, s, th ;

        /* read next guy */
        err = vl_file_meta_get_double (&ifr, &x) ;
        if   (err == VL_ERR_EOF) break;
        else QERR ;
        err = vl_file_meta_get_double (&ifr, &y ) ; QERR ;
        err = vl_file_meta_get_double (&ifr, &s ) ; QERR ;
        err = vl_file_meta_get_double (&ifr, &th) ;
        if   (err == VL_ERR_EOF) break;
        else QERR ;

        /* make enough space */
        if (ikeys_size < nikeys + 1) {
          ikeys_size += 10000 ;
          ikeys       = realloc (ikeys, 4 * sizeof(double) * ikeys_size) ;
        }

        /* add the guy to the buffer */
        ikeys [4 * nikeys + 0]  = x ;
        ikeys [4 * nikeys + 1]  = y ;
        ikeys [4 * nikeys + 2]  = s ;
        ikeys [4 * nikeys + 3]  = th ;

        ++ nikeys ;
      }

      /* now order by scale */
      qsort (ikeys, nikeys, 4 * sizeof(double), korder) ;

      if (verbose) {
        printf ("sift: read %d keypoints from '%s'\n", nikeys, ifr.name) ;
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
      if (out.active) printf("sift: writing all ....... to . '%s'\n", out.name);
      if (frm.active) printf("sift: writing frames .... to . '%s'\n", frm.name);
      if (dsc.active) printf("sift: writing descriptors to . '%s'\n", dsc.name);
      if (met.active) printf("sift: writing meta ...... to . '%s'\n", met.name);
    }

    /* ...............................................................
     *                                                   Make detector
     * ............................................................ */

    {
      VlScaleSpace *gss_scsp;
      VlScaleSpace *grad;
      VlScaleSpace const *resp;
      float const *descrs = 0;

      /*det = vl_siftdet_new (pim.width, pim.height, O, S, omin,
                            VL_TRUE, calc_desc) ;*/
      det = vl_covdet_new_disc_detector(pim.width, pim.height, VL_COVDET_METHOD_DOG,
                                          numOctaves, firstOctave, numLevels, VL_TRUE, VL_TRUE);
      gss_scsp = vl_covdet_get_gss(det);
      grad = vl_covdet_get_gradient(det);
      resp = vl_covdet_get_responses(det);

      if (edgeThreshold >= 0) vl_covdet_set_edge_thresh(det, edgeThreshold) ;
      if (peakThreshold >= 0) vl_covdet_set_peak_thresh(det, peakThreshold) ;
      if (magnif >= 0) vl_covdet_set_magnif(det, magnif) ;

      if (!det) {
        snprintf (err_msg, sizeof(err_msg), "Could not create a covariant frames detector.") ;
        err = VL_ERR_ALLOC ;
        goto done ;
      }

      if (verbose > 1) {
        printf("vl_sift: num octaves = %d\n", (int)vl_scalespace_get_octaves_num(gss_scsp)) ;
        printf("vl_sift: num levels per octave = %d\n", (int)vl_scalespace_get_levels_num(gss_scsp)) ;
        printf("vl_sift: first octave = %d\n", (int)vl_scalespace_get_octave_min(gss_scsp)) ;
        printf("vl_sift: edge threshold = %g\n", vl_covdet_get_edge_thresh(det)) ;
        printf("vl_sift: peak threshold = %g\n", vl_covdet_get_peak_thresh(det)) ;
        printf("vl_sift: norm threshold = %g\n", vl_covdet_get_norm_thresh(det)) ;
        printf("vl_sift: magnif = %g\n", vl_covdet_get_magnif(det)) ;
        printf("vl_sift: window size = %g\n", vl_covdet_get_window_size(det)) ;
        printf("vl_sift: will calculate descriptors? %s\n", VL_YESNO(calc_desc)) ;
        printf("vl_sift: will source frames? %s\n", VL_YESNO(ikeys)) ;
        printf("vl_sift: will force orientations? %s\n", VL_YESNO(force_orientations)) ;
      }

      /* .............................................................
       *                                            Process each octave
       * .......................................................... */
      i = 0 ;
      {
        VlFrameOrientedDisc const *keys = 0 ;
        int                         nkeys ;

        /* run detector ............................................. */
        if (ikeys == 0) {
          vl_covdet_detect(det, fdata);

          keys = vl_covdet_get_oriented_discs(det) ;
          nkeys = vl_covdet_get_frames_num(det) ;
          descrs = vl_covdet_get_descriptors(det) ;

          i = 0 ;

          if (verbose > 1) {
            printf ("sift: detected %d (oriented) keypoints\n", nkeys) ;
          }
        } else {
          vl_scalespace_init(gss_scsp, fdata);
          vl_scalespace_apply(gss_scsp, grad, vl_imgradient_polar_f_callback, 0);
          nkeys = nikeys ;
        }

        /* optionally save GSS */
        if (gss.active) {
          err = save_gss (gss_scsp, &gss, basename, verbose, VL_FALSE) ;
          if (err) {
            snprintf (err_msg, sizeof(err_msg),
                      "Could not write GSS to PGM file.") ;
            goto done ;
          }
        }

        /* for each keypoint ........................................ */
        for (; i < nkeys ; ++i) {
          double                angles [4] ;
          int                   nangles ;
          VlFrameOrientedDisc const *k = 0 ;
          VlFrameOrientedDisc        tk ;
          VlScaleSpaceFrame    ik ;

          /* obtain keypoint orientations ........................... */
          if (ikeys) {
            vl_scalespace_frame_init (resp, &ik,
                                   ikeys [4 * i + 0],
                                   ikeys [4 * i + 1],
                                   ikeys [4 * i + 2]) ;

            tk.x = ik.x;
            tk.y = ik.y;
            tk.sigma = ik.sigma;
            k = &tk;

            /* optionally compute orientations too */
            if (force_orientations) {
              nangles = vl_covdet_calc_ssframe_orientations
                        (grad, angles, &ik) ;
            } else {
              angles [0] = ikeys [4 * i + 3] ;
              nangles    = 1 ;
            }
          } else {
            k = keys + i;
            angles[0] = k->angle;
            nangles = 1;
          }

          /* for each orientation ................................... */
          for (q = 0 ; q < (unsigned) nangles ; ++q) {
            float descr_buf [128] ;
            float const *descr = 0;

            /* compute descriptor (if necessary) */
            if (calc_desc && (out.active || dsc.active)) {
              if (ikeys) {
                vl_covdet_calc_ssframe_descriptor
                    (det, descr_buf, &ik, angles [q]) ;
                descr = descr_buf;
              } else {
                descr = descrs + i * 128;
              }
            }

            if (out.active) {
              int l ;
              vl_file_meta_put_double (&out, k -> x     ) ;
              vl_file_meta_put_double (&out, k -> y     ) ;
              vl_file_meta_put_double (&out, k -> sigma ) ;
              vl_file_meta_put_double (&out, angles [q] ) ;
              if (calc_desc) {
                for (l = 0 ; l < 128 ; ++l) {
                  vl_file_meta_put_uint8 (&out, (vl_uint8) (512.0 * descr [l])) ;
                }
              }
              if (out.protocol == VL_PROT_ASCII) fprintf(out.file, "\n") ;
            }

            if (frm.active) {
              vl_file_meta_put_double (&frm, k -> x     ) ;
              vl_file_meta_put_double (&frm, k -> y     ) ;
              vl_file_meta_put_double (&frm, k -> sigma ) ;
              vl_file_meta_put_double (&frm, angles [q] ) ;
              if (frm.protocol == VL_PROT_ASCII) fprintf(frm.file, "\n") ;
            }

            if (calc_desc && dsc.active) {
              int l ;
              for (l = 0 ; l < 128 ; ++l) {
                double x = 512.0 * descr[l] ;
                x = (x < 255.0) ? x : 255.0 ;
                vl_file_meta_put_uint8 (&dsc, (vl_uint8) (x)) ;
              }
              if (dsc.protocol == VL_PROT_ASCII) fprintf(dsc.file, "\n") ;
            }
          }
        }
      }
    }

    /* ...............................................................
     *                                                       Finish up
     * ............................................................ */

    if (met.active) {
      fprintf(met.file, "<sift\n") ;
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
    if (ikeys) {
      free (ikeys) ;
      ikeys_size = nikeys = 0 ;
      ikeys = 0 ;
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
    vl_file_meta_close (&gss) ;
    vl_file_meta_close (&ifr) ;

    /* if bad print error message */
    if (err) {
      fprintf
        (stderr,
         "sift: err: %s (%d)\n",
         err_msg,
         err) ;
      exit_code = 1 ;
    }
  }

  /* quit */
  return exit_code ;
}
