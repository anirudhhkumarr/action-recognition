/** @author   Andrea Vedaldi
 ** @brief    Support for command line drivers - Definition.
 ** @internal
 **
 ** This file contains support code which is shared by the command
 ** line drivers.
 **/

/*
Copyright (C) 2007-12 Andrea Vedaldi and Brian Fulkerson.
All rights reserved.

This file is part of the VLFeat library and is made available under
the terms of the BSD license (see the COPYING file).
*/

#ifndef VL_SCALESPACE_EXPORT
#define VL_SCALESPACE_EXPORT

#include "generic-driver.h"

#include <vl/generic.h>
#include <vl/pgm.h>
#include <vl/stringop.h>
#include <vl/scalespace.h>

#include <stdio.h>
#include <assert.h>
#include <float.h>

/* ----------------------------------------------------------------- */
/** @brief Save scale space on disk
 ** @internal
 **/
int
save_gss (VlScaleSpace const *gss, VlFileMeta * fm, const char * basename,
          int verbose, vl_bool normalise)
{
  char tmp [1024] ;
  int i ;
  int s, err = 0 ;
  int w, h ;
  int o_min = vl_scalespace_get_octave_min(gss);
  int o_max = vl_scalespace_get_octave_max(gss);
  int s_min = vl_scalespace_get_level_min(gss);
  int s_max = vl_scalespace_get_level_max(gss);
  int o;
  VlPgmImage pim ;
  vl_uint8 *buffer = 0 ;
  vl_size q ;

  float max = FLT_MIN;
  float min = FLT_MAX;
  float norm;
  float shift;

  for(o = o_min; o <= o_max; ++o){
    if (! fm -> active) {
      return VL_ERR_OK ;
    }

    w = vl_scalespace_get_octave_width  (gss, o) ;
    h = vl_scalespace_get_octave_height (gss, o) ;

    pim.width     = w ;
    pim.height    = h ;
    pim.max_value = 255 ;
    pim.is_raw    = 1 ;

    buffer = malloc (sizeof(vl_uint8) * w * h) ;
    if (! buffer) {
      err = VL_ERR_ALLOC ;
      break;
    }

    q = vl_string_copy (tmp, sizeof(tmp), basename) ;
    if (q >= sizeof(tmp)) {
      err = VL_ERR_OVERFLOW ;
      break;
    }

    for (s = s_min ; s <= s_max ; ++s) {
      float const *pt = vl_scalespace_get_octave (gss, o, s) ;

      /* conversion */
      if (normalise) {
        for(i = 0; i < w*h; ++i){
          if (pt[i] > max) max = pt[i];
          if (pt[i] < min) min = pt[i];
        }
        norm = 1./(max - min);
        shift = - min * norm;

        for (i = 0 ; i < w * h ; ++i) {
          buffer [i] = (vl_uint8) ((pt [i] * norm + shift) * 255.) ;
        }
      } else {
        for (i = 0 ; i < w * h ; ++i) {
          buffer [i] = (vl_uint8) (pt [i]) ;
        }
      }

      /* save */
      snprintf(tmp + q, sizeof(tmp) - q, "_%02d_%03d", o, s) ;

      err = vl_file_meta_open (fm, tmp, "wb") ;
      if (err) break ;

      err = vl_pgm_insert (fm -> file, &pim, buffer) ;
      if (err) break ;

      if (verbose) {
        printf("sift: saved gss level to '%s'\n", fm -> name) ;
      }

      vl_file_meta_close (fm) ;
    }
  }

  if (buffer) free (buffer) ;
  vl_file_meta_close (fm) ;
  return err ;
}



/* VL_SCALESPACE_EXPORT */
#endif
