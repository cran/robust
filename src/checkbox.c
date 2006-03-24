#include <S.h>
#include "robust.h"
#include <ctype.h>

void s_whitmm(char **vec, Sint *nv, Sint *n)
{
  Sint i, j;
  for (i=0; i<*nv; i++)
    for (j=0; j<*n; j++)
      vec[i][j] = ' ';
}
  
void s_ctonmm(char **vec, Sint *nv, Sint *n, char **ovec, Sint *nout)
{
  /* vec : nv x n, input  */
  /* ovec: n  x n, output */
  Sint i, j;
  char c;
  Sint preIsDigit = 0, cpos = 0;

  *nout = -1;
  for (i=0; i<*nv; i++) {
    for (j=0; j<*n; j++) {
      if (isdigit(c=vec[i][j])) {
        if (!preIsDigit) {
          preIsDigit = 1;
          ++(*nout);
          cpos = 0;
          ovec[*nout][cpos] = c;
        }
        else {
          ++cpos;
          ovec[*nout][cpos] = c;
        }
      }
      else
        preIsDigit = 0;
    }
  }
}
