#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
  Check these declarations against the C/Fortran source code.
*/
  
  /* .C calls */
  extern void Rbclcov(void *, void *, void *, void *, void *, void *, void *);
extern void Rm2(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void Rm2rasch(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void Rnrbcpln(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void Rnrmlepln(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void Rnrmlerasch(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void Rsimulpln(void *, void *, void *, void *, void *, void *);
extern void Rstartpln(void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
  {"Rbclcov",     (DL_FUNC) &Rbclcov,      7},
  {"Rm2",         (DL_FUNC) &Rm2,         10},
  {"Rm2rasch",    (DL_FUNC) &Rm2rasch,    10},
  {"Rnrbcpln",    (DL_FUNC) &Rnrbcpln,    14},
  {"Rnrmlepln",   (DL_FUNC) &Rnrmlepln,   16},
  {"Rnrmlerasch", (DL_FUNC) &Rnrmlerasch, 16},
  {"Rsimulpln",   (DL_FUNC) &Rsimulpln,    6},
  {"Rstartpln",   (DL_FUNC) &Rstartpln,    5},
  {NULL, NULL, 0}
};

void R_init_pln(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
