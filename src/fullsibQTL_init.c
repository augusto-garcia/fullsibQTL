#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void onemap_calc_genoprob_outbred(void *, void *, void *, void *, void *, void *, void *, void *);
extern void R_char_qtl(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void R_scan_qtl(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"onemap_calc_genoprob_outbred", (DL_FUNC) &onemap_calc_genoprob_outbred,  8},
    {"R_char_qtl",                   (DL_FUNC) &R_char_qtl,                   13},
    {"R_scan_qtl",                   (DL_FUNC) &R_scan_qtl,                   14},
    {NULL, NULL, 0}
};

void R_init_fullsibQTL(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}