#include <R.h>
#include <R_ext/Rdynload.h>
#include "flowFP.h"

static const R_CallMethodDef CallEntries[] = {
    {"bin_level", (DL_FUNC)&bin_level, 6},
    {"tag_events", (DL_FUNC)&tag_events, 5},
    {"count_events", (DL_FUNC)&count_events, 2},
    {NULL, NULL, 0}
};

void R_init_flowFP(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    //R_useDynamicSymbols(dll, FALSE);
}
