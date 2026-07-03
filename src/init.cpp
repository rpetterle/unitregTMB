#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

extern "C" void tmb_custom_init(DllInfo *dll);

extern "C" void R_init_unitregTMB(DllInfo *dll) {
    tmb_custom_init(dll);
    R_useDynamicSymbols(dll, TRUE);
}
