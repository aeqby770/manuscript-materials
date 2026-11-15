// LineaLinkR · src/entrypoint.c
//
// Forwards R_init_LineaLinkR to extendr's auto-generated registrator.
// This .c file exists solely so the linker is forced to materialise a
// call site for R_init_LineaLinkR_extendr and therefore pulls the
// registration symbol out of the static Rust library.  Without this
// forward, GNU ld strips R_init_LineaLinkR_extendr as unused and R
// fails to locate any .Call entry points at package load time.
//
// The Rust side lives in
// src/extendr_bindings.rs → `extendr_module! { mod LineaLinkR; ... }`.

void R_init_LineaLinkR_extendr(void *dll);

void R_init_LineaLinkR(void *dll) {
    R_init_LineaLinkR_extendr(dll);
}
