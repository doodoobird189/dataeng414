
# Create targets `programmer_pdf` and `coding_convention_pdf`
include(UseLATEX)
add_latex_document(programmer.tex INPUTS classes.eps DEFAULT_PDF MANGLE_TARGET_NAMES)
add_latex_document(coding_convention.tex DEFAULT_PDF MANGLE_TARGET_NAMES)
