# Symbols used by the non-standard evaluation of the 'meta' package:
# metagen(TE = TE, seTE = seTE, data = meta_data, studlab = study) refers to
# columns of meta_data by name, which codetools cannot see.
utils::globalVariables(c("TE", "seTE", "study"))
