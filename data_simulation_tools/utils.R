create_named_list_ <- function(...) {
  setNames(list(...), as.character(match.call()[-1]))
}
