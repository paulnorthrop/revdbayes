.onUnload <- function (libpath) {
  library.dynam.unload("revdbayes", libpath)
}
