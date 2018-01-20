int_string <- function(k)
{
  ndig <- nchar(as.character(k))
  ret <- as.character(1:k)

  paste_mod <- function(...)
  {
    paste(sep="", ...)
  }

  pad <- function(v)
  {
    l <- nchar(v)
    do.call("paste_mod", as.list(c(rep("0", ndig-l), v)))
  }

  sapply(ret, pad)
}
