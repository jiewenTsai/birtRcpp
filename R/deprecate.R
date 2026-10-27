# Deprecation warnings: once per session for each deprecated function, so old scripts stay
# readable; options(birtRcpp.deprecation = "always") shows every one (for tests).

deprecation_shown <- new.env(parent = emptyenv())

deprecate_warn <- function(what, instead) {
  always <- identical(getOption("birtRcpp.deprecation"), "always")
  if (!always && isTRUE(deprecation_shown[[what]])) return(invisible(FALSE))
  deprecation_shown[[what]] <- TRUE
  text <- sprintf("%s is deprecated and gets no new features; use %s instead (shown once per session).", what, instead)
  warning(warningCondition(text, class = "deprecatedWarning"))
  invisible(TRUE)
}

deprecate_constructor <- function(constructor, model) {
  what <- sprintf("%s()", constructor)
  instead <- sprintf("ecm(data, model = \"%s\") or gibbs(data, model = \"%s\")", model, model)
  deprecate_warn(what, instead)
}
