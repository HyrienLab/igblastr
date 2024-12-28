library(testthat)
library(igblastr)

is_bioc_build_machine <- isTRUE(as.logical(Sys.getenv("IS_BIOC_BUILD_MACHINE")))
if (!is_bioc_build_machine) {
    ## We temporarily set the cache to a different (non-persistent)
    ## location to prevent the tests from messing up the real cache
    ## located at 'R_user_dir("igblastr", "cache")'.
    test_cache <- file.path(tempdir(), "igblastr_test_cache")
    options(igblastr_cache=test_cache)
}

if (!has_igblast()) install_igblast()
test_check("igblastr")

if (!is_bioc_build_machine) {
    options(igblastr_cache=NULL)
}
