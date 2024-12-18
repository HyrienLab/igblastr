.onLoad <- function(libname, pkgname)
{
    create_all_c_region_dbs(force=TRUE)
    clean_all_germline_dbs()
}

