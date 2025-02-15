### =========================================================================
### Low-level cache utilities
### -------------------------------------------------------------------------
###
### Nothing in this file is exported.
###
### The igblastr cache is divided into various compartments that are
### completely independent one from each other. These compartments are:
###
###   compartment       path
###   ----------------  ------------------------------------
###   IGBLAST_ROOTS     <igblastr-cache>/igblast_roots
###   GERMLINE_DBS      <igblastr-cache>/germline_dbs
###   C_REGION_DBS      <igblastr-cache>/c_region_dbs
###   IMGT_LOCAL_STORE  <igblastr-cache>/store/IMGT-releases
###
### Always use the igblastr_cache() function implemented in this file to
### construct these paths.

IGBLAST_ROOTS    <- "IGBLAST_ROOTS"
GERMLINE_DBS     <- "GERMLINE_DBS"
C_REGION_DBS     <- "C_REGION_DBS"
IMGT_LOCAL_STORE <- "IMGT_LOCAL_STORE"


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### igblastr_cache()
###

.IGBLASTR_CACHES <- c(IGBLAST_ROOTS, GERMLINE_DBS,
                      C_REGION_DBS, IMGT_LOCAL_STORE)

### Returns the path to the cache compartment specified via 'which'.
igblastr_cache <- function(which=NULL)
{
    path <- getOption("igblastr_cache", R_user_dir("igblastr", "cache"))
    if (is.null(which))
        return(path)
    stopifnot(isSingleNonWhiteString(which))
    switch(which,
        IGBLAST_ROOTS   =file.path(path, "igblast_roots"),
        GERMLINE_DBS    =file.path(path, "germline_dbs"),
        C_REGION_DBS    =file.path(path, "c_region_dbs"),
        IMGT_LOCAL_STORE=file.path(path, "store", "IMGT-releases"),
	stop(wmsg("'which' must be one of ",
                  paste0("\"", .IGBLASTR_CACHES, "\"", collapse=", ")))
    )
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### get_internal_igblast_roots()
###

get_internal_igblast_roots <- function() igblastr_cache(IGBLAST_ROOTS)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .reset_igblastr_cache()
###

.reset_igblastr_cache <- function(which=NULL)
{
    path <- igblastr_cache(which)
    nuke_file(path)
}

reset_germline_dbs_cache <- function() .reset_igblastr_cache(GERMLINE_DBS)
reset_c_region_dbs_cache <- function() .reset_igblastr_cache(C_REGION_DBS)

