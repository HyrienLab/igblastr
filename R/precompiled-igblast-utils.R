### =========================================================================
### Low-level utilities for handling NCBI precompiled IgBlast
### -------------------------------------------------------------------------
###
### Nothing in this file is exported.


PRECOMPILED_NCBI_IGBLAST_PREFIX <- "ncbi-igblast-"

.get_precompiled_ncbi_igblast_pattern <- function()
    sprintf("^(%s([.0-9]+)).*$", PRECOMPILED_NCBI_IGBLAST_PREFIX)

### 'name' must be the name of an IgBlast tarball (i.e. *.tar.gz file)
### from NCBI FTP site e.g. "ncbi-igblast-1.22.0-x64-win64.tar.gz".
infer_rootbasename_from_ncbi_igblast_tarball_name <- function(name)
{
    stopifnot(isSingleNonWhiteString(name))
    pattern <- .get_precompiled_ncbi_igblast_pattern()
    sub(pattern, "\\1", name)
}

### 'name' must be the name of an IgBlast tarball (i.e. *.tar.gz file)
### or *.dmg file from NCBI FTP site e.g. "ncbi-igblast-1.22.0+.dmg".
infer_version_from_ncbi_igblast_name <- function(name)
{
    stopifnot(isSingleNonWhiteString(name))
    pattern <- .get_precompiled_ncbi_igblast_pattern()
    sub(pattern, "\\2", name)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### untar2()
###

### A thin wrapper to untar() with more user-friendly error handling.
### 'exdir' should be the path to an existing directory that is
### preferrably empty.
untar2 <- function(tarfile, ncbi_name, exdir=".")
{
    stopifnot(isSingleNonWhiteString(tarfile),
              isSingleNonWhiteString(ncbi_name),
              isSingleNonWhiteString(exdir),
              dir.exists(exdir))
    code <- suppressWarnings(untar(tarfile, exdir=exdir))
    if (code != 0L)
        stop(wmsg("Anomaly: something went wrong during ",
                  "extraction of '", tarfile, "' (the local copy ",
                  "of '", ncbi_name, "') to '", exdir, "'."))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### extract_igblast_dmg()
###

.attach_dmg <- function(dmgfile)
{
    out <- system2("hdiutil", args=c("attach", dmgfile), stdout=TRUE)
}

.detach_dmg <- function(dmg_mounting_point)
{
    args <- c("detach", dmg_mounting_point)
    out <- suppressWarnings(system2("hdiutil", args=args,
                                    stdout=TRUE, stderr=TRUE))
}

.full_expand_pkgfile <- function(pkgfile, expand_dir)
{
    if (file.exists(expand_dir))
        stop(wmsg("Could not unarchive ", pkgfile, " ",
                  "to ", expand_dir, " (destination exists)"))
    args <- c("--expand-full", pkgfile, expand_dir)
    out <- suppressWarnings(system2("pkgutil", args=args,
            stdout=TRUE, stderr=TRUE))
    status <- attr(out, "status")
    if (!(is.null(status) || isTRUE(all.equal(status, 0L))))
        stop(wmsg(out))
}

### 'exdir' should be the path to an existing directory that is
### preferrably empty.
extract_igblast_dmg <- function(dmgfile, ncbi_name, exdir=".")
{
    stopifnot(isSingleNonWhiteString(dmgfile),
              isSingleNonWhiteString(ncbi_name),
              isSingleNonWhiteString(exdir),
              dir.exists(exdir))

    version <- infer_version_from_ncbi_igblast_name(ncbi_name)
    dmg_mounting_point_name <- sub("\\.dmg$", "", ncbi_name)
    dmg_mounting_point <- file.path("/Volumes", dmg_mounting_point_name)
    pkg_basename <- paste0(dmg_mounting_point_name, ".pkg")
    pkg_path <- file.path(dmg_mounting_point, pkg_basename)
    expand_dir <- file.path(exdir, dmg_mounting_point_name)

    .attach_dmg(dmgfile)
    on.exit(.detach_dmg(dmg_mounting_point))
    .full_expand_pkgfile(pkg_path, expand_dir)
    olddir <- file.path(exdir, version)
    newdir <- file.path(expand_dir, "binaries.pkg", "Payload")
    replace_file(olddir, newdir)
    unlink(expand_dir, recursive=TRUE, force=TRUE)
    olddir
}

