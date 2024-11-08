### =========================================================================
### Low-level utilities for handling NCBI precompiled IgBLAST
### -------------------------------------------------------------------------
###
### Nothing in this file is exported.


PRECOMPILED_NCBI_IGBLAST_PREFIX <- "ncbi-igblast-"

.get_precompiled_ncbi_igblast_pattern <- function()
    sprintf("^(%s([.0-9]+)).*$", PRECOMPILED_NCBI_IGBLAST_PREFIX)

### 'name' must be the name of an IgBLAST tarball (i.e. *.tar.gz file)
### or *.dmg file from NCBI FTP site e.g. "ncbi-igblast-1.22.0+.dmg".
infer_igblast_version_from_ncbi_name <- function(name)
{
    stopifnot(isSingleNonWhiteString(name))
    pattern <- .get_precompiled_ncbi_igblast_pattern()
    sub(pattern, "\\2", name)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### extract_igblast_tarball()
###

### 'name' must be the name of an IgBLAST tarball (i.e. *.tar.gz file)
### from NCBI FTP site e.g. "ncbi-igblast-1.22.0-x64-win64.tar.gz".
.infer_igblast_rootbasename_from_ncbi_name <- function(name)
{
    stopifnot(isSingleNonWhiteString(name))
    pattern <- .get_precompiled_ncbi_igblast_pattern()
    sub(pattern, "\\1", name)
}

### 'destdir' should be the path to an existing directory.
### Extracts and installs IgBLAST in <destdir>/<version>/
extract_igblast_tarball <- function(tarfile, ncbi_name, destdir=".")
{
    ## Create <destdir>/tmpexdir/
    tempexdir <- file.path(destdir, "tmpexdir")
    unlink(tempexdir, recursive=TRUE, force=TRUE)
    dir.create(tempexdir)
    on.exit(unlink(tempexdir, recursive=TRUE, force=TRUE))

    ## Untar in <destdir>/tmpexdir/
    untar2(tarfile, ncbi_name, exdir=tempexdir)

    ## Get IgBLAST rootbasename and version.
    rootbasename <- .infer_igblast_rootbasename_from_ncbi_name(ncbi_name)
    version <- infer_igblast_version_from_ncbi_name(ncbi_name)

    ## Move <destdir>/tmpexdir/<rootbasename> (newdir)
    ## to <destdir>/<version> (olddir) after nuking the latter if needed.
    olddir <- file.path(destdir, version)
    newdir <- file.path(tempexdir, rootbasename)
    replace_file(olddir, newdir)
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

### 'destdir' should be the path to an existing directory.
### Extracts and installs IgBLAST in <destdir>/<version>/
extract_igblast_dmg <- function(dmgfile, ncbi_name, destdir=".")
{
    stopifnot(isSingleNonWhiteString(dmgfile),
              isSingleNonWhiteString(ncbi_name),
              isSingleNonWhiteString(destdir),
              dir.exists(destdir))

    version <- infer_igblast_version_from_ncbi_name(ncbi_name)
    dmg_mounting_point_name <- sub("\\.dmg$", "", ncbi_name)
    dmg_mounting_point <- file.path("/Volumes", dmg_mounting_point_name)
    pkg_basename <- paste0(dmg_mounting_point_name, ".pkg")
    pkg_path <- file.path(dmg_mounting_point, pkg_basename)
    expand_dir <- file.path(destdir, dmg_mounting_point_name)

    .attach_dmg(dmgfile)
    on.exit(.detach_dmg(dmg_mounting_point))
    .full_expand_pkgfile(pkg_path, expand_dir)
    olddir <- file.path(destdir, version)
    newdir <- file.path(expand_dir, "binaries.pkg", "Payload")
    replace_file(olddir, newdir)
    unlink(expand_dir, recursive=TRUE, force=TRUE)
}

