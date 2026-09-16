### =========================================================================
### Low-level utilities to retrieve data from the IMGT/V-QUEST download site
### -------------------------------------------------------------------------
###
### Nothing in this file is exported.
###


IMGT_URL <- "https://www.imgt.org"

### Do not remove the trailing slash.
.VQUEST_DOWNLOAD_ROOT_URL <- paste0(IMGT_URL, "/download/V-QUEST/")

### VQUEST_REFERENCE_DIRECTORY
VQUEST_REFERENCE_DIRECTORY <- "IMGT_V-QUEST_reference_directory"

.VQUEST_RELEASE_FILE_URL <-
    paste0(.VQUEST_DOWNLOAD_ROOT_URL, "IMGT_vquest_release.txt")

### Do not remove the trailing slash.
.VQUEST_ARCHIVES_URL <- paste0(.VQUEST_DOWNLOAD_ROOT_URL, "archives/")

IMGT_TERMS_OF_USE <- paste0(
    "CONDITIONS OF USE AND LICENSE: The IMGT data is provided to the ",
    "academic users and NPO's (Not for Profit Organization(s)) under ",
    "the CC BY-NC-ND 4.0 license. ",
    "See https://creativecommons.org/licenses/by-nc-nd/4.0/. ",
    "Any other use of IMGT material, from the private sector, needs ",
    "a financial arrangement with CNRS."
)

get_IMGT_connecttimeout <- function() getOption("IMGT_connecttimeout")


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### get_latest_IMGT_release()
### list_archived_IMGT_zips()
###

.IMGT_cache <- new.env(parent=emptyenv())

.fetch_latest_IMGT_release <- function()
{
    content <- getUrlContent(.VQUEST_RELEASE_FILE_URL, encoding="UTF-8",
                             connecttimeout=get_IMGT_connecttimeout())
    sub("^([^ ]*)(.*)$", "\\1", content)
}

get_latest_IMGT_release <- function(recache=FALSE)
{
    if (!isTRUEorFALSE(recache))
        stop(wmsg("'recache' must be TRUE or FALSE"))
    release <- .IMGT_cache[["LATEST_RELEASE"]]
    if (is.null(release) || recache) {
        release <- .fetch_latest_IMGT_release()
        .IMGT_cache[["LATEST_RELEASE"]] <- release
    }
    release
}

### Returns a data.frame with 3 columns (Name, Last modified, Size)
### and 1 row per .zip file.
.fetch_list_of_archived_IMGT_zips <- function()
{
    scrape_html_dir_index(.VQUEST_ARCHIVES_URL, style="IMGT", suffix=".zip",
                          connecttimeout=get_IMGT_connecttimeout())
}

### If 'as.df' is TRUE then the listing is returned as a data.frame
### with 3 columns (Name, Last modified, Size) and 1 row per .zip file.
list_archived_IMGT_zips <- function(as.df=FALSE, recache=FALSE)
{
    if (!isTRUEorFALSE(as.df))
        stop(wmsg("'as.df' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(recache))
        stop(wmsg("'recache' must be TRUE or FALSE"))
    listing <- .IMGT_cache[["ARCHIVES_TABLE"]]
    if (is.null(listing) || recache) {
        listing <- .fetch_list_of_archived_IMGT_zips()
        .IMGT_cache[["ARCHIVES_TABLE"]] <- listing
    }
    if (!as.df)
        listing <- listing[ , "Name"]
    listing
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### download_and_unzip_IMGT_release()
###

.download_and_unzip_latest_IMGT_zip <- function(exdir, ...)
{
    zip_filename <- paste0(VQUEST_REFERENCE_DIRECTORY, ".zip")

    ## Sometimes, after a new release, the IMGT people forget to make
    ## the zip file of the new release available. We're trying to detect
    ## this and fail graciously when it's the case.
    zip_url <- paste0(.VQUEST_DOWNLOAD_ROOT_URL, zip_filename)
    zip_exists <- urlExists(zip_url, connecttimeout=get_IMGT_connecttimeout())
    if (!zip_exists) {
        release <- get_latest_IMGT_release()
        stop(wmsg("It looks like the zip of the latest IMGT/V-QUEST ",
                  "release (", release, ") is not available (yet?) at ",
                  .VQUEST_DOWNLOAD_ROOT_URL),
             "\n  ",
             wmsg("Please install an older release in the meantime."))
    }

    local_zip <- download_as_tempfile(.VQUEST_DOWNLOAD_ROOT_URL, zip_filename,
                                      ...)
    nuke_file(exdir)
    unzip(local_zip, exdir=exdir)
}

.get_archived_IMGT_zip <- function(release)
{
    stopifnot(isSingleNonWhiteString(release))
    all_zips <- list_archived_IMGT_zips()
    idx <- grep(release, all_zips, fixed=TRUE)
    if (length(idx) == 0L)
        stop(wmsg("Anomaly: no .zip file found at ",
                  .VQUEST_ARCHIVES_URL, " for release ", release))
    if (length(idx) > 1L)
        stop(wmsg("Anomaly: more that one .zip file found at ",
                  .VQUEST_ARCHIVES_URL, " for release ", release))
    all_zips[[idx]]
}

.unzip_archived_IMGT_zip <- function(zipfile, release, exdir)
{
    nuke_file(exdir)
    unzip(zipfile, exdir=exdir, junkpaths=TRUE)
    zip_filename <- paste0(VQUEST_REFERENCE_DIRECTORY, ".zip")
    local_zip <- file.path(exdir, zip_filename)
    unzip(local_zip, exdir=exdir)
    unlink(local_zip)
}

.download_and_unzip_archived_IMGT_zip <- function(release, exdir, ...)
{
    archived_zip_filename <- .get_archived_IMGT_zip(release)
    archived_zipfile <- download_as_tempfile(.VQUEST_ARCHIVES_URL,
                                             archived_zip_filename, ...)
    .unzip_archived_IMGT_zip(archived_zipfile, release, exdir)
}

### Download and unzip in 'exdir'.
download_and_unzip_IMGT_release <- function(release, exdir, ...)
{
    if (dir.exists(exdir))
        nuke_file(exdir)
    if (release == get_latest_IMGT_release()) {
        .download_and_unzip_latest_IMGT_zip(exdir, ...)
    } else {
        .download_and_unzip_archived_IMGT_zip(release, exdir, ...)
    }
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### normalize_IMGT_organism()
###

### List manually scrapped from https://www.imgt.org/IMGTrepertoire/Proteins/
### on September 11, 2026.
### The names on the '.IMGT_ORGANISMS' vector are the common names, and the
### vector elements the corresponding latin names. Common names must be in
### lower case, latin names must be in lower case with their 1st letter
### capitalized. No spaces or dashes! The only allowed non-letter character
### is the underscore (_).
.IMGT_ORGANISMS <- c(

    ## --- Mammalia ---
    human="Homo_sapiens",
    house_mouse="Mus_musculus",
    arabian_camel="Camelus_dromedarius",
    llama="Lama_glama",
    pig="Sus_scrofa",
    sheep="Ovis_aries",
    norway_rat="Rattus_norvegicus",
    rabbit="Oryctolagus_cuniculus",
    rhesus_monkey="Macaca_mulatta",
    platypus="Ornithorhynchus_anatinus",
    alpaca="Vicugna_pacos",
    crab_eating_macaque="Macaca_fascicularis",
    bovine="Bos_taurus", cow="Bos_taurus",
    dog="Canis_lupus_familiaris",
    horse="Equus_caballus",
    western_lowland_gorilla="Gorilla_gorilla_gorilla",
    ring_tailed_lemur="Lemur_catta",
    sumatran_orangutan="Pongo_abelii",
    domestic_ferret="Mustela_putorius_furo",
    bornean_orangutan="Pongo_pygmaeus",
    american_mink="Neogale_vison",
    domestic_cat="Felis_catus",
    goat="Capra_hircus",
    common_bottlenose_dolphin="Tursiops_truncatus",
    chimpanzee="Pan_troglodytes",

    ## --- Sauropsida ---
    spectacled_caiman="Caiman_crocodilus",
    chicken="Gallus_gallus",

    ## --- Amphibia ---
    african_clawed_frog="Xenopus_laevis",

    ## --- Teleostei ---
    japanese_flounder="Paralichthys_olivaceus",
    channel_catfish="Ictalurus_punctatus",
    arctic_char="Salvelinus_alpinus",
    atlantic_cod="Gadus_morhua",
    emerald_rockcod="Trematomus_bernacchii",
    goldfish="Carassius_auratus",
    blackfin_icefish="Chaenocephalus_aceratus",
    ladyfish="Elops_saurus",
    atlantic_salmon="Salmo_salar",
    torafugu="Takifugu_rubripes",
    rainbow_trout="Oncorhynchus_mykiss",
    black_rockcod="Notothenia_coriiceps",
    zebrafish="Danio_rerio",
    spotted_wolffish="Anarhichas_minor",
    grass_carp="Ctenopharyngodon_idella",

    ## --- Chondrichthyes ---
    nurse_shark="Ginglymostoma_cirratum",
    bull_shark="Carcharhinus_leucas",
    clearnose_skate="Raja_eglanteria",
    horn_shark="Heterodontus_francisci",
    little_skate="Leucoraja_erinacea",
    spotted_ratfish="Hydrolagus_colliei",
    sandbar_shark="Carcharhinus_plumbeus",
    spotted_wobbegong="Orectolobus_maculatus",

    ## --- other ---
    marbled_lungfish="Protopterus_aethiopicus"
)

.lookup_imgt_organisms <- function(organism, imgt_organisms)
{
    stopifnot(isSingleNonWhiteString(organism), is.character(imgt_organisms))
    common_names <- names(imgt_organisms)
    stopifnot(!is.null(common_names))

    ## Find latin names or common names equal to 'organism'.
    idx <- which(imgt_organisms == organism | common_names == organism)
    if (length(idx) != 0L)
        return(idx)

    ## Find latin names or common names that have a part equal to 'organism'.
    Lparts <- CharacterList(strsplit(common_names, "_"))
    Rparts <- CharacterList(strsplit(imgt_organisms, "_"))
    idx <- which(any(Lparts == organism) | any(Rparts == organism))
    if (length(idx) != 0L)
        return(idx)

    ## Find latin names or common names that start with 'organism' or
    ## have a part that starts with 'organism'.
    Lsuffixes <- substr(Lparts, 1L, nchar(organism))
    Rsuffixes <- substr(Rparts, 1L, nchar(organism))
    which(startsWith(imgt_organisms, organism) |
          startsWith(common_names, organism) |
          any(Lsuffixes == organism) |
          any(Rsuffixes == organism))
}

.stop_on_ambiguous_imgt_organism <- function(organism, matched_imgt_organisms)
{
    in1string <- paste0(matched_imgt_organisms,
                        " (", names(matched_imgt_organisms), ")",
                        collapse=", ")
    stop(wmsg("Ambiguous organism abbreviation: ", organism),
         "\n  ",
         wmsg("Matches: ", in1string))
}

normalize_IMGT_organism <- function(organism)
{
    if (!isSingleNonWhiteString(organism))
        stop(wmsg("'organism' must be a single (non-empty) string"))
    organism <- chartr(" -", "__", trimws2(organism))
    organism <- tolower(gsub("_+", "_", organism))
    imgt_organisms <- tolower(.IMGT_ORGANISMS)
    idx <- .lookup_imgt_organisms(organism, imgt_organisms)
    idx <- idx[!duplicated(imgt_organisms[idx])]
    if (length(idx) >= 2L)
        .stop_on_ambiguous_imgt_organism(organism, .IMGT_ORGANISMS[idx])
    if (length(idx) == 1L)
        return(.IMGT_ORGANISMS[[idx]])
    paste0(toupper(substr(organism, 1L, 1L)),
           substr(organism, 2L, nchar(organism)))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### find_organism_in_IMGT_store()
### list_organisms_in_IMGT_store()
###

list_organisms_in_IMGT_store <- function(IMGT_store)
{
    stopifnot(isSingleNonWhiteString(IMGT_store), dir.exists(IMGT_store))
    refdir <- file.path(IMGT_store, VQUEST_REFERENCE_DIRECTORY)
    if (!dir.exists(refdir))
        stop(wmsg("Anomaly: directory ", refdir, " not found"))
    sort(list.files(refdir))
}

### 'IMGT_store' must be the path to the local store of a given IMGT release.
### Returns the path to the subdir of 'IMGT_store' that corresponds to the
### specified organism. For example, for IMGT release 202449-1 and Homo
### sapiens, this path is:
###     <igblastr-cache>
###     └── store
###         └── IMGT-releases
###             └── 202449-1
###                 └── IMGT_V-QUEST_reference_directory
###                     └──  Homo_sapiens
find_organism_in_IMGT_store <- function(organism, IMGT_store)
{
    stopifnot(isSingleNonWhiteString(organism))
    imgt_organisms <- list_organisms_in_IMGT_store(IMGT_store)
    idx <- match(tolower(organism), tolower(imgt_organisms))
    if (!is.na(idx)) {
        refdir <- file.path(IMGT_store, VQUEST_REFERENCE_DIRECTORY)
        return(file.path(refdir, imgt_organisms[[idx]]))
    }
    in1string <- paste0("\"", imgt_organisms, "\"", collapse=", ")
    stop(wmsg(organism, ": organism not found in ",
              "IMGT/V-QUEST release ", basename(IMGT_store), "."),
         "\n  ",
         wmsg("Available organisms: ", in1string, "."))
}

