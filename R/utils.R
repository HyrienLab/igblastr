
isSingleNonWhiteString <- function(x) isSingleString(x) && !grepl("^\\s*$", x)

urlExists <- function(url)
{
    response <- try(HEAD(url, user_agent("igblastr")), silent=TRUE)
    if (inherits(response, "try-error"))
        stop(as.character(response), "  Please check your internet connection.")
    response$status_code != 404L
}

getUrlContent <- function(url)
{
    response <- try(GET(url, user_agent("igblastr")), silent=TRUE)
    if (inherits(response, "try-error"))
        stop(as.character(response), "  Please check your internet connection.")
    if (response$status_code == 404L)
        stop(wmsg("Not Found (HTTP 404): ", url))
    stop_for_status(response)
    content(response)
}

system_command_works <- function(command, args=character())
{
    out <- try(suppressWarnings(system2(command, args=args,
                                        stdout=TRUE, stderr=TRUE)),
               silent=TRUE)
    if (inherits(out, "try-error"))
        return(FALSE)
    status <- attr(out, "status")
    is.null(status) || isTRUE(all.equal(status, 0L))
}

