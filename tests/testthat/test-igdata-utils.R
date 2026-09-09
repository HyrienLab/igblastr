read_igdata <- igblastr:::read_igdata
write_igdata <- igblastr:::write_igdata

.create_messy_igdata_file <- function()
{
    text <- c(
        "# comment line 1",
        "# comment line 2",
        "A   -1 a  \t  \t104",
        "B\t   2",
        "A   3 N/A 304",
        "C ",
        "D \t \t\t5  X\t\t504  ",
        "E 6 Z"
    )
    tmp_file <- tempfile()
    writeLines(text, tmp_file)
    tmp_file
}

.expected_df_from_messy_igdata_file <- function()
    data.frame(
        allele_name=c("A", "B", "A", "C", "D", "E"),
        offset=c(NA, 2L, 3L, NA, 5L, 6L),
        alias=c("a", NA, NA, NA, "X", "Z"),
        length=c(104L, NA, 304L, NA, 504L, NA)
    )

test_that("read_igdata()", {
    messy_file <- .create_messy_igdata_file()
    expected_df <- .expected_df_from_messy_igdata_file()
    col2class <- vapply(expected_df, class, character(1))
    df <- read_igdata(messy_file, col2class)
    expect_identical(df, expected_df)
})

test_that("write_igdata()", {
    df0 <- .expected_df_from_messy_igdata_file()
    tmp_file <- tempfile()
    write_igdata(df0, tmp_file)
    text <- readLines(tmp_file)
    expect_equal(length(text), 7)
    expect_equal(text[[1L]], "#allele_name, offset, alias, length")
    expect_equal(text[[2L]], "A\t-1\ta\t104")
    expect_equal(text[[3L]], "B\t2\tN/A\t-1")
    expect_equal(text[[4L]], "A\t3\tN/A\t304")
    expect_equal(text[[5L]], "C\t-1\tN/A\t-1")
    expect_equal(text[[6L]], "D\t5\tX\t504")
    expect_equal(text[[7L]], "E\t6\tZ\t-1")

    col2class <- vapply(df0, class, character(1))
    expect_identical(read_igdata(tmp_file, col2class), df0)

    df <- df0
    df$offset[[2L]] <- 0.5
    expect_error2(write_igdata(df, tmp_file),
                  "can only have character or integer columns")
    df <- df0
    for (s in c(NA_character_, "N/A", "", "x ", "x y")) {
        df$allele_name[[3L]] <- s
        expect_error2(write_igdata(df, tmp_file), "not allowed")
    }
    df <- df0
    df$offset[[3L]] <- -2L
    errmsg <- "negative values"
    expect_error2(write_igdata(df, tmp_file), errmsg)
    df <- df0
    df$alias[[4L]] <- ""
    errmsg <- "the empty string or strings with whitespaces"
    expect_error2(write_igdata(df, tmp_file), errmsg)
})

