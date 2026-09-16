
test_that("normalize_IMGT_organism()", {
    normalize_IMGT_organism <- igblastr:::normalize_IMGT_organism

    expect_identical(normalize_IMGT_organism("human"), "Homo_sapiens")
    expect_identical(normalize_IMGT_organism("hom"), "Homo_sapiens")
    expect_identical(normalize_IMGT_organism("sapiens"), "Homo_sapiens")
    expect_identical(normalize_IMGT_organism("sap"), "Homo_sapiens")

    expect_identical(normalize_IMGT_organism("mouse"), "Mus_musculus")
    expect_identical(normalize_IMGT_organism("mus"), "Mus_musculus")

    expect_identical(normalize_IMGT_organism("rat"), "Rattus_norvegicus")
    expect_identical(normalize_IMGT_organism("nor"), "Rattus_norvegicus")

    expect_identical(normalize_IMGT_organism("mul"), "Macaca_mulatta")

    expect_identical(normalize_IMGT_organism("rab"), "Oryctolagus_cuniculus")
    expect_identical(normalize_IMGT_organism("ory"), "Oryctolagus_cuniculus")
    expect_identical(normalize_IMGT_organism("cun"), "Oryctolagus_cuniculus")

    expect_identical(normalize_IMGT_organism("camel"), "Camelus_dromedarius")
    expect_identical(normalize_IMGT_organism("cat"), "Felis_catus")
    expect_identical(normalize_IMGT_organism("croc"), "Caiman_crocodilus")
    expect_identical(normalize_IMGT_organism("dolphin"), "Tursiops_truncatus")
    expect_identical(normalize_IMGT_organism("fer"), "Mustela_putorius_furo")
    expect_identical(normalize_IMGT_organism("fur"), "Mustela_putorius_furo")
    expect_identical(normalize_IMGT_organism("pig"), "Sus_scrofa")

    expect_identical(normalize_IMGT_organism("crab"), "Macaca_fascicularis")
    expect_identical(normalize_IMGT_organism("eating"), "Macaca_fascicularis")
    expect_identical(normalize_IMGT_organism("macaque"), "Macaca_fascicularis")

    expect_identical(normalize_IMGT_organism("chimp"), "Pan_troglodytes")
    expect_identical(normalize_IMGT_organism("pan"), "Pan_troglodytes")
    expect_identical(normalize_IMGT_organism("trog"), "Pan_troglodytes")

    errmsg <- paste0("Pan_troglodytes \\(chimpanzee\\), ",
                     "Gallus_gallus \\(chicken\\)")
    expect_error2(normalize_IMGT_organism("chi"), errmsg)
    errmsg <- paste0("Pan_troglodytes \\(chimpanzee\\), ",
                     "Oncorhynchus_mykiss \\(rainbow_trout\\)")
    expect_error2(normalize_IMGT_organism("tro"), errmsg)

    expect_identical(normalize_IMGT_organism("abelii"), "Pongo_abelii")
    expect_identical(normalize_IMGT_organism("abe"), "Pongo_abelii")
    expect_identical(normalize_IMGT_organism("pygmaeus"), "Pongo_pygmaeus")
    expect_identical(normalize_IMGT_organism("pyg"), "Pongo_pygmaeus")

    errmsg <- paste0("Pongo_abelii \\(sumatran_orangutan\\), ",
                     "Pongo_pygmaeus \\(bornean_orangutan\\)")
    expect_error2(normalize_IMGT_organism("pongo"), errmsg)
    expect_error2(normalize_IMGT_organism("pon"), errmsg)
    expect_error2(normalize_IMGT_organism("orangutan"), errmsg)
    expect_error2(normalize_IMGT_organism("oran"), errmsg)

    expect_identical(normalize_IMGT_organism("goa"), "Capra_hircus")
    expect_identical(normalize_IMGT_organism("gol"), "Carassius_auratus")
    expect_identical(normalize_IMGT_organism("gor"), "Gorilla_gorilla_gorilla")
    expect_identical(normalize_IMGT_organism("gal"), "Gallus_gallus")

    expect_identical(normalize_IMGT_organism("fra"), "Heterodontus_francisci")
    expect_identical(normalize_IMGT_organism("fro"), "Xenopus_laevis")
    expect_identical(normalize_IMGT_organism("af"), "Xenopus_laevis")
    expect_identical(normalize_IMGT_organism("cla"), "Xenopus_laevis")

    expect_identical(normalize_IMGT_organism("rat"), "Rattus_norvegicus")
    expect_identical(normalize_IMGT_organism("ratf"), "Hydrolagus_colliei")

    expect_identical(normalize_IMGT_organism("gi"), "Ginglymostoma_cirratum")

    expect_identical(normalize_IMGT_organism("cab"), "Equus_caballus")

    expect_identical(normalize_IMGT_organism("horn"), "Heterodontus_francisci")
    expect_identical(normalize_IMGT_organism("hors"), "Equus_caballus")

    errmsg <- paste0("Vicugna_pacos \\(alpaca\\), ",
                     "Salvelinus_alpinus \\(arctic_char\\)")
    expect_error2(normalize_IMGT_organism("al"), errmsg)

    ## Case and outer whitespace are ignored, inner non-letter blocks (i.e.
    ## inner blocks made of whitespace and underscores) are treated as a
    ## single underscore:
    expect_identical(normalize_IMGT_organism("hOmO  SaPiens "), "Homo_sapiens")
    expect_identical(normalize_IMGT_organism(" RAt"), "Rattus_norvegicus")
    expect_identical(normalize_IMGT_organism("  NoRW "), "Rattus_norvegicus")
    expect_identical(normalize_IMGT_organism("   bOv "), "Bos_taurus")
    expect_identical(normalize_IMGT_organism(" CoW   "), "Bos_taurus")
    expect_identical(normalize_IMGT_organism("  TaUr "), "Bos_taurus")
    expect_identical(normalize_IMGT_organism("   boS"), "Bos_taurus")
    expect_identical(normalize_IMGT_organism("  bOs  tA "), "Bos_taurus")
    expect_identical(normalize_IMGT_organism("  bOs__tA "), "Bos_taurus")
    expect_identical(normalize_IMGT_organism("  bOs _ _tA "), "Bos_taurus")

    ## No error on unknown organism:
    expect_identical(normalize_IMGT_organism("  bIG _ _fOOt "), "Big_foot")
})

