\name{compute_auxdata}

\alias{compute_auxdata}
\alias{validate_and_complete_auxdata_with_ref}
\alias{solve_cdr3_ends_using_fwr4_aa_comparisons}
\alias{solve_cdr3_ends_using_fwr4_dna_comparisons}
\alias{solve_cdr3_ends_using_fwr4_dna_PWM}

\title{Compute IgBLAST auxiliary data}

\description{
  Utility functions to generate \emph{auxiliary data} for a set of
  germline J gene allele sequences.

  See \code{?\link{load_auxdata}} for more information about IgBLAST
  \emph{auxiliary data}.
}

\usage{
## The primary utility for computing "auxiliary data":

compute_auxdata(J_alleles, codon_starts=NULL, no.warnings=FALSE)

## Additional utilities that can be used to complete the "auxiliary data"
## returned by compute_auxdata():

validate_and_complete_auxdata_with_ref(auxdata, ref_auxdata)

solve_cdr3_ends_using_fwr4_aa_comparisons(auxdata, J_alleles,
                                          max.dist=2, standout.by=5)

solve_cdr3_ends_using_fwr4_dna_comparisons(auxdata, J_alleles,
                                           max.dist=4, standout.by=15)

solve_cdr3_ends_using_fwr4_dna_PWM(auxdata, J_alleles,
                                   min.score=0.80, standout.by=0.40)
}

\arguments{
  \item{J_alleles}{
    A \link[Biostrings]{DNAStringSet} object containing germline J
    gene allele sequences.
  }
  \item{codon_starts}{
    When supplied, must be a named integer vector \emph{parallel} to
    \code{J_alleles}, that is, it must have the same length and names
    as \code{J_alleles}. For each allele in \code{J_alleles},
    \code{codon_starts} should indicate the start position (1-based) of
    the first codon in the allele. This should be a number >= 1 and <= 3.
    \code{NA}s are allowed when the coding frame is not known.

    Note that the "coding frame start" is 0-based so can simply be inferred
    from the "codon start" by subtracting 1.
  }
  \item{no.warnings}{
    Set to \code{TRUE} to suppress the warnings that \code{compute_auxdata()}
    otherwise issues when the returned data.frame has \code{NA}s in its
    \code{coding_frame_start} column or negative values in its \code{cdr3_end}
    column.
  }
  \item{auxdata}{
    A data.frame as returned by \code{\link{compute_auxdata}()} that
    contains annotations for the alleles in \code{J_alleles}.
  }
  \item{ref_auxdata}{
    A data.frame containing "reference auxiliary data". Typically the
    \emph{auxiliary data} provided by IgBLAST for a given organism.
  }
  \item{max.dist}{
    The maximum Hamming distance allowed between the content of the sliding
    window along the J allele to solve and the set of known FWR4 sequences.
    By default, \code{solve_cdr3_ends_using_fwr4_aa_comparisons()} uses
    \code{max.dist=2} and \code{solve_cdr3_ends_using_fwr4_dna_comparisons()}
    uses \code{max.dist=4}.
    See Details section below for the details.
  }
  \item{standout.by}{
    See Details section below for the details.
  }
  \item{min.score}{
    The minimum score allowed for the sliding Position Weight Matrix
    along the the J allele to solve. See Details section below for the
    details.
  }
}

\details{
  \code{compute_auxdata()} is meant to be used first, to generate
  \emph{auxiliary data} from scratch for a set of germline J gene
  allele sequences. Note that the data.frame returned by
  \code{compute_auxdata()} can contain \code{NA}s in its
  \code{cdr3_end} column. We call these \emph{unsolved CDR3 ends}.

  \code{validate_and_complete_auxdata_with_ref()} and the
  \code{solve_cdr3_ends_using_fwr4_*()} functions are meant to be
  used on the \emph{auxiliary data} returned by \code{compute_auxdata()}
  to complete it, that is, to eliminate \code{NA}s from it when possible.

  \subsection{The \code{compute_auxdata()} function}{
    \code{compute_auxdata()} uses a motif-based approach to identify
    the CDR3/FWR4 junction in the supplied J allele sequences.

    For human BCR germline J alleles, the FWR4 region is expected to start
    with the following amino acid motifs (X represents any amino acid):
    \itemize{
      \item "WGXG" on the heavy chain (IGH locus);
      \item "FGXG" on the light chain (IGK and IGL loci).
    }
    Things might deviate slightly for other organisms, or for TCR germline
    J alleles. For example, for TCR germline J alleles, the FWR4 region is
    expected to start with the "FGXG" motif, regardless of the locus.

    One complication is that the supplied J allele sequences are nucleotide
    sequences for which the coding frame might not be known.
    \code{compute_auxdata()} uses two distinct strategies to search for
    canonical motif "WGXG" or "FGXG" in a given J allele sequence, depending
    on whether the coding frame is known or not:
    \enumerate{
      \item If the coding frame is known (i.e. a vector of codon starts
            is passed thru the \code{codon_starts} argument and it reports
            a non-\code{NA} value for the allele), then this information is
            used to translate the coding region of the allele to an amino
            acid sequence. The function then tries to locate canonical
            motif "WGXG" or "FGXG" on the translated sequence.

            Note that if the allele is from the IGH locus, then canonical
            motif "WGXG" is used. For all other loci (IGK, IGL, TRA, TRB,
            TRG, TRD loci), canonical motif "FGXG" is used.

      \item If the coding frame is not known (i.e. the \code{codon_starts}
            argument is not set or it's set but the codon start reported
            for the allele is \code{NA}), then the search above is performed
            for all 3 possible coding frames, that is, for all 3 possible
            codon start values (1, 2, and 3).
    }
    Next, the CDR3/FWR4 junction is identified based on the result of the
    motif search (and the \code{cdr3_end} set accordingly):
    \itemize{
      \item If the motif search didn't return any match, then the CDR3/FWR4
            junction is considered to remain unknown. In this case, the
            \code{cdr3_end} is set to \code{NA}.

      \item If the motif search returned exactly one match, then the match
            is considered to mark the start of the FWR4. This means that
            the CDR3/FWR4 junction is now identified, and so the
            \code{cdr3_end} is set accordingly.

      \item If the motif search returned more than one match, then only the
            first match is considered to mark the start of the FWR4.
            All downstream matches are ignored.
    }

    Finally, the \code{coding_frame_start} is obtained as follow. For a given
    J allele:
    \enumerate{
      \item If a vector of codon starts was passed thru the \code{codon_starts}
            argument and it reports a non-\code{NA} value for the allele,
            then the \code{coding_frame_start} (0-based) is set to the
            reported codon start (1-based) minus 1.

      \item If the above cannot be performed, but if the CDR3/FWR4 junction
            of the allele was successfully identified by the motif search,
            then the \code{coding_frame_start} is inferred from the
            \code{cdr3_end}.

      \item If none of the two operations above can be performed, then the
            \code{coding_frame_start} is set to \code{NA} with a warning.
    }
  }

  \subsection{The \code{validate_and_complete_auxdata_with_ref()} function}{
    \code{validate_and_complete_auxdata_with_ref()} does two things:
    \enumerate{
      \item Validation: The function validates the \emph{auxiliary data}
            passed thru its first argument (\code{auxdata}) against the
            \emph{auxiliary data} passed thru its second argument
            (\code{ref_auxdata}). This is done by comparing the entries in
            the former with the entries in the latter for the J alleles that
            are annotated in both, and by raising an error if they disagree.
      \item Completion: The function replaces each \code{NA} in \code{auxdata}
            with the corresponding value in \code{ref_auxdata}, if the latter
            is available.
    }
  }

  \subsection{The \code{solve_cdr3_ends_using_fwr4_*_comparisons()} functions}{
    For each unsolved J allele (i.e. for each J allele with an unknown
    CDR3/FWR4 junction), \code{solve_cdr3_ends_using_fwr4_aa_comparisons()}
    tries to identify the CDR3/FWR4 junction as follow:
    \itemize{
      \item The function moves a sliding window of 10 amino acids along the
            allele sequence.
      \item For each position p of the window, it computes the Hamming distance
            D(p) between the sequence in the window and the set of known FWR4
            sequences.
      \item The "best window" is the window that minimizes D(p).
      \item The position Bp of the "best window" must satisfy the two
            following criteria in order to be considered the start of
            the FWR4 of the unsolved J allele:
            \enumerate{
              \item D(Bp) <= \code{max.dist}.
              \item D(p) >= \code{max.dist} + \code{standout.by} for any
                    p != Bp. In other words, all other windows must have
                    \code{max.dist} + \code{standout.by} or more mismatches
                    with any of the known FWR4 sequences.
            }
    }

    \code{solve_cdr3_ends_using_fwr4_dna_comparisons()} works the same way
    but at the nucleotide level i.e. it uses a sliding window of 30 nucleotides
    along the allele sequence. Only the positions of the window that are
    compatible with the coding frame are considered.
  }

  \subsection{The \code{solve_cdr3_ends_using_fwr4_dna_PWM()} function}{
    \code{solve_cdr3_ends_using_fwr4_dna_PWM()} uses the first 30 nucleotides
    of each known FWR4 to compute a 4 x 30 Position Weight Matrix (PWM) for
    the whole set of known FWR4 regions. Then, for each unsolved J allele (i.e.
    for each J allele with an unknown CDR3/FWR4 junction), it tries to
    identify the CDR3/FWR4 junction as follow:
    \itemize{
      \item The function moves a sliding window of 30 nucleotides along
            the allele sequence. Only the positions of the window that
            are compatible with the coding frame are considered.
      \item For each position p of the window, it computes the score S(p)
            of the sequence in the window against the PWM.
      \item The "best window" is the window that maximizes S(p).
      \item The position Bp of the "best window" must satisfy the two
            following criteria in order to be considered the start of
            the FWR4 of the unsolved J allele:
            \enumerate{
              \item S(Bp) >= \code{min.score}.
              \item S(p) <= \code{min.score} - \code{standout.by} for any
                    p != Bp.
            }
    }
    Note that the \code{\link[Biostrings]{PWM}()} function from the
    \pkg{Biostrings} package is used to compute the Position Weight Matrix.
    See \code{?Biostrings::\link[Biostrings]{PWM}} for more information.
  }
}

\value{
  All the functions documented in this man page return a data.frame
  with 1 row per supplied J allele and the same columns as the data.frame
  returned by \code{read_auxdata()}.
  See \code{?\link{read_auxdata}} for more information.
}

\seealso{
  \itemize{
    \item \code{\link{print_J_alleles}} to display a set of J allele
          sequences where the sequences are justified w.r.t. their
          CDR3/FWR4 junction.

    \item \code{\link{load_auxdata}} to access IgBLAST \emph{auxiliary data}.

    \item \code{\link{load_germline_sequences}} to load the nucleotide
          sequences of the gene alleles stored in a cached germline db.

    \item \code{\link{compute_V_gene_delineations}} to annotate a set of
          germline V gene allele sequences.

    \item \url{https://ncbi.github.io/igblast/cook/How-to-set-up.html}
          for important information about IgBLAST \emph{auxiliary data}.

    \item \link[Biostrings]{DNAStringSet} and \link[Biostrings]{AAStringSet}
          objects in the \pkg{Biostrings} package.

    \item The \code{\link{igblastn}} function to run the \code{igblastn}
          \emph{standalone executable} included in IgBLAST from R. This
          is the main function in the \pkg{igblastr} package.

    \item IgBLAST is described at
          \url{https://pubmed.ncbi.nlm.nih.gov/23671333/}.
  }
}

\examples{
## ---------------------------------------------------------------------
## 1. BASIC compute_auxdata() EXAMPLE
## ---------------------------------------------------------------------

## Let's load the J allele sequences from _OGRDB.human.IGH+IGK+IGL.202605:
db_name <- "_OGRDB.human.IGH+IGK+IGL.202605"
J_alleles <- load_germline_sequences(db_name, region_types="J")
J_alleles  # DNAStringSet object

computed_auxdata <- compute_auxdata(J_alleles)
head(computed_auxdata)

## The CDR3/FWR4 junction was found for all the alleles:
anyNA(computed_auxdata$cdr3_end)  # FALSE

## This is also reflected by print_J_alleles() which is able to display
## the CDR3/FWR4 junction for all alleles:
print_J_alleles(J_alleles, computed_auxdata, translate=TRUE)

## Note that 'computed_auxdata' is identical to the "auxiliary data"
## included in the germline db (which is provided by AIRR/OGRDB):
stopifnot(identical(computed_auxdata, load_auxdata(db_name)))

## Also note that 'computed_auxdata' is in agreement with the "auxiliary
## data" provided by IgBLAST for human, except for alleles IGHJ6*02 and
## IGHJ6*03 (the 'extra_bps' column contains incorrect values for these
## two alleles, 0's instead of 1's):
human_auxdata0 <- load_auxdata("human", which="original")
bad_alleles <- c("IGHJ6*02", "IGHJ6*03")
subset(human_auxdata0, allele_name \%in\% bad_alleles)

## We manually correct this:
fixme <- human_auxdata0[ , "allele_name"] \%in\% bad_alleles
human_auxdata0[fixme, "extra_bps"] <- 1L  # replace 0L with 1L

## Now the data in 'computed_auxdata' matches exactly the corresponding
## data in 'human_auxdata0':
m <- match(computed_auxdata[ , "allele_name"],
           human_auxdata0[ , "allele_name"])
human_auxdata <- human_auxdata0[m, ]
rownames(human_auxdata) <- NULL
stopifnot(identical(computed_auxdata, human_auxdata))

## ---------------------------------------------------------------------
## 2. ADVANCED compute_auxdata() EXAMPLE
## ---------------------------------------------------------------------

## In this example we're going to annotate some mouse J alleles that we
## retrieve directly from IMGT/V-QUEST.

## We start by downloading the IMGT FASTA files for mouse to a
## temporary directory:
dir.create(mouse_destdir <- tempfile("fasta_dir_"))
download_IMGT_germline_sequences("202631-1",
                                 organism="Mus_musculus",
                                 destdir=mouse_destdir)
## List the downloaded files:
list.files(mouse_destdir)  # 7 FASTA files

## Let's load the J allele sequences for the lambda locus (IGLJ.fasta
## file):
J_alleles <- readDNAStringSet(file.path(mouse_destdir, "IGLJ.fasta"))
J_alleles  # only 8 sequences

## First attempt at annotating the J alleles with compute_auxdata():
auxdata1 <- compute_auxdata(J_alleles)  # emits warning!

## Two alleles could not be annotated:
auxdata1

## We can improve the situation by supplying the codon starts to
## compute_auxdata(). These are provided by IMGT in the IMGT FASTA
## headers:
parsed_headers <- parse_imgt_fasta_headers(names(J_alleles))
parsed_headers[1:3, ]
codon_starts <- as.integer(parsed_headers[ , "codon_start"])

## Set the allele names on the 'codon_starts' vector:
names(codon_starts) <- parsed_headers[ , "allele_name"]
codon_starts

## Second attempt at annotating the J alleles with compute_auxdata():
auxdata2 <- compute_auxdata(J_alleles, codon_starts=codon_starts)
auxdata2

## Supplying the codon starts had the effect of eliminating the NAs
## that we saw in the "coding_frame_start" and "extra_bps" columns
## of 'auxdata1'.

## However, alleles IGLJ2P*01 and IGLJ3P*01 are still unsolved.
## This is reflected by print_J_alleles() which cannot display
## the CDR3/FWR4 junction for these two alleles:
print_J_alleles(J_alleles, auxdata2, translate=TRUE)

## In the next examples below, we're going to use
## validate_and_complete_auxdata_with_ref() and the
## solve_cdr3_ends_using_fwr4_*_comparisons() functions to try
## and solve these alleles.

## ---------------------------------------------------------------------
## 3. validate_and_complete_auxdata_with_ref()
## ---------------------------------------------------------------------

igblast_auxdata <- load_auxdata("mouse", which="original")

## Note that IGLJ2P*01 is not annotated in 'igblast_auxdata',
## but 'IGLJ3P*01' is:
"IGLJ2P*01" \%in\% igblast_auxdata[ , "allele_name"]  # FALSE
"IGLJ3P*01" \%in\% igblast_auxdata[ , "allele_name"]  # TRUE

## validate_and_complete_auxdata_with_ref() replaces each NA in 'auxdata2'
## with the corresponding value in 'igblast_auxdata', if the latter is
## available:
auxdata3 <-
    validate_and_complete_auxdata_with_ref(auxdata2, igblast_auxdata)
auxdata3  # IGLJ3P*01 is solved!

## Note that the function also compares the data in 'auxdata2' with the
## data in 'igblast_auxdata' (for the J alleles that are annotated in
## both), and raises an error if they disagree.

## ---------------------------------------------------------------------
## 4. solve_cdr3_ends_using_fwr4_aa_comparisons() and
##    solve_cdr3_ends_using_fwr4_dna_comparisons()
## ---------------------------------------------------------------------

## Finally IGLJ2P*01 can be solved with:
auxdata4a <-
    solve_cdr3_ends_using_fwr4_aa_comparisons(auxdata3, J_alleles)
auxdata4a
stopifnot(!anyNA(auxdata4a$cdr3_end))

## or with:
auxdata4b <-
    solve_cdr3_ends_using_fwr4_dna_comparisons(auxdata3, J_alleles)
stopifnot(identical(auxdata4a, auxdata4b))

print_J_alleles(J_alleles, auxdata4a, translate=TRUE,
                igblast_organism="mouse")

## ---------------------------------------------------------------------
## 5. solve_cdr3_ends_using_fwr4_dna_PWM()
## ---------------------------------------------------------------------

## In this example we're going to annotate some ferret J alleles that we
## retrieve directly from IMGT/V-QUEST.

dir.create(ferret_destdir <- tempfile("fasta_dir_"))
download_IMGT_germline_sequences("202631-1",
                                 organism="Mustela_putorius_furo",
                                 destdir=ferret_destdir)
## List the downloaded files:
list.files(ferret_destdir)  # 7 FASTA files

## Let's load the J allele sequences for the kappa locus (IGKJ.fasta
## file):
J_alleles <- readDNAStringSet(file.path(ferret_destdir, "IGKJ.fasta"))
J_alleles  # only 5 sequences

## Extract their codon starts:
parsed_headers <- parse_imgt_fasta_headers(names(J_alleles))
codon_starts <- as.integer(parsed_headers[ , "codon_start"])
names(codon_starts) <- parsed_headers[ , "allele_name"]

## First attempt at annotating the J alleles with compute_auxdata():
auxdata1 <- compute_auxdata(J_alleles, codon_starts=codon_starts)
auxdata1

## Allele IGKJ2*01 didn't get solved:
print_J_alleles(J_alleles, auxdata1, translate=TRUE)

## Sill doesn't solve IGKJ2*01:
solve_cdr3_ends_using_fwr4_aa_comparisons(auxdata1, J_alleles)

## However, the following calls are able to solve IGKJ2*01:
auxdata2a <-
    solve_cdr3_ends_using_fwr4_dna_comparisons(auxdata1, J_alleles)
auxdata2a
stopifnot(!anyNA(auxdata2a$cdr3_end))

auxdata2b <- solve_cdr3_ends_using_fwr4_dna_PWM(auxdata1, J_alleles)
auxdata2b
stopifnot(identical(auxdata2a, auxdata2b))

print_J_alleles(J_alleles, auxdata2a, translate=TRUE)
}

\keyword{utilities}
