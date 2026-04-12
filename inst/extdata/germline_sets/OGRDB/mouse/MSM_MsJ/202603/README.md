This is the `extdata/germline_sets/OGRDB/mouse/MSM_MsJ/202603/` folder
in the **igblastr** package.

This folder contains germline sequences (in FASTA format) obtained from
AIRR-community/OGRDB for mouse strain MSM/MsJ. The sequences are coming
from the following germline sets:
- `MSM/MsJ IGH` version 1
- `MSM/MsJ IGKV` version 1
- `MSM/MsJ IGLV` version 2
- `IGKJ (all strains)` version 1
- `IGLJ (all strains)` version 1

Note that these FASTA files were obtained programmatically by running the
following code in the folder on April 11, 2026:
```r
library(igblastr)
germline_sets <- c(`MSM/MsJ IGH`=1, `MSM/MsJ IGKV`=1, `MSM/MsJ IGLV`=2,
                   `IGKJ (all strains)`=1, `IGLJ (all strains)`=1)
download_OGRDB_germline_sequences("Mus musculus", germline_sets)
```

