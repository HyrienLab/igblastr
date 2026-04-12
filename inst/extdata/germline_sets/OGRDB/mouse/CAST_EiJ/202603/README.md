This is the `extdata/germline_sets/OGRDB/mouse/CAST_EiJ/202603/` folder
in the **igblastr** package.

This folder contains germline sequences (in FASTA format) obtained from
AIRR-community/OGRDB for mouse strain CAST/EiJ. The sequences are coming
from the following germline sets:
- `CAST/EiJ IGH` version 2
- `CAST/EiJ IGKV` version 1
- `CAST/EiJ IGLV` version 2
- `IGKJ (all strains)` version 1
- `IGLJ (all strains)` version 1

Note that these FASTA files were obtained programmatically by running the
following code in the folder on April 11, 2026:
```r
library(igblastr)
germline_sets <- c(`CAST/EiJ IGH`=2, `CAST/EiJ IGKV`=1, `CAST/EiJ IGLV`=2,
                   `IGKJ (all strains)`=1, `IGLJ (all strains)`=1)
download_OGRDB_germline_sequences("Mus musculus", germline_sets)
```

