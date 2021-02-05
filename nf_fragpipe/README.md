# FragPipe Nextflow Workflow

This nextflow workflow tries to replicate some aspects
of FragPipe.

The advantage of this workflow compared to FragPipe is that it

  A) Runs in a highly parallel fashion
  B) Can resume failed tasks

## Planned tasks

  * Add decoys and contaminants to FASTA file
  * Run search using MSFragger (open or closed)
  * FDR Filtering using philosopher
  * PTM profile using PTMShepherd (in open mod search)
  * Spectral library generation using EasyPQP

## Current issues

MSFragger cannot be shipped together with this pipeline due to
license restrictions. Currently, the path is hard-coded
in `runSearch`.

This also applies to the Thermo library that' shipped with MSFragger.
This is currently set in `params.thermo_ext_dir`.

The **PTMShepherd.jar** is also referenced externally, currently in
`params.ptm_shepherd`. I'm unsure whether we'd be allowed to ship it.

The `grosenberger/easypqp` container doesn't use any versions at the moment.

Some containers cause conflicts with the Unix user ids. (Failed write operations).
I solved this (in some local versions) using the `containerOptions = "--user 1000"`
directive (as part of the process). Probably, it's better to move this to a global
config of the whole workflow.

The pipeline is theoretically designed to support closed and
open modification search. 