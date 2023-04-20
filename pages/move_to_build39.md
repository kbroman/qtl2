---
layout: page
title: Resources to move to mouse genome build 39
description: "Resources for moving to mouse genome build GRCm39 (aka mm11), with a particular focus on mouse genetics projects related to the Collaborative Cross and Diversity Outbred"
---


The mouse genome build [GRCm39 was released on
2020-06-24](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001635.27/).
Here we describe available resources for using build 39 with R/qtl2,
and how to move your projects from build 38 to build 39. We focus
particularly on projects related to the Collaborative Cross or
Diversity Outbred (DO) mice, and on the use of the [MUGA SNP
arrays](https://doi.org/10.1534/g3.115.022087), such as the MegaMUGA
and
[GigaMUGA](https://www.neogen.com/categories/genotyping-arrays/gigamuga/)
arrays.

### Build 39 resources

We provide a variety of resources to support the analysis of mouse
genetic data in R/qtl2, particularly with the MUGA SNP arrays. These
have build updated to use mouse genome build GRCm39.

- The mouse genetic map of [Cox et al
  (2009)](https://doi.org/10.1534/genetics.109.105486) has been
  revised. See [the github repository
  CoxMapV3](https://github.com/kbroman/CoxMapV3).

- The MUGA array annotations have been revised to use positions from
  the GRCm39 genome build. See the [github repository
  MUGAarrays](https://github.com/kbroman/MUGAarrays).

- The data files with founder genotypes, organized for use with data
  in R/qtl2 format, have been revised for build GRCm39.

  - [`MM_processed_files_build39.zip`](https://doi.org/10.6084/m9.figshare.22666336.v1)
  - [`GM_processed_files_build39.zip`](https://doi.org/10.6084/m9.figshare.22666504.v1)
  - [`MMnGM_processed_files_build39.zip`](https://doi.org/10.6084/m9.figshare.22666510.v1)

- The [mmconvert package](https://github.com/rqtl/mmconvert), which is
  [available on CRAN](https://cran.r-project.org/package=mmconvert),
  can be used to convert map positions between build GRCm39 and the
  revised Cox genetic map and also includes a function
  `cross2_to_grcm39()` for converting a `cross2` object (created with
  `read_cross2()`, and for a cross using the MegaMUGA and/or GigaMUGA
  arrays) to build GRCm39.

- A new SQLite database with CC/DO founder variants from Sanger along
  with ensembl genes is [available on
  figshare](https://doi.org/10.6084/m9.figshare.22630642.v1).

  Download this ~10.2 GB database as `fv.2021.snps.db3`:

  ```r
  download.file("https://figshare.com/ndownloader/files/40157572", "fv.2021.snps.db3")
  ```

  You can then use `create_variant_query_func()` just as before:

  ```r
  qvf <- create_variant_query_func("fv.2021.snps.db3")
  ```

  The genes database uses different columns than the default in
  `create_gene_query_func()`, and so use the following:

  ```r
  qgf <- create_gene_query_func("fv.2021.snps.db3", chr="chromosome", start="start_position", stop="end_position")
  ```