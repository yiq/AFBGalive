AFBGAlive
=========

The streaming service that annotates vcf records with background allele frequencies found in 1000G and ExAC projects

Data file
=========
This utility requires a database file pre-filled with allele frequencies in 1000G and ExAC to work properly. The file need to be at afdata/af.db relative to the current directory. The utility 'vcfload' can be used to generate this file

Usage
=====

afbgalive
---------
```
afbgalive [-z] 
```

vcf records are read directly from stdin. If '-z' is specified, the output is then in bgzip compressed from

vcfload
-------
This utility is used to read allele frequency information from 1000G or ExAC vcf files, and populate a sqlite database
```
vcfload <db_filename> <tablename> < file
```
`db_filename` is the database file that will be eventually provided to afbgalive.
`tablename` can either be 'g1k' or 'exac', depending on which vcf file is provided as input
`file` is the published vcf file from 1000G or ExAC project

After importing all vcf files, a composite (chrom,pos) index needs to be created on both `g1k` and `exac` table manually to achieve reasonable performance.

