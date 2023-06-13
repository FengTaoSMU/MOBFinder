#!/usr/bin/bash
# makeblastdb: 2.12.0+

makeblastdb -in ./02.relaxase_blastp/MOBB_homologous.fasta -dbtype prot -out ./03.relaxase_db/MOBB/MOBB_DB
makeblastdb -in ./02.relaxase_blastp/MOBC_homologous.fasta -dbtype prot -out ./03.relaxase_db/MOBC/MOBC_DB
makeblastdb -in ./02.relaxase_blastp/MOBF_homologous.fasta -dbtype prot -out ./03.relaxase_db/MOBF/MOBF_DB
makeblastdb -in ./02.relaxase_blastp/MOBH_homologous.fasta -dbtype prot -out ./03.relaxase_db/MOBH/MOBH_DB
makeblastdb -in ./02.relaxase_blastp/MOBL_homologous.fasta -dbtype prot -out ./03.relaxase_db/MOBL/MOBL_DB
makeblastdb -in ./02.relaxase_blastp/MOBM_homologous.fasta -dbtype prot -out ./03.relaxase_db/MOBM/MOBM_DB
makeblastdb -in ./02.relaxase_blastp/MOBP_homologous.fasta -dbtype prot -out ./03.relaxase_db/MOBP/MOBP_DB
makeblastdb -in ./02.relaxase_blastp/MOBQ_homologous.fasta -dbtype prot -out ./03.relaxase_db/MOBQ/MOBQ_DB
makeblastdb -in ./02.relaxase_blastp/MOBT_homologous.fasta -dbtype prot -out ./03.relaxase_db/MOBT/MOBT_DB
makeblastdb -in ./02.relaxase_blastp/MOBV_homologous.fasta -dbtype prot -out ./03.relaxase_db/MOBV/MOBV_DB
