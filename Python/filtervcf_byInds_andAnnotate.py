#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 16 21:58:28 2022

@author: nataliazajac
"""

import pandas as pd
import gffpandas.gffpandas as gffpd
import io
import sys

def filter_vcf(vcf,list_of_samples,gff):
    """Main function."""
    
    usage = """
    %python filtervcf.py [vcf] [list_of_samples] [gff]
    Note:
     list of samples has a single sample per line
    """

    def read_vcf(path):
        with open(path, 'r') as f:
            lines = [l for l in f if not l.startswith('##')]
        return pd.read_csv(
            io.StringIO(''.join(lines)),
            dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},
            sep='\t'
        ).rename(columns={'#CHROM': 'CHROM'})

    your_vcf=vcf
    list_of_samples=list_of_samples
    annotation=gff

    table=read_vcf(your_vcf)
    with open(list_of_samples) as f:
        samples = f.read().splitlines()
    annotation = gffpd.read_gff3(annotation)
    
    snps_per_Sample=dict()

    for i in samples:
        m = table[table[i] != "./.:.:.:.:.:.:."]["POS"].tolist()    
        snps_per_Sample[i]=m
    
    filtered_df = annotation.filter_feature_of_type(['gene'])

    gene_per_snp=dict()

    for sample,snps in snps_per_Sample.items():
        for i in snps:
            df = filtered_df.df[(filtered_df.df['start'] <= i) & (filtered_df.df['end'] >= i)]
            if not df.empty:
                gene_per_snp[i] = df['attributes'].str.split(';', expand=True)[1].tolist()
   
    snps_table = pd.melt(pd.DataFrame.from_dict(snps_per_Sample, orient='index').reset_index(), id_vars="index")[["index", "value"]]
    snps_table = snps_table.rename({"index" : "sample", "value" : "pos"}, axis = 1)
    genes_snps = pd.DataFrame.from_dict(gene_per_snp, orient='index').reset_index().rename({"index": "pos", 0: "Gene name"}, axis = 1)
    final = pd.merge(snps_table, genes_snps, on = "pos")

    final.to_csv("my_snps_of_interest.csv", sep = "\t")


if __name__=="__main__":
    filter_vcf(sys.argv[1],sys.argv[2],sys.argv[3])
