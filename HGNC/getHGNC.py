#!/usr/bin/env python

import io
import os
import click
import re
import pandas as pd


class HGNCextractor:
  def __init__(self,path):
    with open(path, 'r') as f:
      lines = [l for l in f if not l.startswith('##')]
    self.df = pd.read_csv(
          io.StringIO(''.join(lines)),
          dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
                 'QUAL': str, 'FILTER': str, 'INFO': str},
          sep='\t'
      ).rename(columns={'#CHROM': 'CHROM'})

  def printHGNC(self):
    for idx,row in self.df.iterrows():
      info = row['INFO']
      anns = info.split("CSQ") # we are splitting snpEff and VEP
      vep_items = anns[1].split("|")
      print(row['CHROM'],row['POS'], vep_items[3],vep_items[23])

# This is the surrogate for main(): everything happens here

@click.command(context_settings = dict( help_option_names = ['-h', '--help'] ))
@click.option('--vcf',      '-v', type=str, help='extract HGNC info from this VCF file', required=True)

def getHGNC(vcf):
    extractor = HGNCextractor(vcf)
    extractor.printHGNC()

if __name__ == "__main__":
  getHGNC()
