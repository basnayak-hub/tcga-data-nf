#!/usr/bin/env python3

import click
import pandas as pd
import yaml

@click.group()
def cli():
    pass

@click.command()
@click.argument('filename')
def linear(filename):
    print("fitting a linear model from: ", filename)
    
@click.command()
@click.option('-e', '--expression_file', 'expression_file', type=str, required=True,
              help='Path to file that has the sample names')
@click.option('-m', '--methylation_file', 'clinical_file', type=str, required=True,
              help='Path to file that has the methylayion data')
@click.option('-o', '--output_yaml', 'output_yaml', type=str, required=True,
              help='Path to file that has the clinical variables')
@click.option('--sample_col', type=str, show_default=True, default='bcr_patient_barcode',
              help='Column that stores the sample names')
@click.option('--group_col', type=str, show_default=True, default='gender',
              help='Column that stores the group categories')
@click.option('--is_patient_clinical', is_flag=True, show_default=True,
              help='Whether the clinical file is the patient_clinical csv')  
@click.option('--separate_tissues', is_flag=True, show_default=True,
              help='Group by tumor vs healthy tissues')  
def get_group_yaml(expression_file, clinical_file, output_yaml, sample_col='bcr_patient_barcode', group_col='gender', is_patient_clinical=False, separate_tissues=False):
    """Get data from clinical file and generate a yaml file of groups.
    This is useful for aggropanda and possibly ligress/bonobo"""
    
    print("Getting samples from: ", expression_file)
    with open(expression_file, 'r') as f:
        tab = pd.read_csv(f, sep='\t', index_col = 0, nrows = 10)

    print('There are %d samples in the expression data' %(len(tab.columns)))
    
    sample_df = pd.DataFrame()
    sample_df['sample'] =  tab.columns
    sample_df['patient'] = [i[:12] for i in  tab.columns]
    sample_df['tissue'] = [tissue_group_from_barcode(b) for b in sample_df['sample'].tolist()]
    
    
    print("Getting groups from: %s, variable: %s" %(clinical_file, group_col))
    with open(clinical_file,'r') as f:
        clinical = pd.read_csv(f, index_col = 0)
    
    if is_patient_clinical:
        print('Removing first two rows of the clinical table, they contain extra info that we do not need')
        clinical = clinical.iloc[2:,:]
        
    group_df = clinical.loc[:,[sample_col, group_col]]
    
    
    group_df = group_df.fillna(value = {group_col:'other'})

    group_df = sample_df.merge(group_df, left_on = 'patient', right_on = sample_col, how = 'left')

    diz_group={} 
    
    if separate_tissues:
        diz_group = {('-').join(g): t['sample'].tolist() for g,t in group_df.groupby(by = [group_col, 'tissue']) }
    else:
        diz_group = {g: t['sample'].tolist() for g,t in group_df.groupby(by = group_col) }

    for k,t in diz_group.items():
        print('In group %s there are %d elements ' %(str(k), len(t)))

    with open(output_yaml, 'w') as f:
        documents = yaml.dump(diz_group, f)

cli.add_command(linear)
cli.add_command(get_group_yaml)

if __name__ == '__main__':
    cli()