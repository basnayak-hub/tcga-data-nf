#!/usr/bin/env python

from netZooPy.dragon.dragon import *
import numpy as np
import os
import sys
import click
import pandas as pd
import yaml
from netZooPy.lioness.lioness_for_dragon import LionessDragon
from netZooPy.dragon.dragon import remove_zero_variance_preds, get_shrunken_covariance_dragon

@click.group()
def cli():
    pass

@click.command()
@click.argument('filename')
def linear(filename):
    print("fitting a linear model from: ", filename)


@click.command()
@click.option('-m', '--meth_file', 'meth_file', type=str, required=True,
              help='Path to the methylation file')
@click.option('-e', '--expr_file', 'expr_file', type=str, required=True,
              help='Path to expression file')
@click.option('-i', '--input_dragon', 'input_dragon', type=str, required=True,
              help='TSV, The input data for dragon is written here')
@click.option('-o', '--output_dragon', 'output_dragon', type=str, required=True,
              help='TSV, The output results for dragon are written here')
@click.option('--methylation_barcode', type=str, show_default=True, default='TCGAbarcode',
              help='Column that stores the sample names in methylation file')
@click.option('--expression_barcode', type=str, show_default=True, default='TCGAbarcode',
              help='Column that stores the sample names in expression file')
@click.option('--use_full_barcode', is_flag=True, show_default=True,
              help='If False, only the first 16 characters of the barcode are used')
def dragon(meth_file, expr_file, input_dragon, output_dragon, methylation_barcode='TCGAbarcode', expression_barcode='TCGAbarcode', use_full_barcode=False):
    """Run dragon"""

    #meth_file = sys.argv[1]
    #expr_file = sys.argv[2]
    #out_dir = sys.argv[3]
    #uuid = sys.argv[4]

    #read methylation
    meth = pd.read_table(meth_file,header=0,index_col=0,sep=",")
    print("Shape of the methylation table:")
    print(meth.shape)
    # read expression
    expr = pd.read_table(expr_file,header=0,index_col=0,sep=",")
    print("Shape of the expression table:")
    print(expr.shape)
    if use_full_barcode:
        pass
    else:
        expr[expression_barcode] = [i[:16] for i in expr[expression_barcode]]
        meth[methylation_barcode] = [i[:16] for i in meth[methylation_barcode]]
    
    # Add suffixes to the column names and merge tables
    meth = meth.add_suffix('_methylation')
    meth = meth.rename(index=str, columns={'TCGAbarcode_methylation':'TCGAbarcode'})

    expr = expr.add_suffix('_expression')
    expr = expr.rename(index=str, columns={'TCGAbarcode_expression':'TCGAbarcode'})


    # which ids are in both?                                                                   
    #all_data = pd.merge(meth,expr,on="TCGAbarcode",how="inner")
    #all_data.set_index(all_data["TCGAbarcode"], inplace = True, drop=True)
    
    all_data = meth.merge(expr, on = 'TCGAbarcode', how = 'inner').set_index('TCGAbarcode')
    

    # Need to check variance
    layer_vars = all_data.var( axis = 0)
    if np.any(layer_vars.values == 0):
        print('Removing zero variance predictors')
        layer_mask = layer_vars[layer_vars == 0].index.tolist()
        print('There are %d zero variance columns out of %d columns' %(len(layer_mask), len(all_data.columns)))
        print('These are the zero variance columns: %s' %layer_mask)
        all_data = all_data.loc[:,~all_data.columns.isin(layer_mask)]


    print('Writing all data to disk')
    print('Data Shape:')
    print(all_data.shape)
    all_data.to_csv(input_dragon,sep="\t")

    # subset only the ones we want for dragon
    meth_data = all_data.filter(regex='methylation')
    exp_data = all_data.filter(regex='expression')

    # run dragon                                                                               
    lambdas, lambdas_landscape = estimate_penalty_parameters_dragon(meth_data,exp_data)
    print("Lambdas:")
    print(lambdas)

    newnames = sum([meth_data.columns.tolist(),exp_data.columns.tolist()],[])
    #try except block to catch singular matrix error
    try: 
        r = get_partial_correlation_dragon(meth_data,exp_data,lambdas)
    except np.linalg.LinAlgError:
        print('Matrix is singular, skipping this one')
        r = np.zeros((len(newnames), len(newnames)))

    df = pd.DataFrame(r,columns=newnames,index=newnames)

    df.to_csv(output_dragon,sep="\t")
    

    n = exp_data.shape[0]
    p1 = meth_data.shape[1]
    p2 = exp_data.shape[1]

    print("n, p1, p2:")
    print(n)
    print(p1)
    print(p2)

    # adj_p_vals, p_vals = estimate_p_values_dragon(r, n, p1, p2, lambdas)
    # p_vals_mc = estimate_p_values_mc(r, n, p1, p2, lambdas)

    #df_adj = pandas.DataFrame(adj_p_vals,columns=newnames,index=newnames)
    #df_adj.to_csv(out_dir_long + "/" + uuid + "_dragon_adj_p.tsv",sep="\t")

    #df_raw = pandas.DataFrame(p_vals,columns=newnames,index=newnames)
    #df_raw.to_csv(out_dir_long + "/" + uuid + "_dragon_raw_p.tsv",sep="\t")

    # df_mc = pandas.DataFrame(p_vals_mc,columns=newnames,index=newnames)
    # df_mc.to_csv(out_dir_long + "/" + uuid + "_dragon_mc_p.tsv",sep="\t")




@click.command()
@click.option('-m', '--meth_file', 'meth_file', type=str, required=True,
              help='Path to the methylation file')
@click.option('-e', '--expr_file', 'expr_file', type=str, required=True,
              help='Path to expression file')
@click.option('-o', '--output_dir', 'output_dir', type=str, required=True,
              help='TSV, The output results for dragon are written here')
@click.option('--barcode', type=str, show_default=True, default='TCGAbarcode',
              help='barcode column')
def lioness_dragon(meth_file, expr_file, output_dir, barcode='TCGAbarcode'):
    """Run dragon"""

    # Read methylation data
    aaa = pd.read_csv(meth_file,sep=',', header=0,index_col=0)
    # Read expression data
    bbb = pd.read_csv(expr_file,sep=',',header=0,index_col=0)
    
    # Add suffixes to the column names and merge tables
    aaa = aaa.add_suffix('_meth')
    aaa = aaa.rename(index=str, columns={'TCGAbarcode_meth':'TCGAbarcode'})

    bbb = bbb.add_suffix('_expr')
    bbb = bbb.rename(index=str, columns={'TCGAbarcode_expr':'TCGAbarcode'})
    # merge data
    all_data = aaa.merge(bbb, on = 'TCGAbarcode', how = 'inner').set_index('TCGAbarcode')
    
    
    # Need to check variance
    layer_vars = all_data.var( axis = 0)
    if np.any(layer_vars.values == 0):
        print('Removing zero variance predictors')
        layer_mask = layer_vars[layer_vars == 0].index.tolist()
        print('There are %d zero variance columns out of %d columns' %(len(layer_mask), len(all_data.columns)))
        print('These are the zero variance columns: %s' %layer_mask)
        all_data = all_data.loc[:,~all_data.columns.isin(layer_mask)]

    
    print("Running LIONESS-DRAGON for\n - methylation: %s\n - expression: %s" %(meth_file,expr_file))

    print('Data head')
    print(all_data.head())
    

    try: 
        s = LionessDragon(all_data = all_data, output_dir = output_dir, merge_col =  "TCGAbarcode",ext1 = "_meth",ext2="_expr")
        s.lioness_loop()
    except np.linalg.LinAlgError:
        print('ERROR: Matrix is singular, skipping this one')
        os.mkdir(output_dir)
    except Exception as err:
        print(err)
        print('ERROR: one column has variance zero, skipping this one')
        #os.mkdir(output_dir)
            

cli.add_command(dragon)
cli.add_command(lioness_dragon)

if __name__ == '__main__':
    cli()