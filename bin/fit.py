#!/usr/bin/env python

import click
import pandas as pd

@click.group()
def cli():
    pass

@click.command()
@click.argument('filename')
def linear(filename):
    print("fitting a linear model from: ", filename)

@click.command()
@click.option('-o', '--output_table', 'output_table', type=str, required=True,
              help='Output table that merges all the tables')
@click.option('-t', '--tables', 'tables', type=str, required=True,
              help='Input tables')
def merge_tables(output_table, tables):
    """Merge multiple csv into one

    Args:
        output_table (str): output csv
        tables (str): csv tables to be merged
    """
    df = pd.read_csv(tables.split(',')[0])
    for ts in tables.split(',')[1:]:
        temp = pd.read_csv(ts)
        df = pd.concat([df,temp], axis = 0)
    
    # Saving without index
    df.to_csv(output_table, index=False)



cli.add_command(linear)
cli.add_command(merge_tables)

if __name__ == '__main__':
    cli()
