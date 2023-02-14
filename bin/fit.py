#!/usr/bin/env python

import click

@click.command()
@click.argument('filename')
def linear(filename):
    print("fitting a linear model from: ", filename)




if __name__ == '__main__':
    linear()
