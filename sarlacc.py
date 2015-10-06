#!/usr/bin/env python3
import click

@click.group()
@click.option('--debug/--no-debug', default=False)
def cli(debug):
  click.echo('SARLACC')

@cli.command()
def sync():
  click.echo('Syncing')
