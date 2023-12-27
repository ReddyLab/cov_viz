# Archived

This repository has been archived in favor of [ccgr_portal_cov_viz](https://github.com/ReddyLab/ccgr_portal_cov_viz)

## Description

Generate experiment coverage visualization data from the CCGR Portal database.

## Usage

    cov_viz <output directory> <experiment accession id> <assembly name ("GRCH37" or "GRCH38")> [bucket size (2,000,000 default)] [chromosome]

The database connection URL is set using the `DATABASE_URL` environment variable, matching the django environment this may be running in.

## Build

Run `cargo build`

## Installation

Run `cargo install --path .`
