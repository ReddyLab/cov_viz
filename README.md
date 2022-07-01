## Description
Generate experiment coverage visualization data from CCGR Portal database.

## Usage

    ./cov_viz <output directory> <experiment accession id> <assembly name ("GRCH37" or "GRCH38")> [bucket size (2,000,000 default)] [chromosome]


The database connection URL is set using the `DATABASE_URL` environment variable, matching the django environment this may be running in.	
