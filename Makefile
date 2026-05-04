SHELL := /bin/bash

ENV_NAME := cchfv-seg-env
THREADS ?= 4
PYTHON ?= python3

.PHONY: env run dry test smoke clean help

help:
	@echo "make env   - create conda env"
	@echo "make run   - run full Snakemake workflow"
	@echo "make dry   - dry run to show planned steps"
	@echo "make test  - install Python requirements and run smoke checks"
	@echo "make smoke - run checks without installing packages"
	@echo "make clean - remove work and results"

env:
	conda env create -f env/environment.yml || echo "Env may already exist"
	@echo "Activate with: conda activate $(ENV_NAME)"

run:
	@if [ -z "$$NCBI_EMAIL" ]; then echo "Set NCBI_EMAIL before running"; exit 1; fi
	snakemake -s workflow/Snakefile -c $(THREADS) --printshellcmds

dry:
	snakemake -s workflow/Snakefile -n -c $(THREADS)

test:
	$(PYTHON) -m pip install -r env/requirements.txt
	$(MAKE) smoke

smoke:
	$(PYTHON) -m py_compile analysis/scripts/*.py
	$(PYTHON) tests/validate_example_inputs.py
	$(PYTHON) analysis/scripts/example_qc_plot.py --in data-example/example_counts.tsv --out results-example/example_plot.png
	@echo "Wrote results-example/example_plot.png"

clean:
	rm -rf work results results-example logs .snakemake
