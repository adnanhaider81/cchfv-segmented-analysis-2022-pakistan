# Contributing

Thanks for improving this analysis.

## Ground rules
- Do not commit patient data or restricted datasets.
- Keep parameters in `config/config.yaml` for reproducibility.
- Match tool versions in `env/environment.yml` where possible.
- Use neutral sample names in examples, documentation, and public outputs.

## Workflow
1. Fork and create a feature branch.
2. Make focused changes with clear commit messages.
3. Dry run before opening a PR:
   ```bash
   conda activate cchfv-seg-env
   export NCBI_EMAIL="your@email"
   make dry
   ```
4. Verify `make test` works.
5. Open a PR describing your change and parameter updates.

## Code style
- Python scripts under `analysis/scripts` should use argparse and `--help` text.
- Prefer TSV for tables and FASTA for sequences.
- Avoid hard coded paths. Use flags and config.
