# Targeted palette-refresh submission fix

The former submission helper counted the TSV header as a task when renv wrote
a startup status message before the header. This produced arrays of 7, 42 and
4 tasks instead of 6 references, 41 manifest samples and 3 pooled objects.
The updated parser accepts only numeric task IDs in TSV rows and requires a
contiguous 1..N map. It was tested with both renv and here() startup messages.
No R or package upgrades are required, and no resource requests change.

The helper also gives a clearer error for a VZ object that still contains the
retired `Mature UBCs` VZ identity. It does not reinterpret or relabel cells.
Such an object must be checked against the current reviewed-label/mapping
workflow. RL and residual Mature UBCs remain valid identities.

The reference errors remain undiagnosed until their stderr is inspected.
The 40984272 sample logs show partial completed consensus snapshots followed
by stale VZ identities; do not overwrite or remove those snapshots. Use a new
tag after the input and reference failures have been resolved.

Upload this archive to the project's scripts directory, then extract there:

```bash
cd /home/acflint/R/Projects/XeniumFCProject/scripts
tar -xzf plot_palette_refresh_20261020_array_fix.tar.gz
```

Extraction and read-only log inspection should take seconds. Changed R parses
are expected under 10 seconds, shell syntax checks under 5 seconds, and the
supported path-only dry-run under 1 minute. Local syntax and parser tests pass;
the local dry-run cannot verify HPC RDS inputs, which are not available locally.
No HPC jobs were submitted or analysis objects loaded by this patch.

Read the actual reference and pooled-VZ errors before attempting a rerun:

```bash
cd /home/acflint/R/Projects/XeniumFCProject
tail -n +1 \
  logs/Xen_PaletteRefresh_40984267_1.err \
  logs/Xen_PaletteRefresh_40984267_2.err \
  logs/Xen_PaletteRefresh_40984267_6.err \
  logs/Xen_PaletteRefresh_40984277_1.err
```

The scale-replacement messages are plotting notices; the actionable failure
is the subsequent `Error` message. Reference errors cannot be inferred from
the samples' VZ taxonomy error or from accounting state alone.
