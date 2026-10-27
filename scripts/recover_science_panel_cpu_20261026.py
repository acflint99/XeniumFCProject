from pathlib import Path
import os
import shutil
import subprocess
import tempfile

root = Path('/home/acflint/R/Projects/XeniumFCProject')
os.chdir(root)
run_dir = root / 'outputs/reference_label_refresh_runs/updated_reference_labels_20261022'
ledger = run_dir / 'jobs.tsv'
original = ledger.read_text()
rows = [line.split('\t') for line in original.splitlines()]
expected = {
    'Science_standardize': '41049984',
    'Science_panel': '41049985',
    'reference_validation': '41015200',
}
for label, job in expected.items():
    matches = [r for r in rows if len(r) == 4 and r[:2] == ['annotation', label]]
    if len(matches) != 1 or matches[0][2:] != [job, '1']:
        raise SystemExit('Ledger differs from the reviewed panel-recovery state; stop and inspect: ' + label)

def run(*args):
    result = subprocess.run(args, text=True, capture_output=True)
    if result.returncode:
        raise SystemExit(result.stderr or result.stdout or 'Command failed: ' + str(args))
    return result.stdout.strip()

standardize_state = run('sacct', '-X', '-n', '-P', '-j', '41049984',
                        '--format=JobID,State,ExitCode')
if not any(line.split('|')[:3] == ['41049984', 'COMPLETED', '0:0']
           for line in standardize_state.splitlines()):
    raise SystemExit('Science standardization is not the expected successful job; stop and inspect.')

panel_state = run('sacct', '-X', '-n', '-P', '-j', '41049985',
                  '--format=JobID,State,ExitCode')
if not any(line.split('|')[:3] == ['41049985', 'FAILED', '2:0']
           for line in panel_state.splitlines()):
    raise SystemExit('Science panel is not the expected failed dry-run job; stop and inspect.')

if run('squeue', '-h', '-j', '41015200', '-o', '%T') != 'PENDING':
    raise SystemExit('Reference validation is no longer pending; stop and inspect.')

queue = run('squeue', '-h', '-u', 'acflint', '-o', '%i|%j')
for line in queue.splitlines():
    job, name = line.split('|', 1)
    if name.startswith('Science'):
        raise SystemExit('Inspect the other Science job before panel recovery: ' + line)

recovery = Path(tempfile.mkdtemp(prefix='science_panel_cpu_recovery_', dir=run_dir))
shutil.copy2(ledger, recovery / 'jobs.before.tsv')
(recovery / 'runtime_decision.txt').write_text(
    'Science panel job 41049985 failed during yaml config parsing on an unconstrained node. '
    'Replacement uses the corrected launcher with an Intel-node constraint.\n'
)
print('Recovery record:', recovery, flush=True)

args = [
    'sbatch', '--parsable', '--hold', '--export=ALL,SCIENCE_REFERENCE_OVERWRITE=true',
    '--constraint=intel', 'scripts/run_science_02_subset_gene_panel.slurm'
]
panel = run(*args).split(';', 1)[0]
if not panel.isdigit():
    raise SystemExit('Unexpected submission response: ' + panel)
(recovery / 'replacement_jobs.tsv').write_text('Science_panel\t' + panel + '\n')
print('Science_panel:', panel + ' (held)', flush=True)

run('scontrol', 'update', 'JobId=41015200', 'Dependency=afterok:' + panel)
if ledger.read_text() != original:
    raise SystemExit('Ledger changed concurrently. Replacement panel remains held; inspect recovery record.')
for row in rows:
    if len(row) == 4 and row[:2] == ['annotation', 'Science_panel']:
        row[2] = panel
with tempfile.NamedTemporaryFile(mode='w', dir=run_dir, prefix='.jobs_panel_recovery_', delete=False) as handle:
    handle.write('\n'.join('\t'.join(row) for row in rows) + '\n')
    staged_ledger = Path(handle.name)
shutil.copymode(ledger, staged_ledger)
os.replace(staged_ledger, ledger)
run('scontrol', 'release', panel)
print('Recovery submitted. Updated ledger:', ledger)
print(run('squeue', '-u', 'acflint', '-o', '%.24i %.36j %.10T %.12M %.100E %R'))
