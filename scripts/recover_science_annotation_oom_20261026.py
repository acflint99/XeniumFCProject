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
    'Science_standardize': '41049023',
    'Science_panel': '41049024',
    'reference_validation': '41015200',
}
for label, job in expected.items():
    matches = [r for r in rows if len(r) == 4 and r[:2] == ['annotation', label]]
    if len(matches) != 1 or matches[0][2:] != [job, '1']:
        raise SystemExit('Ledger differs from the reviewed recovery state; stop and inspect: ' + label)

def run(*args):
    result = subprocess.run(args, text=True, capture_output=True)
    if result.returncode:
        raise SystemExit(result.stderr or result.stdout or 'Command failed: ' + str(args))
    return result.stdout.strip()

failed = run('sacct', '-X', '-n', '-P', '-j', '41049023',
             '--format=JobID,State,ExitCode')
if not any(line.split('|')[:3] == ['41049023', 'OUT_OF_MEMORY', '0:125']
           for line in failed.splitlines()):
    raise SystemExit('The prior standardization job is not the expected OOM failure; stop and inspect.')

panel_state = run('sacct', '-X', '-n', '-P', '-j', '41049024',
                  '--format=JobID,State,ExitCode')
if not any(line.split('|')[:3] == ['41049024', 'CANCELLED by 13047', '0:0']
           for line in panel_state.splitlines()):
    raise SystemExit('The prior Science panel job is not the expected canceled dependency job; stop and inspect.')
if run('squeue', '-h', '-j', '41015200', '-o', '%T') != 'PENDING':
    raise SystemExit('Reference validation is no longer pending; stop and inspect.')

queue = run('squeue', '-h', '-u', 'acflint', '-o', '%i|%j')
for line in queue.splitlines():
    job, name = line.split('|', 1)
    if name.startswith('Science') and job != '41049024':
        raise SystemExit('Inspect the other Science job before recovery: ' + line)

recovery = Path(tempfile.mkdtemp(prefix='science_oom_recovery_', dir=run_dir))
shutil.copy2(ledger, recovery / 'jobs.before.tsv')
(recovery / 'resource_decision.txt').write_text(
    'Science standardization job 41049023 was OUT_OF_MEMORY at 64G. '
    'Replacement requests 128G, matching the completed Aldinger standardization job 41015193.\n'
)
print('Recovery record:', recovery, flush=True)
def submit(label, script, dependency=None, memory=None):
    args = ['sbatch', '--parsable', '--hold', '--export=ALL,SCIENCE_REFERENCE_OVERWRITE=true']
    if memory:
        args.extend(['--mem=' + memory])
    if dependency:
        args.append('--dependency=afterok:' + dependency)
    job = run(*(args + [script])).split(';', 1)[0]
    if not job.isdigit():
        raise SystemExit('Unexpected submission response: ' + job)
    with (recovery / 'replacement_jobs.tsv').open('a') as handle:
        handle.write(label + '\t' + job + '\n')
    print(label + ': ' + job + ' (held)', flush=True)
    return job

standardize = submit('Science_standardize', 'scripts/run_science_01_standardize_cell_types.slurm', memory='128G')
panel = submit('Science_panel', 'scripts/run_science_02_subset_gene_panel.slurm', standardize)
run('scontrol', 'update', 'JobId=41015200', 'Dependency=afterok:' + panel)

if ledger.read_text() != original:
    raise SystemExit('Ledger changed concurrently. Replacement jobs remain held; inspect recovery record.')
replacements = {'Science_standardize': standardize, 'Science_panel': panel}
for row in rows:
    if len(row) == 4 and row[0] == 'annotation' and row[1] in replacements:
        row[2] = replacements[row[1]]
with tempfile.NamedTemporaryFile(mode='w', dir=run_dir, prefix='.jobs_oom_recovery_', delete=False) as handle:
    handle.write('\n'.join('\t'.join(row) for row in rows) + '\n')
    staged_ledger = Path(handle.name)
shutil.copymode(ledger, staged_ledger)
os.replace(staged_ledger, ledger)

run('scontrol', 'release', panel)
run('scontrol', 'release', standardize)
print('Recovery submitted. Updated ledger:', ledger)
print(run('squeue', '-u', 'acflint', '-o', '%.24i %.36j %.10T %.12M %.100E %R'))
