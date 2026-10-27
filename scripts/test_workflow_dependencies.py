"""Exercise the actual submission case with a recorder; never call Slurm or R."""
from pathlib import Path
import subprocess

source = Path(__file__).with_name('submit_updated_references_to_merges.sh').read_text()
body = source[source.index('case "$phase" in\n  annotation)\n    a='):]
for phase in ('postqc', 'merges'):
    for action in (('remove', 'none') if phase == 'postqc' else ('none',)):
        prefix = f'''set -euo pipefail
phase={phase}
mode=--dry-run
run_id=synthetic_test
RESIDUAL_QC_ACTION={action}
queue_stage() {{
  printf 'RECORD|%s|%s|%s|%s\\n' "$1" "$3" "$4" "$5" >&2
  printf 'plan_%s' "$1"
}}
'''
        result = subprocess.run(['bash'], input=prefix + body, text=True, capture_output=True, check=True)
        jobs = {}
        for line in result.stderr.splitlines():
            if line.startswith('RECORD|'):
                _, label, dependency, array, exports = line.split('|')
                jobs[label] = (dependency, array, exports)
        if phase == 'postqc':
            for branch in ('vz', 'rl'):
                assert jobs[f'{branch}_resolution_review'][0] == f'plan_{branch}_postqc'
                assert jobs[f'{branch}_marker_review'][0] == f'plan_{branch}_resolution_review'
            assert jobs['residual_postqc'][0] == ('plan_residual_qc_apply' if action == 'remove' else '')
            assert jobs['residual_resolution_review'][0] == 'plan_residual_postqc'
            assert jobs['residual_marker_review'][0] == 'plan_residual_resolution_review'
            assert f'RESIDUAL_QC_ACTION={action}' in jobs['residual_postqc'][2].split(',')
            assert ('residual_qc_apply' in jobs) == (action == 'remove')
        else:
            assert jobs['spearman_validation'][0] == 'plan_nonspatial_merge'
            assert jobs['metaneighbor_validation'][:2] == ('plan_nonspatial_merge', '1-3%1')
            assert jobs['spatial_merge'][0] == 'plan_corrected_broad'
            assert 'cluster_validation/spearman' in jobs['spearman_validation'][2]
            assert 'cluster_validation/metaneighbor' in jobs['metaneighbor_validation'][2]
print('PASS: all three resolution-review dependencies, both residual QC choices, and post-merge Spearman/MetaNeighbor scheduling; no jobs submitted.')
