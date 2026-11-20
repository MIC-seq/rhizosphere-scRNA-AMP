#!/usr/bin/env python3
"""Freeze source labels and rebuild count matrices without RNA-name substring loss.

No sequences are exported. All original cells/states are retained. Existing source
files are read only; outputs are confined to the configured analysis directory.
"""
import csv
import gzip
import hashlib
import json
import re
import shutil
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CFG = json.loads((ROOT / 'config/analysis.json').read_text())
SRC = Path(CFG['source_root'])
OUT = Path(CFG['output_root'])
assert ROOT == OUT, 'Run the script from its configured analysis root'
for sub in ['inputs', 'inputs/matrices', 'audit', 'tables', 'figures', 'reports', 'logs']:
    (OUT / sub).mkdir(parents=True, exist_ok=True)

manifest = []
def opener(path):
    return gzip.open(str(path), 'rt') if str(path).endswith('.gz') else open(str(path))

def read_table(path):
    with opener(path) as handle:
        return list(csv.DictReader(handle, delimiter='\t'))

def record(path, role):
    h = hashlib.sha256()
    with open(str(path), 'rb') as handle:
        for block in iter(lambda: handle.read(1048576), b''):
            h.update(block)
    manifest.append({'source_path': str(path), 'role': role,
                     'bytes': path.stat().st_size, 'sha256': h.hexdigest()})

def write_table(path, rows, fields=None):
    if not fields:
        fields = list(rows[0]) if rows else ['status']
    with open(str(path), 'w') as handle:
        w = csv.DictWriter(handle, fields, delimiter='\t', lineterminator='\n', extrasaction='ignore')
        w.writeheader()
        w.writerows(rows)

def copy_input(rel, name):
    p = SRC / rel
    record(p, name)
    shutil.copy2(str(p), str(OUT / 'inputs' / name))
    return OUT / 'inputs' / name

statefile = copy_input('results/15_final/source_data/Fig5b_state_UMAP_cells.tsv.gz', 'final_states.tsv.gz')
states = read_table(statefile)
assert len(states) == 1708 and len({r['cell_id'] for r in states}) == 1708
copy_input('results/15_final/source_data/Fig5b_state_key.tsv', 'state_key.tsv')
copy_input('results/15_final/source_data/Fig5e_original_candidate_annotations.tsv', 'existing_candidate_panel.tsv')
copy_input('results/12_integrated_amp_function/figure5_optimized/amp_state_prediction/target_cells_individual_AMP_counts.tsv.gz', 'candidate_counts.tsv.gz')
copy_input('config/samples.tsv', 'original_samples.tsv')
copy_input('results/10_function_annotation/excluded_amp_similar_genes.tsv', 'excluded_amp_similar_genes.tsv')
excluded = {r['gene_id'] for r in read_table(OUT / 'inputs/excluded_amp_similar_genes.tsv')}
states_by_cell = {r['cell_id']: r for r in states}
species_by_cell = {r['cell_id']: r['species'] for r in states}

annfile = SRC / 'results/10_function_annotation/unified_gene_cog_annotation.tsv'
record(annfile, 'full source annotation (no protein sequences exported)')
annotations = {}
exclusion_audit = []
for row in read_table(annfile):
    if row['species'] not in CFG['species']:
        continue
    bio = row['gene_biotype'].strip().lower()
    product = row['product'].strip()
    text = (bio + ' ' + product).lower()
    legacy_exclude = any(s in text for s in ['rrna', 'trna', 'ribosomal rna', 'transfer rna'])
    # Exact biotype, with a narrowly anchored fallback for untyped RNA features.
    is_rna = bio in {'rrna', 'trna', 'ribosomal_rna', 'transfer_rna'}
    if bio in {'', '-', 'na', 'unknown'}:
        is_rna = bool(re.match(r'^(?:\d+(?:\.\d+)?s\s+ribosomal rna|trna[- (]|transfer rna\b)', product.lower()))
    reason = 'AMP_similarity' if row['gene_id'] in excluded else ('RNA_gene' if is_rna else 'retained')
    clean = {k: row.get(k, '') for k in ['gene_id', 'species', 'locus_tag', 'assembly_accession', 'gene_biotype', 'product', 'COG_category', 'COG_description']}
    clean['legacy_retained'] = int(not legacy_exclude and reason == 'retained')
    clean['restored_coding_gene'] = int(legacy_exclude and not is_rna and bio == 'protein_coding' and reason == 'retained')
    exclusion_audit.append({'gene_id': row['gene_id'], 'species': row['species'], 'biotype': bio,
                            'product': product, 'decision': reason,
                            'legacy_substring_excluded': int(legacy_exclude),
                            'restored_coding_gene': clean['restored_coding_gene']})
    if reason == 'retained':
        annotations[row['gene_id']] = clean

emapper = SRC / 'results/10_function_annotation/eggnog/target_reference.emapper.annotations'
record(emapper, 'existing eggNOG KO/pathway annotations')
with emapper.open() as handle:
    header = None
    for line in handle:
        if line.startswith('##') or not line.strip():
            continue
        if line.startswith('#'):
            header = line.lstrip('#').rstrip('\n').split('\t')
            continue
        row = dict(zip(header, line.rstrip('\n').split('\t')))
        gene = row['query']
        if gene in annotations:
            annotations[gene].update({k: row.get(k, '-') for k in ['Preferred_name', 'KEGG_ko', 'KEGG_Pathway', 'KEGG_Module', 'EC', 'eggNOG_OGs']})
for row in annotations.values():
    for key in ['Preferred_name', 'KEGG_ko', 'KEGG_Pathway', 'KEGG_Module', 'EC', 'eggNOG_OGs']:
        row.setdefault(key, '-')

counts = Counter()
for sample in sorted({r['sample'] for r in states}):
    path = SRC / 'results/07_directed_mapping' / sample / (sample + '_gene_cell_counts.tsv.gz')
    record(path, 'source UMI counts')
    with gzip.open(str(path), 'rt') as handle:
        for line in handle:
            gene, barcode, umi = line.rstrip('\n').split('\t')
            cell = sample + '__' + barcode
            if cell in states_by_cell and gene in annotations and annotations[gene]['species'] == species_by_cell[cell]:
                counts[(gene, cell)] += int(umi)

summaries, combined_cells, memberships = [], [], []
curated = {p: group for group, paths in CFG['pathway_groups'].items() for p in paths}
pathnames = {}
with (OUT / 'reference/kegg_pathway_names.tsv').open() as f:
    for line in f:
        parts = line.rstrip('\n').split('\t', 1)
        if len(parts) == 2:
            pathnames[re.sub(r'\D', '', parts[0])] = parts[1]
assert set(curated) <= set(pathnames), 'Missing curated KEGG pathway labels'
write_table(OUT / 'inputs/pathway_catalogue.tsv',
            [{'pathway_id': 'map' + p, 'pathway': pathnames[p], 'group': group} for p, group in curated.items()])

for sp in CFG['species']:
    slug = re.sub(r'[^A-Za-z0-9_.-]+', '_', sp)
    directory = OUT / 'inputs/matrices' / slug
    directory.mkdir(exist_ok=True)
    legacy_dir = SRC / 'results/12_integrated_amp_function/species_matrices' / slug
    metadata_path = legacy_dir / 'cells.tsv'
    record(metadata_path, 'original cell QC and AMP metadata')
    metadata = {r['cell_id']: r for r in read_table(metadata_path)}
    for filename in ['genes.tsv', 'gene_by_cell.mtx.gz']:
        record(legacy_dir / filename, 'legacy matrix sensitivity')
        shutil.copy2(str(legacy_dir / filename), str(directory / ('legacy_' + filename)))
    cell_ids = sorted(r['cell_id'] for r in states if r['species'] == sp)
    genes = sorted({g for g, c in counts if species_by_cell[c] == sp})
    gi, ci = {g: i+1 for i, g in enumerate(genes)}, {c: i+1 for i, c in enumerate(cell_ids)}
    entries = sorted((gi[g], ci[c], n) for (g, c), n in counts.items() if c in ci)
    with gzip.open(str(directory / 'gene_by_cell.mtx.gz'), 'wt') as f:
        f.write('%%MatrixMarket matrix coordinate integer general\n% fixed original Fig5b cells, exact RNA feature exclusion\n')
        f.write('{} {} {}\n'.format(len(genes), len(cell_ids), len(entries)))
        for row in entries:
            f.write('{} {} {}\n'.format(*row))
    write_table(directory / 'genes.tsv', [annotations[g] for g in genes])
    totals, detected = Counter(), Counter()
    for g, c in counts:
        if c in ci:
            totals[c] += counts[g, c]
            detected[c] += 1
    cells = []
    for cell in cell_ids:
        row = dict(metadata[cell])
        row.update(states_by_cell[cell])
        row.update({'retained_gene_umi': totals[cell], 'retained_detected_genes': detected[cell]})
        cells.append(row)
    write_table(directory / 'cells.tsv', cells)
    combined_cells.extend(cells)
    for g in genes:
        ids = {re.sub(r'\D', '', p) for p in annotations[g]['KEGG_Pathway'].split(',')}
        for pid in sorted(ids & set(curated)):
            memberships.append({'species': sp, 'gene_id': g, 'pathway_id': 'map' + pid,
                                'pathway': pathnames[pid], 'group': curated[pid],
                                'gene_name': annotations[g]['Preferred_name'], 'product': annotations[g]['product'],
                                'restored_coding_gene': annotations[g]['restored_coding_gene']})
    # Require exact reproduction of the old matrix on the originally retained genes.
    legacy_genes = read_table(legacy_dir / 'genes.tsv')
    legacy_cells = read_table(legacy_dir / 'cells.tsv')
    mismatches = 0
    with gzip.open(str(legacy_dir / 'gene_by_cell.mtx.gz'), 'rt') as f:
        lines = (line for line in f if not line.startswith('%'))
        dims = next(lines)
        for line in lines:
            a, b, n = line.split()
            key = (legacy_genes[int(a)-1]['gene_id'], legacy_cells[int(b)-1]['cell_id'])
            mismatches += int(counts.get(key, 0) != int(n))
    assert mismatches == 0, '{} legacy count mismatch'.format(sp)
    summaries.append({'species': sp, 'cells': len(cells), 'samples': len({r['sample'] for r in cells}),
                      'observed_genes': len(genes), 'legacy_observed_genes': len(legacy_genes),
                      'restored_coding_observed_genes': sum(annotations[g]['restored_coding_gene'] for g in genes),
                      'genes_with_KO': sum(annotations[g]['KEGG_ko'] not in ['', '-'] for g in genes),
                      'genes_with_pathway': sum(annotations[g]['KEGG_Pathway'] not in ['', '-'] for g in genes),
                      'legacy_nonzero_count_mismatches': mismatches,
                      'total_retained_umi': sum(totals.values())})
write_table(OUT / 'inputs/cells.tsv', combined_cells)
write_table(OUT / 'inputs/pathway_gene_membership.tsv', memberships)
write_table(OUT / 'audit/gene_filter_audit.tsv', exclusion_audit)
write_table(OUT / 'audit/input_summary.tsv', summaries)
write_table(OUT / 'audit/source_manifest.tsv', manifest)
(OUT / 'audit/user_design_confirmation.txt').write_text(CFG['sample_independence'] + '\n')
print(json.dumps(summaries, indent=2))
