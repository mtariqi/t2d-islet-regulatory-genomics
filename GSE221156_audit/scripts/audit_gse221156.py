#!/usr/bin/env python3
# ==============================================================================
# Project:      T2D Islet Regulatory Genomics: donor-level beta-cell transcriptomics in type 2 diabetes
# Script:       GSE221156_audit/scripts/audit_gse221156.py
# Description:  Audit GEO GSE221156 records: libraries, donors, pooled-library mappings and supplementary files.
#
# Author:       Md Tariqul Islam
# Affiliation:  Northeastern University
# ORCID:        https://orcid.org/0009-0009-6545-8040
# Repository:   https://github.com/mtariqi/t2d-islet-regulatory-genomics
# License:      See LICENSE in the repository root
# Citation:     See CITATION.cff in the repository root
# Version:      Tracked in git; cite the tagged release used for analysis
# ==============================================================================

"""Audit saved GEO records. Does not download expression data or infer cell donors."""
import argparse
import csv
import hashlib
import json
import re
from collections import Counter
from pathlib import Path


def write_tsv(path, rows, fields):
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def audit(source, out):
    out.mkdir(parents=True, exist_ok=True)
    text = (source / 'samples.txt').read_text()
    libraries, links, files, donors = [], [], [], {}
    stable_fields = {'disease_status_source': 'disease state', 'sex': 'sex',
                     'age': 'age', 'BMI': 'bmi', 'HbA1c': 'hba1c',
                     'reported_race_source': 'race'}
    for block in text.split('^SAMPLE = ')[1:]:
        lines = block.splitlines()
        gsm = lines[0].strip()
        values, char = {}, {}
        for line in lines[1:]:
            if ' = ' not in line:
                continue
            key, value = line.split(' = ', 1)
            values.setdefault(key, []).append(value)
            if key == '!Sample_characteristics_ch1':
                ck, cv = value.split(':', 1)
                char[ck.lower().strip()] = cv.strip()
        title = values['!Sample_title'][0]
        if not re.fullmatch(r'Islet\d+(?:/\d+)?(?:_MS\d+)?', title):
            raise ValueError(f'Unexpected donor title: {title}')
        donor_ids = ['Islet' + n for n in re.findall(r'\d+', title.split('_')[0])]
        pooled = len(donor_ids) > 1
        library = char['library id']
        supplementary = [v for k, vs in values.items()
                         if k.startswith('!Sample_supplementary_file') for v in vs]
        types = []
        for url in supplementary:
            kind = next((k for k in ('barcodes', 'features', 'matrix')
                         if ('_' + k + '.') in url), 'other')
            types.append(kind)
            files.append({'gsm': gsm, 'library_id': library, 'file_type': kind,
                          'file_name': url.rsplit('/', 1)[-1],
                          'url': url.replace('ftp://', 'https://'),
                          'byte_size': '', 'sha256': '',
                          'validation_status': 'listed_only_not_downloaded'})
        if sorted(types) != ['barcodes', 'features', 'matrix']:
            raise ValueError(f'Incomplete triplet: {gsm}')
        libraries.append({'gsm': gsm, 'library_id': library, 'title': title,
                          'donor_ids': ';'.join(donor_ids), 'pooled': str(pooled).lower(),
                          'chemistry': char['chemistry'],
                          'disease_status_source': char['disease state'],
                          'center_source': char.get('center', ''),
                          'assembly_source': ';'.join(v.removeprefix('Assembly: ')
                              for v in values.get('!Sample_data_processing', [])
                              if v.startswith('Assembly: ')),
                          'sra_relation': ';'.join(v for v in values.get('!Sample_relation', [])
                                                  if v.startswith('SRA: '))})
        for i, donor in enumerate(donor_ids):
            row = {'donor_id': donor}
            for field, key in stable_fields.items():
                parts = char[key].split('/')
                if len(parts) != len(donor_ids):
                    raise ValueError(f'Cannot resolve {key} for {gsm}')
                row[field] = parts[i]
            row['disease_status'] = {'NonDiabetic': 'ND', 'PreDiabetic': 'PD',
                                    'Type2Diabetic': 'T2D'}[row['disease_status_source']]
            if donor in donors and donors[donor] != row:
                raise ValueError(f'Conflicting donor metadata: {donor}')
            donors[donor] = row
            links.append({'gsm': gsm, 'library_id': library, 'donor_id': donor,
                          'mapping_level': 'library_membership_only',
                          'cell_assignment_status': 'requires_demultiplexing_map' if pooled
                          else 'single_donor_library_pending_cell_QC'})
    groups = Counter(d['disease_status'] for d in donors.values())
    if len(libraries) != 54 or len(donors) != 48 or groups != Counter(ND=17, PD=14, T2D=17):
        raise ValueError(f'Unexpected cohort: {len(libraries)}, {len(donors)}, {groups}')
    for name, rows in [('GSE221156_libraries.tsv', libraries),
                       ('GSE221156_donors.tsv', sorted(donors.values(), key=lambda d: int(d['donor_id'][5:]))),
                       ('GSE221156_library_donor_links.tsv', links),
                       ('GSE221156_geo_file_inventory.tsv', files)]:
        write_tsv(out / name, rows, list(rows[0]))
    assets = []
    collection = json.loads((source / 'cellxgene.json').read_text())
    for dataset in collection['datasets']:
        for asset in dataset.get('assets', []):
            assets.append({'dataset_id': dataset['dataset_id'], 'title': dataset['title'],
                           'cell_count': dataset['cell_count'], 'file_type': asset['filetype'],
                           'byte_size': asset['filesize'], 'url': asset['url'],
                           'payload_status': 'not_downloaded_or_inspected'})
    write_tsv(out / 'GSE221156_cellxgene_assets.tsv', assets, list(assets[0]))
    summary = {'geo_libraries': len(libraries), 'unique_donors': len(donors),
               'donor_groups': dict(groups), 'pooled_libraries': sum(l['pooled']=='true' for l in libraries),
               'single_donor_libraries': sum(l['pooled']=='false' for l in libraries),
               'listed_geo_files': len(files), 'chemistry_library_counts': dict(Counter(l['chemistry'] for l in libraries)),
               'expression_payloads_inspected': False,
               'source_sha256': {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(source.iterdir()) if p.is_file()}}
    (out / 'GSE221156_audit_summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps(summary, indent=2))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--source-dir', type=Path, required=True)
    parser.add_argument('--outdir', type=Path, required=True)
    args = parser.parse_args()
    audit(args.source_dir, args.outdir)
