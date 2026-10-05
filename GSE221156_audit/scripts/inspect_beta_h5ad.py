#!/usr/bin/env python3
"""Inspect H5AD metadata and matrix layout without loading expression into memory.
Requires h5py. Does not certify integer counts or infer donor assignments.
"""
import argparse
import json
import h5py


def scalar(value):
    if isinstance(value, bytes):
        return value.decode()
    if hasattr(value, 'tolist'):
        return scalar(value.tolist())
    if isinstance(value, (list, tuple)):
        return [scalar(v) for v in value]
    return value


def layout(node):
    result = {'encoding': scalar(node.attrs.get('encoding-type', ''))}
    if isinstance(node, h5py.Dataset):
        result.update(shape=list(node.shape), dtype=str(node.dtype))
    else:
        result.update(keys=list(node.keys()), shape=scalar(node.attrs.get('shape', [])))
        if 'data' in node:
            result['data_dtype'] = str(node['data'].dtype)
            result['stored_values'] = node['data'].shape[0]
    return result


p = argparse.ArgumentParser()
p.add_argument('input')
a = p.parse_args()
with h5py.File(a.input, 'r') as f:
    result = {'top_level_keys': list(f.keys()), 'obs_columns': list(f['obs'].keys()),
              'var_columns': list(f['var'].keys()), 'matrices': {}}
    for key in ('X', 'raw/X'):
        if key in f:
            result['matrices'][key] = layout(f[key])
    if 'layers' in f:
        for key in f['layers']:
            result['matrices']['layers/'+key] = layout(f['layers'][key])
    result['candidate_obs_fields'] = {}
    for key, node in f['obs'].items():
        if any(t in key.lower() for t in ('donor', 'sample', 'library', 'disease', 'chemistry', 'batch', 'cell_type')):
            if isinstance(node, h5py.Group) and 'categories' in node:
                result['candidate_obs_fields'][key] = {'categories': scalar(node['categories'][:])}
            elif isinstance(node, h5py.Dataset):
                result['candidate_obs_fields'][key] = {'first_10': scalar(node[:10])}
    print(json.dumps(result, indent=2))
