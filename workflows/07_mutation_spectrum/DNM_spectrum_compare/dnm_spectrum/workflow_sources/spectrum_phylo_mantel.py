#!/usr/bin/env python3
"""
spectrum_phylo_mantel.py
Mantel test of pairwise somatic mutational-spectrum distance vs pairwise
phylogenetic distance across species. A non-significant Mantel correlation
justifies merging subsocial somatic spectra into a single ancestral baseline
(see baseline_somatic_test.py): apparent social-vs-subsocial shifts are then
not plausibly confounded with phylogenetic similarity.

Inputs:
    --species NAME PATH        one per species (somatic spectrum TSV)
    --tree_newick STR          unrooted topology in Newick (no branch lengths)
    --output_tsv PATH
    --output_pdf PATH
"""
import argparse
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

MUT_TYPES = ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']
BASES = ['A', 'C', 'G', 'T']
CAT7_SHORT = ['C>A', 'C>G', 'C>T_nonCpG', 'C>T_CpG', 'T>A', 'T>C', 'T>G']


# --- Readers ---------------------------------------------------------------

def read_somatic(path):
    data = {}
    with open(path) as f:
        header = f.readline().strip().split('\t')
        for line in f:
            fields = line.strip().split('\t')
            row = dict(zip(header, fields))
            data[(row['mutation_type'], row['context'])] = int(row['count'])
    return data


def build_96vector(counts):
    vec = np.zeros(96, dtype=float)
    i = 0
    for mut in MUT_TYPES:
        for five in BASES:
            for three in BASES:
                vec[i] = counts.get((mut, f'{five}{mut[0]}{three}'), 0)
                i += 1
    return vec


def collapse_to_7(vec96):
    d = vec96.reshape(6, 4, 4)
    ct = MUT_TYPES.index('C>T')
    g = BASES.index('G')
    out = np.zeros(7)
    out[0] = d[MUT_TYPES.index('C>A')].sum()
    out[1] = d[MUT_TYPES.index('C>G')].sum()
    ct_cpg = d[ct, :, g].sum()
    ct_all = d[ct].sum()
    out[2] = ct_all - ct_cpg
    out[3] = ct_cpg
    out[4] = d[MUT_TYPES.index('T>A')].sum()
    out[5] = d[MUT_TYPES.index('T>C')].sum()
    out[6] = d[MUT_TYPES.index('T>G')].sum()
    return out


def normalize(vec):
    s = vec.sum()
    return vec / s if s > 0 else vec


# --- Newick ----------------------------------------------------------------

def parse_newick(s):
    """Minimal Newick parser (no branch lengths, no comments).
    Returns nested tuple ((children, label)). Leaves: ([], label)."""
    s = s.strip().rstrip(';').strip()
    pos = [0]

    def _clade():
        children = []
        if pos[0] < len(s) and s[pos[0]] == '(':
            pos[0] += 1
            children.append(_clade())
            while s[pos[0]] == ',':
                pos[0] += 1
                children.append(_clade())
            if s[pos[0]] != ')':
                raise ValueError(f'Expected ) at pos {pos[0]}: {s}')
            pos[0] += 1
        label_start = pos[0]
        while pos[0] < len(s) and s[pos[0]] not in ',()':
            pos[0] += 1
        label = s[label_start:pos[0]].strip()
        return (children, label)

    return _clade()


def leaf_paths(clade, node_id=[0], path=()):
    children, label = clade
    my_id = node_id[0]; node_id[0] += 1
    here = path + (my_id,)
    if not children:
        return {label: here}
    out = {}
    for c in children:
        out.update(leaf_paths(c, node_id, here))
    return out


def phylo_distance_matrix(newick, species_order):
    tree = parse_newick(newick)
    paths = leaf_paths(tree)
    missing = [sp for sp in species_order if sp not in paths]
    if missing:
        raise ValueError(f'Species missing from tree: {missing}. '
                         f'Tree leaves: {sorted(paths.keys())}')
    n = len(species_order)
    D = np.zeros((n, n), dtype=float)
    for i, a in enumerate(species_order):
        for j, b in enumerate(species_order):
            if i == j:
                continue
            pa, pb = paths[a], paths[b]
            common = 0
            for x, y in zip(pa, pb):
                if x == y:
                    common += 1
                else:
                    break
            D[i, j] = len(pa) + len(pb) - 2 * common
    return D


# --- Distance + Mantel -----------------------------------------------------

def cosine_distance_matrix(prop_matrix):
    """prop_matrix: (n_species, n_features). Returns (n, n) cosine distance."""
    n = prop_matrix.shape[0]
    D = np.zeros((n, n), dtype=float)
    norms = np.linalg.norm(prop_matrix, axis=1)
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            denom = norms[i] * norms[j]
            if denom == 0:
                D[i, j] = 0.0
            else:
                D[i, j] = 1.0 - float(np.dot(prop_matrix[i], prop_matrix[j]) / denom)
    return D


def mantel(D1, D2, n_perm=9999, seed=0):
    rng = np.random.default_rng(seed)
    iu = np.triu_indices_from(D1, k=1)
    u1 = D1[iu]
    u2 = D2[iu]
    r_obs = float(np.corrcoef(u1, u2)[0, 1])
    n = D1.shape[0]
    idx = np.arange(n)
    ge = 0
    for _ in range(n_perm):
        perm = rng.permutation(idx)
        D2p = D2[np.ix_(perm, perm)]
        u2p = D2p[iu]
        r = float(np.corrcoef(u1, u2p)[0, 1])
        if abs(r) >= abs(r_obs):
            ge += 1
    p = (ge + 1) / (n_perm + 1)
    return r_obs, p


# --- Output ----------------------------------------------------------------

def write_matrix_block(f, title, mat, species):
    f.write(f'# {title}\n')
    f.write('species\t' + '\t'.join(species) + '\n')
    for i, sp in enumerate(species):
        row = '\t'.join(f'{mat[i, j]:.6f}' for j in range(len(species)))
        f.write(f'{sp}\t{row}\n')
    f.write('\n')


def write_tsv(path, species, D_phylo, D96, D7, r96, p96, r7, p7, n_perm):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w') as f:
        f.write(f'# Species order: {",".join(species)}\n\n')
        write_matrix_block(f, 'Phylogenetic distance (tree topology, edges)',
                           D_phylo, species)
        write_matrix_block(f, 'Somatic spectrum cosine distance (96-cat proportions)',
                           D96, species)
        write_matrix_block(f, 'Somatic spectrum cosine distance (7-cat proportions)',
                           D7, species)
        f.write('# Mantel test (permutation; Pearson r on upper-triangle vectors)\n')
        f.write('metric\tr\tp_perm\tn_permutations\n')
        f.write(f'96-cat\t{r96:.4f}\t{p96:.4g}\t{n_perm}\n')
        f.write(f'7-cat\t{r7:.4f}\t{p7:.4g}\t{n_perm}\n')


def plot_figure(species, D_phylo, D96, r96, p96, r7, p7, out_pdf):
    fig, axes = plt.subplots(1, 3, figsize=(16, 5.2))

    im0 = axes[0].imshow(D_phylo, cmap='viridis')
    axes[0].set_xticks(range(len(species)))
    axes[0].set_xticklabels(species, rotation=45, ha='right')
    axes[0].set_yticks(range(len(species)))
    axes[0].set_yticklabels(species)
    axes[0].set_title('Phylogenetic distance\n(tree topology, edges)')
    for i in range(len(species)):
        for j in range(len(species)):
            axes[0].text(j, i, f'{D_phylo[i,j]:.0f}',
                         ha='center', va='center',
                         color='white' if D_phylo[i, j] > D_phylo.max() / 2 else 'black',
                         fontsize=8)
    fig.colorbar(im0, ax=axes[0], shrink=0.75)

    im1 = axes[1].imshow(D96, cmap='magma')
    axes[1].set_xticks(range(len(species)))
    axes[1].set_xticklabels(species, rotation=45, ha='right')
    axes[1].set_yticks(range(len(species)))
    axes[1].set_yticklabels(species)
    axes[1].set_title('Somatic spectrum\ncosine distance (96-cat)')
    for i in range(len(species)):
        for j in range(len(species)):
            axes[1].text(j, i, f'{D96[i,j]:.2f}',
                         ha='center', va='center',
                         color='white' if D96[i, j] > D96.max() / 2 else 'black',
                         fontsize=8)
    fig.colorbar(im1, ax=axes[1], shrink=0.75)

    ax = axes[2]
    iu = np.triu_indices_from(D_phylo, k=1)
    xs = D_phylo[iu]
    ys = D96[iu]
    ax.scatter(xs, ys, s=40, color='#333333')
    for k in range(len(iu[0])):
        i, j = iu[0][k], iu[1][k]
        ax.annotate(f'{species[i]}-{species[j]}', (xs[k], ys[k]),
                    fontsize=7, xytext=(3, 3), textcoords='offset points')
    ax.set_xlabel('Phylogenetic distance')
    ax.set_ylabel('Spectrum cosine distance (96-cat)')
    ax.set_title(f'Mantel 96-cat: r={r96:.3f}, p={p96:.3g}\n'
                 f'Mantel  7-cat: r={r7:.3f}, p={p7:.3g}')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    fig.suptitle('Somatic spectrum distance vs phylogenetic distance', fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    os.makedirs(os.path.dirname(out_pdf), exist_ok=True)
    fig.savefig(out_pdf, dpi=200, bbox_inches='tight')
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(
        description='Mantel test: somatic spectrum distance vs phylogenetic distance')
    parser.add_argument('--species', action='append', nargs=2,
                        metavar=('NAME', 'PATH'), required=True)
    parser.add_argument('--tree_newick', required=True)
    parser.add_argument('--output_tsv', required=True)
    parser.add_argument('--output_pdf', required=True)
    parser.add_argument('--n_permutations', type=int, default=9999)
    parser.add_argument('--seed', type=int, default=0)
    args = parser.parse_args()

    species = [p[0] for p in args.species]
    paths = {p[0]: p[1] for p in args.species}

    vec96 = np.stack([build_96vector(read_somatic(paths[sp])) for sp in species])
    vec7 = np.stack([collapse_to_7(v) for v in vec96])
    prop96 = np.stack([normalize(v) for v in vec96])
    prop7 = np.stack([normalize(v) for v in vec7])

    D_phylo = phylo_distance_matrix(args.tree_newick, species)
    D96 = cosine_distance_matrix(prop96)
    D7 = cosine_distance_matrix(prop7)

    r96, p96 = mantel(D_phylo, D96, args.n_permutations, args.seed)
    r7, p7 = mantel(D_phylo, D7, args.n_permutations, args.seed + 1)

    write_tsv(args.output_tsv, species, D_phylo, D96, D7,
              r96, p96, r7, p7, args.n_permutations)
    plot_figure(species, D_phylo, D96, r96, p96, r7, p7, args.output_pdf)

    print(f'TSV -> {args.output_tsv}')
    print(f'PDF -> {args.output_pdf}')
    print(f'Mantel 96-cat: r={r96:.4f}, p={p96:.4g}')
    print(f'Mantel  7-cat: r={r7:.4f}, p={p7:.4g}')


if __name__ == '__main__':
    main()
