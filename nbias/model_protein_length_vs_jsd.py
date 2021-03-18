"""
The goal is to verify that protein length is not a meaningful confounder
for the J-S distance of proteins to their proteome mean, especially for
the top results (e.g. tRNA synthetase proteins).

Previous work has shown (see notebook "AA distance vs length") that length
and distance are correlated on average. The idea is to check whether the
correlation still holds for particular classes of proteins.

The script works as follow:
- For each assembly, train a model predicting distance from protein length
- Collect accuracy metrics (RMSE):
    - Overall (one per assembly)
    - For various protein types (e.g. tRNA synthetases & controls)

The accuracy metrics are to be compared against overall & controls to conclude.
"""
import gzip
import logging
import os
from os.path import join
from pathlib import Path

import pandas as pd
import numpy as np
from Bio import SeqIO
import torch
from sklearn.metrics import r2_score


logger = logging.getLogger(__name__)


def main():
    logging.basicConfig(level=logging.INFO, format="%(asctime)s (%(levelname)s) %(message)s")

    logger.info('Retrieve assembly information')
    cwd = os.getcwd()
    assemblies_path = join(os.getcwd(), 'data/assemblies.csv')
    assembly_df = pd.read_csv(assemblies_path, index_col='assembly_accession')
    complete_genomes_df = assembly_df[
        assembly_df['assembly_level'] == 'Complete Genome'
    ]

    logger.info('Collect protein lengths')
    protein_lengths = collect_protein_lengths(complete_genomes_df)

    logger.info('Determine protein groups from TIGR results')
    n = 50
    top_tigr_ids, bottom_tigr_ids = get_protein_groups(n)

    logger.info('Train a model for each assembly with complete genome')

    # Data holders
    metrics = {
        'assembly_accession': [],
        'overall_r2': [],
        'top_50_r2': [],
        'bottom_50_r2': [],
        'rand_control_r2': [],
        'model_slope': [],
        'model_intercept': [],
        'model_std': [],
    }

    # Process each assembly
    for i, assembly_id in enumerate(complete_genomes_df.index):
        if i == 0 or (i+1) % 100 == 0:
            logger.info(f'Assembly {i + 1} / {len(complete_genomes_df)}')

        # Get the data
        assembly_lengths = protein_lengths[assembly_id]
        distances_path = os.path.join(
            cwd, f'data/sequences/{assembly_id}/{assembly_id}_amino_acid_distance_to_mean.csv'
        )
        tigr_data_path = os.path.join(
            cwd, f'data/sequences/{assembly_id}/{assembly_id}_tigr.csv.gz'
        )
        distances_df = pd.read_csv(distances_path, index_col='protein_id')
        tigr_data = pd.read_csv(tigr_data_path, index_col='query')

        if len(tigr_data) < n:
            continue

        x, y = make_dataset(assembly_lengths, distances_df)

        max_tries = 5
        for train_try in range(max_tries):
            # Define model
            model = LinearNormal()
            optimizer = torch.optim.Adam(model.parameters(), lr=0.01, weight_decay=5e-4)

            # Split into train / test sets
            train_ix, test_ix = split_train_test(x, y)
            x_train, y_train = x[train_ix], y[train_ix]
            x_test, y_test = x[test_ix], y[test_ix]

            # Train the model
            train(
                model, 
                optimizer, 
                x_train, 
                x_test, 
                y_train, 
                y_test, 
                n_steps=500, 
                batch_size=64, 
            )

            # Collect overall metrics
            overall_r2 = compute_r_squared(model, x_test, y_test)

            if pd.isnull(overall_r2) or overall_r2 < 0.2:
                # Training ocassionally fails, restarts in that case
                logger.warning((
                    f'Training failed for assembly {assembly_id} (try {train_try + 1}): '
                    f'overall R2 = {overall_r2:.4f}'
                ))
                continue

            # Collect metrics on top 50 TIGR results
            x_top_50, y_top_50 = make_dataset(
                assembly_lengths, 
                distances_df,
                tigr_data,
                top_tigr_ids,
            )
            if len(x_top_50) > 1:
                top_50_r2 = compute_r_squared(model, x_top_50, y_top_50)
            else:
                top_50_r2 = np.nan

            # Collect metrics on bottom 50 TIGR results
            x_bottom_50, y_bottom_50 = make_dataset(
                assembly_lengths, 
                distances_df,
                tigr_data,
                bottom_tigr_ids,
            )
            if len(x_bottom_50) > 1:
                bottom_50_r2 = compute_r_squared(model, x_bottom_50, y_bottom_50)
            else:
                bottom_50_r2 = np.nan

            # Collect metrics on random proteins
            rand_control_r2 = np.nan
            if len(x_top_50) > 1:
                random_control_tigr_ids = set(np.random.choice(
                    tigr_data.index.tolist(),
                    size=min(n, len(x_top_50)),
                    replace=False,
                ))
                x_rand_control, y_rand_control = make_dataset(
                    assembly_lengths, 
                    distances_df,
                    tigr_data,
                    random_control_tigr_ids,
                )
                if len(x_rand_control) > 1:
                    rand_control_r2 = compute_r_squared(model, x_rand_control, y_rand_control)

            # Retrieve model parameters
            α = float(model.α.detach().numpy())
            β = float(model.β.detach().numpy())
            σ = float(torch.log(1 + torch.exp(model.s)).detach().numpy())

            # Store results
            metrics['assembly_accession'].append(assembly_id)
            metrics['overall_r2'].append(np.round(overall_r2, 6))
            metrics['top_50_r2'].append(np.round(top_50_r2, 6))
            metrics['bottom_50_r2'].append(np.round(bottom_50_r2, 6))
            metrics['rand_control_r2'].append(np.round(rand_control_r2, 6))
            metrics['model_slope'].append(np.round(α, 6))
            metrics['model_intercept'].append(np.round(β, 6))
            metrics['model_std'].append(np.round(σ, 6))
            break

    logger.info('Training completed')

    logger.info('Exporting results')
    output_path = join(cwd, 'data/confounders/length_model_output.csv')
    output_df = pd.DataFrame.from_dict(metrics)
    output_df.to_csv(output_path, index=False)
    logger.info(f'Results exported to {output_path}')

    logger.info('DONE')


class LinearNormal(torch.nn.Module):
    def __init__(self):
        """
        Define trainable paramters.
        """
        super().__init__()
        self.α = torch.nn.Parameter(torch.randn(()))
        self.β = torch.nn.Parameter(torch.randn(()))
        self.s = torch.nn.Parameter(torch.randn(()))
        
    def forward(self, x):
        m = self.α * x + self.β
        σ = torch.log(1 + torch.exp(self.s))  # softplus to ensure σ > 0 
        return torch.distributions.Normal(m, σ)


def train_one_step(model, optimizer, x_batch, y_batch):
    model.train()
    optimizer.zero_grad()
    loss = compute_loss(model, x_batch, y_batch)
    loss.backward()
    optimizer.step()
    return loss


def compute_loss(model, x, y):
    out_dist = model(x)
    return torch.mean(-out_dist.log_prob(y))


def compute_rmse(model, x_test, y_test):
    model.eval()
    out_dist = model(x_test)
    pred = out_dist.mean
    squared = (pred - y_test)**2
    rmse = torch.sqrt(torch.mean(squared))
    return float(rmse.detach().numpy())


def compute_r_squared(model, x_test, y_test):
    model.eval()
    out_dist = model(x_test)
    pred = out_dist.mean
    return r2_score(
        y_test.detach().numpy().flatten(),
        pred.detach().numpy().flatten(),
    )


def predict(model, x):
    model.eval()
    out_dist = model(x)
    return out_dist.mean, out_dist.stddev


def train(
    model, 
    optimizer, 
    x_train, 
    x_val, 
    y_train, 
    y_val, 
    n_steps, 
    batch_size=64, 
    batch_seed=None,
):
    rs = np.random.RandomState(batch_seed)
    train_ix = np.arange(len(x_train))
    train_losses = []
    val_losses = []
    for step in range(n_steps):
        batch_ix = rs.choice(train_ix, size=batch_size, replace=False)
        train_one_step(model, optimizer, x_train[batch_ix], y_train[batch_ix])


def make_dataset(assembly_lengths, distances_df, tigr_data=None, tigr_ids=None):
    if tigr_data is not None and tigr_ids is not None:
        ids = set(tigr_data.index) & set(tigr_ids)
        if len(ids) == 0:
            return [], []

        relevant_protein_ids = tigr_data.loc[ids]['protein_id'].unique()
        union_ix = sorted(
            set(relevant_protein_ids) &
            set(distances_df.index.tolist()) & 
            set(assembly_lengths.keys())
        )

    else:
        union_ix = sorted(
            set(distances_df.index.tolist()) & 
            set(assembly_lengths.keys())
        )

    if len(union_ix) == 0:
        return [], []

    x = np.log([assembly_lengths[protein_id] for protein_id in union_ix])
    y = np.log(distances_df.loc[union_ix]['distance'].values.tolist())
    
    return (
        torch.from_numpy(x[:,np.newaxis]), 
        torch.from_numpy(y[:,np.newaxis]),
    )


def split_train_test(x, y, test_ratio=0.2, seed=444):
    rs = np.random.RandomState(seed)
    
    ix = np.arange(len(x))
    test_size = int(np.ceil(test_ratio * len(ix)))
    test_ix = np.random.choice(ix, size=test_size, replace=False)
    test_ix_set = set(test_ix)
    
    train_ix = np.array([i for i in ix if i not in test_ix])
    
    return train_ix, test_ix


def get_protein_groups(n=50):
    path = join(os.getcwd(), 'data/tigr_aa_probability_left.csv')
    tigr_results = pd.read_csv(path, index_col='id')

    top_tigr_ids = set(tigr_results.index[:n])
    bottom_tigr_ids = set(tigr_results.index[-n:])

    return top_tigr_ids, bottom_tigr_ids


def collect_protein_lengths(df):
    """
    Returns a dictionary with an entry for each assembly id. 
    This entry is itself a dictionary mapping protein id to protein length.
    """
    cwd = os.getcwd()
    output = {
        assembly_id: {}
        for assembly_id in df.index
    }
    for i, assembly_id in enumerate(df.index):
        if i == 0 or (i+1) % 100 == 0:
            print(f'Assembly {i+1} / {len(df)}')
            
        assembly_output = {}
        fasta_path = join(cwd, f'data/sequences/{assembly_id}/{assembly_id}_protein.faa.gz')
        with gzip.open(fasta_path, 'rt') as f:
            for seq in SeqIO.parse(f, 'fasta'):
                assembly_output[seq.id] = len(seq)
        
        output[assembly_id] = assembly_output
        
    return output


if __name__ == '__main__':
    main()
