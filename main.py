# -*- coding: utf-8 -*-
"""
Created on Wed Sep 30 17:18:28 2026

@author: Sean
"""

#%%     IMPORT MODULES

# Standard modules
import random
import copy
import matplotlib.pyplot as plt
import numpy as np

#%%     DEFINE TISSUE

# Define possible cell types
cell_types = {
    'A': {
        'mut_rate': 0.05,
        'grow_rate': 0.10,
        'plasticity': 0.20,
        'death_rate': 0.05,
        'weight': 0.20
    },
    'B': {
        'mut_rate': 0.05,
        'grow_rate': 0.20,
        'plasticity': 0.10,
        'death_rate': 0.03,
        'weight': 0.50
    },
    'C': {
        'mut_rate': 0.08,
        'grow_rate': 0.15,
        'plasticity': 0.30,
        'death_rate': 0.05,
        'weight': 0.25
    },
    'D': {
        'mut_rate': 0.10,
        'grow_rate': 0.05,
        'plasticity': 0.50,
        'death_rate': 0.10,
        'weight': 0.05
    }
}

#%%     CLASSES

class Cell:
    def __init__(self, celltype, ):
        """
        Initialize a new cell.
        """
        
        # Lineage identification for plotting
        self.cell_id = None
        self.parent_id = None
        self.lineage_root = None
        
        # Validate cell type
        if celltype not in cell_types:
            raise ValueError(f"Unknown cell type: {celltype}")
        
        self.celltype = celltype
        
        # Create a constant reference for WT parameters of this cell type
        self.base_params = {
            'mut_rate': cell_types[celltype]['mut_rate'],
            'grow_rate': cell_types[celltype]['grow_rate'],
            'death_rate': cell_types[celltype]['death_rate'],
            'plasticity': cell_types[celltype]['plasticity']
        }
        
        # Define variables for this particular cell's parameters
        self.mut_rate = cell_types[celltype]['mut_rate']
        self.grow_rate = cell_types[celltype]['grow_rate']
        self.death_rate = cell_types[celltype]['death_rate']
        self.plasticity = cell_types[celltype]['plasticity']

        # Each cell gets its own collection of genes
        self.genes = {}

        # Gene assignments
        self.genes['BRCA1'] = Gene('BRCA1', 'driver', 10^-3,
            gene_effects = {
                'death_rate': -0.02,
                'grow_rate': 0.4,
                'plasticity': 0.10,
                'mut_rate': 0.15
            }
        )
        self.genes['KRAS'] = Gene('KRAS', 'driver', 10^-3,
            gene_effects = {
                'grow_rate': 0.10,
                'plasticity': 0.05
            }
        )

        self.genes['ALB'] = Gene('ALB', 'neutral')
        self.genes['PAH'] = Gene('PAH', 'neutral')
        self.genes['APOA1'] = Gene('APOA1', 'neutral')
        self.genes['TTR'] = Gene('TTR', 'neutral')
        self.genes['TF'] = Gene('TF', 'neutral')
        self.genes['SERPINA1'] = Gene('SERPINA1', 'neutral')
        self.genes['CFTR'] = Gene('CFTR', 'neutral')
        self.genes['SLC2A2'] = Gene('SLC2A2', 'neutral')
        self.genes['SLC5A1'] = Gene('SLC5A1', 'neutral')
        self.genes['KCNJ10'] = Gene('KCNJ10', 'neutral')
        self.genes['AQP1'] = Gene('AQP1', 'neutral')
        self.genes['CA2'] = Gene('CA2', 'neutral')
        self.genes['ALDH2'] = Gene('ALDH2', 'neutral')
        self.genes['ADH1B'] = Gene('ADH1B', 'neutral')
        self.genes['UGT1A1'] = Gene('UGT1A1', 'neutral')
        self.genes['CYP2C19'] = Gene('CYP2C19', 'neutral')
        self.genes['CYP2D6'] = Gene('CYP2D6', 'neutral')
        self.genes['OTC'] = Gene('OTC', 'neutral')


class Gene:
    def __init__(self, gene_name: str, 
                 gene_purpose: str, 
                 gene_mut_rate = 10^-6,
                 gene_effects = None):
        """
        Initialize a new gene.
        """
        
        self.gene_name = gene_name
        self.gene_purpose = gene_purpose
        self.gene_mut_rate = gene_mut_rate
        self.mutated = False
        self.gene_effects = gene_effects or {}

class Tissue:
    def __init__(self, tissue_name: str):
        """
        Initialize a new tissue.
        """
        self.tissue_name = tissue_name

def mutator(cell_list):
    """
    Mutate genes using a cell-level multiplier and
    gene-specific baseline mutation probabilities.
    """
    
    # Initialize recalculation condition
    recalc = False
    
    # Loop through every cell in the population
    for cell in cell_list:
        
        # Loop through each gene for the cell
        for gene in cell.genes.values():

            # Do not mutate a gene more than once
            if gene.mutated:
                continue

            # RNG roll to determine whether the gene mutates
            if random.random() < cell.mut_rate*gene.gene_mut_rate:
                gene.mutated = True
           
                # Recalculation required if mutated gene has effects
                if gene.effects:
                    recalc = True

        # Recalculate cell parameters for the cell if any genes with effects were mutated
        if recalc:
            recalculate_parameters(cell)

def cell_alternator(cell_list, cell_types):
    """
    Switch cells to a different cell type based on their
    current plasticity.
    """
    
    celltype_names = list(cell_types.keys())

    for cell in cell_list:

        # Determine whether this cell changes type
        if random.random() < cell.plasticity:

            # Exclude the current cell type from possible type-switching
            possible_types = [
                ct for ct in celltype_names
                if ct != cell.celltype
            ]

            # Weight possible destination types
            weights = [
                cell_types[ct]['weight']
                for ct in possible_types
            ]

            # Randomly select a different cell type
            new_type = random.choices(possible_types, weights=weights, k=1)[0]

            # Update cell type
            cell.celltype = new_type
            params = cell_types[new_type]
    
            # Update the baseline parameters for the new cell type
            cell.base_params = {
                key: params[key]
                for key in (
                    'mut_rate',
                    'grow_rate',
                    'death_rate',
                    'plasticity'
                )
            }
    
            # Reapply all accumulated mutation effects
            recalculate_parameters(cell)

def recalculate_parameters(cell):
    """
    Recalculate effective cell parameters using percentage
    effects from all accumulated gene mutations.
    """

    # Reset to baseline parameters
    for parameter, baseline in cell.base_params.items():
        setattr(cell, parameter, baseline)

    # Apply the effects of all mutated genes in the cell
    for gene in cell.genes.values():
        if gene.mutated:
            for parameter, effect in gene.effects.items():
                current_value = getattr(cell, parameter)
                setattr(
                    cell,
                    parameter,
                    current_value*(1 + effect)
                )

    # Force all probability parameters to be bounded [0, 1]
    for parameter in cell_types[cell.celltype]:
        if parameter != 'weight':
            value = getattr(cell, parameter)
            setattr(cell, parameter, max(0.0, min(1.0, value)))
            
def mitosis(cell_list):
    """
    Replace every parent cell with 2-5 daughter cells depending on growth rate.
    """

    new_cell_list = []

    # Continue IDs after the current maximum
    next_id = max(
        (cell.cell_id for cell in cell_list if cell.cell_id is not None),
        default=-1
    ) + 1

    daughter_counts = [2, 3, 4, 5]

    for parent in cell_list:
        growth = max(0.0, parent.grow_rate)
        weights = [1 + growth * (n - 2) for n in daughter_counts]

        n_daughters = random.choices(
            daughter_counts, weights=weights, k=1
        )[0]

        for _ in range(n_daughters):
            daughter = copy.deepcopy(parent)

            daughter.parent_id = parent.cell_id
            daughter.cell_id = next_id
            next_id += 1

            # Preserve the founding ancestor
            daughter.lineage_root = parent.lineage_root

            new_cell_list.append(daughter)

    cell_list[:] = new_cell_list
    return cell_list

def reaper(cell_list):
    """
    Remove cells based on their individual death rates.
    """
    surviving_cells = []

    for cell in cell_list:

        # Determine whether the cell survives
        if random.random() >= cell.death_rate:
            surviving_cells.append(cell)

    # Replace the population with surviving cells
    cell_list[:] = surviving_cells

    return cell_list

#%%     VISUALIZATION

def plot_cell_population(cell_list, timestep=None, seed=42):
    """
    Plot cells in 2D, clustering descendants by founding ancestor.
    Color indicates cell type.
    """
    if not cell_list:
        print("The population is empty.")
        return

    rng = np.random.default_rng(seed)

    cell_types = sorted({cell.celltype for cell in cell_list})
    colors = {
        ct: plt.get_cmap("tab10")(i)
        for i, ct in enumerate(cell_types)
    }

    # Group cells by their founding ancestor
    lineages = {}
    for cell in cell_list:
        root = cell.lineage_root
        lineages.setdefault(root, []).append(cell)

    fig, ax = plt.subplots(figsize=(11, 8))

    # Place each lineage in its own cluster
    for cluster_index, (root, cells) in enumerate(
        sorted(lineages.items())
    ):
        # Arrange lineages on a grid
        columns = max(1, int(np.ceil(np.sqrt(len(lineages)))))
        center_x = (cluster_index % columns) * 4.0
        center_y = -(cluster_index // columns) * 4.0

        # Arrange descendants around their lineage center
        for cell in cells:
            x = center_x + rng.normal(0, 0.55)
            y = center_y + rng.normal(0, 0.55)

            ax.scatter(
                x, y,
                color=colors[cell.celltype],
                s=65,
                edgecolors="black",
                linewidths=0.4,
                alpha=0.85
            )

    # Legend for cell types
    handles = [
        plt.Line2D(
            [0], [0],
            marker="o",
            linestyle="",
            color=colors[ct],
            markeredgecolor="black",
            markersize=8,
            label=f"Type {ct}"
        )
        for ct in cell_types
    ]

    ax.legend(handles=handles, title="Cell type")
    ax.set_title(
        "Cell Population" if timestep is None
        else f"Cell Population — Timestep {timestep}"
    )
    ax.set_aspect("equal", adjustable="datalim")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_frame_on(False)

    plt.tight_layout()
    plt.show()

#%% INITIALIZE POPULATION

# Initial population count
init_pop = 10

# Define a list for every cell in the population
cell_list = []

# Create lists of cell type names and their corresponding weights
celltype_names = list(cell_types.keys())
celltype_weights = [
    cell_types[celltype]['weight']
    for celltype in celltype_names
]

# Cell population initialization loop
for i in range(init_pop):

    # Choose one cell type randomly according to its weight
    celltype = random.choices(
        celltype_names,
        weights=celltype_weights,
        k=1
    )[0]

    # Create the cell from the celltype and assign unique cell IDs
    cell = Cell(celltype)
    cell.cell_id = i
    cell.lineage_root = i

    # Add new cell to the population list
    cell_list.append(cell)

#%%     CHECK INITIAL POPULATION

# Output tissue population demographics
for celltype in celltype_names:

    count = sum(
        cell.celltype == celltype
        for cell in cell_list
    )

    print(f"Cell type {celltype}: {count}")


#%%     SIMULATION LOOP

# Define total number of timesteps for simulation
timesteps = 3

# Timestep marching
for t in range(timesteps):
    reaper(cell_list)
    mutator(cell_list)
    cell_alternator(cell_list, cell_types)
    mitosis(cell_list)
    plot_cell_population(cell_list, timestep=t)