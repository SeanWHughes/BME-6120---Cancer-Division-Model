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
import os
from collections import defaultdict

#%%     DEFINE TISSUE

# Define possible cell types
cell_types = {
    'Epithelial A': {
        'mut_rate': 0.50,
        'grow_rate': 0.05263,
        'plasticity': 0.005,
        'death_rate': 0.05,
        'weight': 0.20
    },
    'Epithelial B': {
        'mut_rate': 0.50,
        'grow_rate': 0.03093,
        'plasticity': 0.002,
        'death_rate': 0.03,
        'weight': 0.50
    },
    'Mesenchymal': {
        'mut_rate': 0.80,
        'grow_rate': 0.05263,
        'plasticity': 0.010,
        'death_rate': 0.05,
        'weight': 0.25
    },
    'Satellite': {
        'mut_rate': 1.00,
        'grow_rate': 0.11111,
        'plasticity': 0.020,
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
        self.genes['BRCA1'] = Gene('BRCA1', 'driver', 1*10**-4,
            gene_effects = {
                'death_rate': -0.02,
                'grow_rate': 4.0,
                'plasticity': 2.0,
                'mut_rate': 1.0
            }
        )
        self.genes['KRAS'] = Gene('KRAS', 'driver', 3*10**-4,
            gene_effects = {
                'grow_rate': 2.0,
                'plasticity': 0.4
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
                 gene_mut_rate = 5*10**-5,
                 gene_effects = None):
        """
        Initialize a new gene.
        """

        self.gene_name = gene_name
        self.gene_purpose = gene_purpose
        self.gene_mut_rate = gene_mut_rate
        self.mutated = False
        self.effects = gene_effects or {}

class Tissue:
    def __init__(self, tissue_name: str):
        """
        Initialize a new tissue.
        """
        self.tissue_name = tissue_name
        
#%%     FUNCTIONS

def mutator(cell_list):
    """
    Mutate genes using a cell-level multiplier and
    gene-specific baseline mutation probabilities.
    """

    # Loop through every cell in the population
    for cell in cell_list:

        # Initialize recalculation condition for each cell
        recalc = False

        # Loop through each gene for the cell
        for gene in cell.genes.values():

            # Do not mutate a gene more than once
            if gene.mutated:
                continue

            # RNG roll to determine whether the gene mutates
            if random.random() < cell.mut_rate * gene.gene_mut_rate:
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
    Recalculate a cell's parameters from its current cell type
    and all accumulated mutations.
    """

    # Reset parameters to the current cell type's baseline
    cell.mut_rate = cell.base_params['mut_rate']
    cell.grow_rate = cell.base_params['grow_rate']
    cell.death_rate = cell.base_params['death_rate']
    cell.plasticity = cell.base_params['plasticity']

    # Apply every accumulated mutation
    for gene in cell.genes.values():
        if not gene.mutated or not gene.effects:
            continue

        for parameter, effect in gene.effects.items():
            if hasattr(cell, parameter):
                current_value = getattr(cell, parameter)
                setattr(
                    cell,
                    parameter,
                    current_value * (1 + effect)
                )

    # Keep probability parameters between 0 and 1
    for parameter in (
        'mut_rate',
        'grow_rate',
        'death_rate',
        'plasticity'
    ):
        value = getattr(cell, parameter)
        setattr(cell, parameter, max(0.0, min(1.0, value)))

def mitosis(cell_list):
    """
    Cells divide with probability equal to their effective growth rate.
    Daughters are placed near the parent using random directions
    and variable distances to produce irregular spatial clustering.
    Nondividing cells retain their exact positions.
    """

    new_cell_list = []

    next_id = max(
        (cell.cell_id for cell in cell_list
         if cell.cell_id is not None),
        default=-1
    ) + 1

    divisions = 0

    # Spatial parameters
    mean_offset = 0.03
    min_offset = 0.0
    max_offset = 0.1

    for parent in cell_list:

        # Nondividing cells retain their exact coordinates
        if random.random() >= parent.grow_rate:
            new_cell_list.append(parent)
            continue

        divisions += 1

        parent_x = parent.plot_x
        parent_y = parent.plot_y

        # Give each division a random orientation
        base_angle = random.uniform(0, 2 * np.pi)

        # The daughters tend to separate, but not perfectly symmetrically
        for daughter_index in range(2):

            daughter = copy.deepcopy(parent)

            daughter.parent_id = parent.cell_id
            daughter.cell_id = next_id
            daughter.lineage_root = parent.lineage_root

            # Opposite general directions with random angular variation
            angle = (
                base_angle
                + daughter_index * np.pi
                + random.gauss(0, 0.65)
            )

            # Random distance: most daughters remain close to the parent
            distance = random.triangular(
                min_offset,
                max_offset,
                mean_offset
            )

            # Add small independent spatial noise
            daughter.plot_x = (
                parent_x
                + distance * np.cos(angle)
                + random.gauss(0, 0.025)
            )

            daughter.plot_y = (
                parent_y
                + distance * np.sin(angle)
                + random.gauss(0, 0.025)
            )

            next_id += 1
            new_cell_list.append(daughter)

    cell_list[:] = new_cell_list

    print(
        f"Divisions: {divisions} | "
        f"Population after mitosis: {len(cell_list)}"
    )

    return cell_list

#%%     THERAPY

# Two neutral genes act as drug-resistance genes
RESISTANCE_GENES = ['CYP2D6', 'UGT1A1']

def is_resistant(cell):
    """
    A cell is resistant if either resistance gene is mutated.
    """
    return any(cell.genes[g].mutated for g in RESISTANCE_GENES)


def acquire_resistance(cell_list, rate=2e-4):
    """
    Each resistance gene in each cell may mutate once per timestep.
    """
    for cell in cell_list:
        for g in RESISTANCE_GENES:
            gene = cell.genes[g]
            if not gene.mutated and random.random() < rate:
                gene.mutated = True

def therapy(cell_list, kill=0.40, protection=0.95):
    """
    Kill cells during a treatment step. Resistant cells are mostly protected.
    """
    survivors = []
    for cell in cell_list:
        if is_resistant(cell):
            p_kill = kill * (1 - protection)
        else:
            p_kill = kill
        if random.random() >= p_kill:
            survivors.append(cell)
    cell_list[:] = survivors
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

def report_driver_cells(cell_list, timestep):
    driver_cells = [
        cell for cell in cell_list
        if any(
            gene.mutated and gene.gene_purpose == 'driver'
            for gene in cell.genes.values()
        )
    ]

    if not driver_cells:
        print(f"Timestep {timestep}: No driver-mutated cells.")
        return

    print(
        f"\nTimestep {timestep}: "
        f"{len(driver_cells)} driver-mutated cells"
    )

    for cell in driver_cells[:10]:
        mutated_drivers = [
            gene.gene_name
            for gene in cell.genes.values()
            if gene.mutated and gene.gene_purpose == 'driver'
        ]

        print(
            f"  Cell {cell.cell_id}, "
            f"parent={cell.parent_id}, "
            f"drivers={mutated_drivers}, "
            f"grow_rate={cell.grow_rate:.4f}"
        )

#%%     ANIMATION

plt.ion()

# One persistent figure for the live animation
plot_fig, plot_ax = plt.subplots(figsize=(11, 8))

# Folder for saved snapshots
snapshot_folder = "cell_population_snapshots"
os.makedirs(snapshot_folder, exist_ok=True)

def plot_cell_population(
    cell_list,
    timestep=None,
    tumor_start_timestep=None
):
    """
    Update one live plot every timestep and save a separate
    snapshot every 10 timesteps.
    """

    global plot_fig, plot_ax

    if not cell_list:
        print("The population is empty.")
        return
    
    # Enforce a max cell count for the animation (different than sim max)
    if len(cell_list) > 1000:
        return

    ax = plot_ax
    fig = plot_fig

    # Clear and redraw the live plot
    ax.clear()

    # Assign colors to cell types
    cell_types_present = sorted({
        cell.celltype for cell in cell_list
    })

    colors = {
        ct: plt.get_cmap("tab10")(i)
        for i, ct in enumerate(cell_types_present)
    }

    # Plot cells
    for cell in cell_list:

        if not hasattr(cell, "plot_x") or not hasattr(cell, "plot_y"):
            raise AttributeError(
                f"Cell {cell.cell_id} has no stored plotting coordinates."
            )

        driver_mutated = any(
            gene.mutated and gene.gene_purpose == "driver"
            for gene in cell.genes.values()
        )

        if driver_mutated:
            facecolor = "black"
            edgecolor = colors[cell.celltype]
            size = 180
            linewidth = 2.0
        else:
            facecolor = colors[cell.celltype]
            edgecolor = "black"
            size = 140
            linewidth = 0.6

        ax.scatter(
            cell.plot_x,
            cell.plot_y,
            facecolors=facecolor,
            edgecolors=edgecolor,
            s=size,
            linewidths=linewidth,
            alpha=0.9
        )

    # Cell type legend
    handles = [
        plt.Line2D(
            [0], [0],
            marker="o",
            linestyle="",
            color="white",
            markerfacecolor=colors[ct],
            markeredgecolor="black",
            markersize=10,
            label=f"Type {ct}"
        )
        for ct in cell_types_present
    ]

    handles.extend([
        plt.Line2D(
            [0], [0],
            marker="o",
            linestyle="",
            color="white",
            markerfacecolor="black",
            markeredgecolor="black",
            markersize=10,
            label="Driver-mutated"
        ),
        plt.Line2D(
            [0], [0],
            marker="o",
            linestyle="",
            color="white",
            markerfacecolor="white",
            markeredgecolor="black",
            markersize=10,
            label="Normal cell"
        )
    ])

    ax.legend(
        handles=handles,
        title="Cell type / mutation status"
    )

    # Title
    title = (
        "Cell Population"
        if timestep is None
        else f"Cell Population — Timestep {timestep}"
    )

    ax.set_title(
        f"{title}\nTotal cells: {len(cell_list):,}"
    )

    # Tumor-start annotation
    if (
        tumor_start_timestep is not None
        and timestep is not None
        and timestep >= tumor_start_timestep
    ):
        ax.text(
            0.5, 0.97,
            "TUMOR STARTED",
            transform=ax.transAxes,
            ha="center",
            va="top",
            fontsize=16,
            fontweight="bold",
            color="red",
            bbox=dict(
                facecolor="white",
                edgecolor="red",
                boxstyle="round,pad=0.4",
                alpha=0.9
            )
        )

    # Plot formatting
    ax.set_aspect("equal", adjustable="datalim")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_frame_on(False)

    fig.tight_layout()

    # Save a separate permanent snapshot every 10 timesteps
    if timestep is not None and timestep % 10 == 0:
        filename = os.path.join(
            snapshot_folder,
            f"cell_population_timestep_{timestep:03d}.png"
        )

        fig.savefig(filename, dpi=200, bbox_inches="tight")
        print(f"Saved snapshot: {filename}")

    # Refresh the live animation window
    fig.canvas.draw()
    fig.canvas.flush_events()
    plt.pause(0.1)

#%% INITIALIZE POPULATION

# Initial population count
init_pop = 100

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
    
    cell.plot_x = random.uniform(-0.5, 0.5)
    cell.plot_y = random.uniform(-0.5, 0.5)

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

#%%     SIMULATION SETTINGS

timesteps = 2000            # maximum total steps
therapy_delay = 3          # timesteps between detection and therapy
therapy_duration = 20      # duration of therapy
detect = 200               # clinically detectable driver-mutated cell count
max_cells = 10000           # safety stop for run time

# Randomly select 10 founding lineages
lineage_roots = random.sample(
    range(init_pop),
    min(20, init_pop)
)

# Store lineage sizes and newest-descendant mutation counts over time
lineage_history = {
    root: {'size': [], 'mutations': []}
    for root in lineage_roots
}

total_population_history = []
mean_mutation_history = []

type_counts = {ct: [] for ct in cell_types}         # sensitive cells per type, per step
resistant = []                                      # resistant cells per step
total = []                                          # all cells per step
driver_counts = []                                  # driver-mutated cells per step
type_percentages = {ct: [] for ct in cell_types}    # percentage of each cell type

therapy_start = None
therapy_end = None
therapy2_start = None
therapy2_end = None

#%%     SIMULATION LOOP WITH THERAPY

for t in range(timesteps):
    # Animation
    plot_cell_population(cell_list, t)

    # Kill cells naturally
    reaper(cell_list)

    # Count driver-mutated cells
    driver_count = sum(
        any(
            gene.mutated and gene.gene_purpose == 'driver'
            for gene in c.genes.values()
        )
        for c in cell_list
    )

    # Schedule first therapy delayed after detection
    if therapy_start is None and driver_count >= detect:
        therapy_start = t + therapy_delay
        therapy_end = therapy_start + therapy_duration
        print(
            f"Detection threshold reached at timestep {t}. "
            f"First therapy scheduled for timesteps "
            f"{therapy_start} to {therapy_end - 1}."
        )

    # Schedule second therapy if the tumor reaches detection again
    if (
        therapy_start is not None
        and therapy_end is not None
        and therapy2_start is None
        and t >= therapy_end
        and driver_count >= detect
    ):
        therapy2_start = t + therapy_delay
        therapy2_end = therapy2_start + therapy_duration
        print(
            f"Detection threshold reached again at timestep {t}. "
            f"Second therapy scheduled for timesteps "
            f"{therapy2_start} to {therapy2_end - 1}."
        )

    # Apply therapy during either treatment window
    on = (
        (therapy_start is not None
         and therapy_start <= t < therapy_end)
        or
        (therapy2_start is not None
         and therapy2_start <= t < therapy2_end)
    )

    if on:
        therapy(cell_list)

    # Acquire random mutations
    mutator(cell_list)
    acquire_resistance(cell_list)
    
    # Cell plasticity causes type switching
    cell_alternator(cell_list, cell_types)
    
    # Cell division allows growth
    mitosis(cell_list)
    
    # Track selected lineages
    total_population_history.append(len(cell_list))

    timestep_mutations = []

    for root in lineage_roots:
        lineage_cells = [
            c for c in cell_list
            if c.lineage_root == root
        ]

        lineage_history[root]['size'].append(len(lineage_cells))

        if lineage_cells:
            # Calculate average number of mutated genes in the lineage
            mutation_counts = [
                sum(gene.mutated for gene in cell.genes.values())
                for cell in lineage_cells
            ]

            average_mutations = np.mean(mutation_counts)

            lineage_history[root]['mutations'].append(average_mutations)
        else:
            lineage_history[root]['mutations'].append(None)

    mean_mutation_history.append(
        np.mean(timestep_mutations)
        if timestep_mutations else None
    )

    resistant.append(sum(is_resistant(c) for c in cell_list))

    total.append(sum(
        any(
            gene.mutated and gene.gene_purpose == 'driver'
            for gene in c.genes.values()
        )
        for c in cell_list
    ))
    
    driver_counts.append(driver_count)

    for ct in cell_types:
        type_counts[ct].append(sum(
            c.celltype == ct
            for c in cell_list
        ))

        type_percentages[ct].append(
            100 * sum(c.celltype == ct for c in cell_list)
            / len(cell_list) if cell_list else 0
        )

    if therapy_start is not None and t == therapy_start - 1:
        print(
            f"Before therapy: {resistant[-1]} resistant "
            f"of {total[-1]} cells"
        )

    if not cell_list:
        print(f"Tumor eradicated at step {t}")
        break

    if len(cell_list) > max_cells:
        print(f"Stopped at step {t}: over {max_cells} cells")
        break

#%%     RECURRENCE ANALYSIS

if therapy_start is None:
    print("Detection threshold was never reached; therapy was not started.")

if therapy_start is not None:
    during = total[
        max(0, therapy_start):min(len(total), therapy_end)
    ]
    after = total[min(len(total), therapy_end):]

    remission = min(during, default=total[-1] if total else 0) < detect
    recurred = remission and max(after, default=0) >= detect

    if not remission:
        print("No remission: tumor kept growing through therapy")
    elif recurred:
        print("Remission, then recurrence")
    else:
        print("Remission, no recurrence")

    print(
        f"Resistant cells at end: {resistant[-1]} "
        f"of {total[-1]}"
    )
else:
    print("Recurrence analysis skipped because therapy was not started.")

#%%     PLOT THERAPY RESPONSE

steps = range(len(total))                                       # x-axis: 0, 1, 2, ...
sensitive = [total[i] - resistant[i] for i in steps]            # all non-resistant cells

def floor(values):
    """Replace zeros with 0.5 so the log scale doesn't break."""
    return [max(n, 0.5) for n in values]

def add_therapy(ax):
    """Shade both therapy windows and set a log y-axis."""
    if therapy_start is not None:
        ax.axvspan(
            therapy_start,
            min(therapy_end, len(total)) - 1,
            color='gray',
            alpha=0.15,
            label='First Therapy'
        )

    if therapy2_start is not None:
        ax.axvspan(
            therapy2_start,
            min(therapy2_end, len(total)) - 1,
            color='red',
            alpha=0.12,
            label='Second therapy'
        )

    ax.set_yscale('log')
    ax.set_ylim(bottom=0.5)
    ax.set_ylabel('Cells (log)')

fig, axes = plt.subplots(4, 1, figsize=(8, 13), sharex=True)    # 4 panels, one x-axis

# Panel 1: Tumor Size
ax = axes[0]
ax.plot(steps, floor(total), '-o', color='tab:purple', linewidth=2, label='Total tumor')
ax.axhline(detect, color='gray', linestyle=':', label='Detection limit')
add_therapy(ax)
ax.set_title('1. Total tumor size')
ax.legend(loc='upper left', fontsize=8)

# Panel 2: Therapeutic Resistance
ax = axes[1]
ax.plot(steps, floor(sensitive), '-o', color='tab:blue', label='Sensitive')
ax.plot(steps, floor(resistant), 'k-', linewidth=3, label='Resistant')
add_therapy(ax)
ax.set_title('2. Sensitive vs resistant cells')
ax.legend(loc='upper left', fontsize=8)

# Panel 3: all cells broken down by type
ax = axes[2]
for ct, counts in type_counts.items():
    ax.plot(steps, counts, '-o', markersize=4, label=f'Type {ct}')
add_therapy(ax)
ax.set_title('3. All cells by type')
ax.set_xlabel('Timestep')
ax.legend(loc='upper left', fontsize=8)

# Panel 4: Percentage of total population by cell type
ax = axes[3]
for ct, percentages in type_percentages.items():
    ax.plot(steps, percentages, '-o', markersize=4, label=f'Type {ct}')
add_therapy(ax)
ax.set_yscale('linear')
ax.set_ylim(0, 100)
ax.set_ylabel('Population (%)')
ax.set_title('4. Percentage of total population by cell type')
ax.set_xlabel('Timestep')
ax.legend(loc='upper left', fontsize=8)

plt.tight_layout()
plt.show()

#%%     PLOT LINEAGE SIZE VS AVERAGE MUTATION COUNT

fig_lineage, ax_lineage = plt.subplots(figsize=(11, 8))

# Plot each selected lineage
for root, history in lineage_history.items():
    x = np.array(history['mutations'], dtype=float)
    y = np.array(history['size'], dtype=float)

    # Exclude extinct lineages and invalid values
    valid = np.isfinite(x) & np.isfinite(y) & (y > 0)

    x = x[valid]
    y = y[valid]

    if len(x) >= 2 and len(np.unique(x)) >= 2:
        # Fit a straight line to the lineage data
        slope, intercept = np.polyfit(x, y, 1)

        # Generate points for the fitted line
        x_fit = np.linspace(x.min(), x.max(), 100)
        y_fit = slope * x_fit + intercept

        ax_lineage.plot(
            x_fit,
            y_fit,
            linewidth=1.5,
            alpha=0.7,
            label=f'Lineage {root}'
        )

        # Plot the original observations faintly
        ax_lineage.scatter(
            x,
            y,
            s=10,
            alpha=0.15
        )

ax_lineage.set_xlabel('Average number of mutated genes per cell')
ax_lineage.set_ylabel('Number of cells in lineage')
ax_lineage.set_title('Lineage Size vs Average Mutation Count')
ax_lineage.set_ylim(bottom=0)

# Keep the legend manageable with 100 lineages
ax_lineage.legend(
    loc='upper left',
    bbox_to_anchor=(1.02, 1),
    fontsize=6,
    ncol=2
)

plt.tight_layout()
plt.show()
