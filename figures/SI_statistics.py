# -*- coding: utf-8 -*-
"""
Created on Fri Oct 31 01:50:37 2025

@author: ryanp
"""
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu, kruskal
import scikit_posthocs as sp
import numpy as np
import pandas as pd
import scipy.stats as stats

bmi_gdf = pd.read_csv("C:/Users/ryanp/Downloads/RESEARCH/Scripts/BugModelingData.csv")
network_subset = bmi_gdf[(bmi_gdf.Small ==1)&(bmi_gdf.season=='spring')&(bmi_gdf.DEV<=20)&(bmi_gdf.has_aml==0)]
network_subset['U_bool'] = (network_subset['unconventional_density'] > 0).astype(int)
network_subset['C_bool'] = (network_subset['conventional_density'] > 0).astype(int)

#################################################################

aml_palette = {'0': '#66c2a5', '1': 'red'}

ecoregions = [
    'Central Appalachians',
    'Ridge and Valley',
    'North Central Appalachians',
    'Northern Allegheny Plateau',
    'Erie Drift Plain',
    'Western Allegheny Plateau'
]
hex_colors = [
    '#8dd3c7',
    '#ffffb3',
    '#bebada',
    '#fb8072',
    '#80b1d3',
    '#fdb462'
]
eco_palette = dict(zip(ecoregions, hex_colors))

##############################################################################

def plot_wilcoxon_panel(
    df,
    group_col,
    metrics=['IBI', 'Richness', 'richEPT', 'SHANdivers'],
    metric_labels=None,
    suptitle=None,
    palette=None,
    order=None,
    group_labels=None,
    xlabel=None,
    figsize=(12, 10)
):
    """
    Create a 2x2 panel plot for 4 metrics with Wilcoxon tests
    """
    if metric_labels is None:
        metric_labels = metrics
    
    fig, axes = plt.subplots(2, 2, figsize=figsize)
    axes = axes.flatten()
    
    for idx, (metric, label) in enumerate(zip(metrics, metric_labels)):
        ax = axes[idx]
        
        # Filter data
        plot_df = df[[group_col, metric]].dropna()
        
        # Create boxplot
        sns.boxplot(x=group_col, y=metric, data=plot_df, palette=palette, 
                   order=order, fliersize=0.5, ax=ax)
        
        if group_labels and order:
            ax.set_xticklabels(group_labels)
        
        ax.set_xlabel(xlabel or group_col, fontsize=11)
        ax.set_ylabel(label, fontsize=11)
        ax.set_title(label, fontsize=12, fontweight='bold')
        
        # Perform Mann–Whitney U test
        groups = plot_df[group_col].unique()
        if len(groups) == 2:
            g1, g2 = groups
            data1 = plot_df.loc[plot_df[group_col] == g1, metric]
            data2 = plot_df.loc[plot_df[group_col] == g2, metric]
            
            stat, pval = mannwhitneyu(data1, data2, alternative='two-sided')
            print(f"{label} - Mann–Whitney U: U={stat:.3f}, p={pval:.3g}")
            
            # Annotate significance
            y_max = plot_df[metric].max()
            y_min = plot_df[metric].min()
            y_offset = (y_max - y_min) * 0.05
            y = y_max + y_offset
            
            x1, x2 = 0, 1
            ax.plot([x1, x1, x2, x2], [y, y + y_offset, y + y_offset, y], 
                   lw=1.5, c='black')
            
            if pval < 0.001:
                text = '***'
            elif pval < 0.01:
                text = '**'
            elif pval < 0.05:
                text = '*'
            else:
                text = 'ns'
            ax.text((x1+x2)/2, y + 1.2*y_offset, text, ha='center', 
                   va='bottom', fontsize=11)
    
    if suptitle:
        fig.suptitle(suptitle, fontsize=14, fontweight='bold', y=0.995)
    
    plt.tight_layout()
    plt.show()

##############################################################################
#differences across categorical groups

def plot_kruskal_panel(
    df,
    group_col,
    metrics=['IBI', 'Richness', 'richEPT', 'SHANdivers'],
    metric_labels=None,
    suptitle=None,
    palette=None,
    order=None,
    group_labels=None,
    xlabel=None,
    rotate_xticks=False,
    figsize=(14, 10)
):
    
    if metric_labels is None:
        metric_labels = metrics
    
    fig, axes = plt.subplots(2, 2, figsize=figsize)
    axes = axes.flatten()
    
    dunn_results = {}
    
    for idx, (metric, label) in enumerate(zip(metrics, metric_labels)):
        ax = axes[idx]
        
        # Filter data
        plot_df = df[[group_col, metric]].dropna()
        
        # Create boxplot
        sns.boxplot(x=group_col, y=metric, data=plot_df, palette=palette, 
                   order=order, fliersize=0.5, ax=ax)
        
        if group_labels and order:
            ax.set_xticklabels(group_labels)
        
        ax.set_xlabel(xlabel or group_col, fontsize=11)
        ax.set_ylabel(label, fontsize=11)
        ax.set_title(label, fontsize=12, fontweight='bold')
        
        if rotate_xticks:
            ax.tick_params(axis='x', rotation=45)
            for tick in ax.get_xticklabels():
                tick.set_ha('right')
        
        # Kruskal-Wallis test
        groups = [group[metric].values for _, group in plot_df.groupby(group_col)]
        stat, p = kruskal(*groups)
        print(f"\n{label} - Kruskal–Wallis: H={stat:.3f}, p={p:.3g}")
        
        # Dunn's test
        dunn = sp.posthoc_dunn(plot_df, val_col=metric, group_col=group_col, 
                               p_adjust='bonferroni')
        dunn_results[metric] = dunn
        print(f"Dunn's test for {label}:")
        print(dunn)
        
        # Find significant pairs
        sig_pairs = np.argwhere(dunn.values < 0.05)
        sig_pairs = [(dunn.index[i], dunn.columns[j], dunn.values[i, j]) 
                     for i, j in sig_pairs if i < j]
        sig_pairs = sorted(sig_pairs, key=lambda x: x[2])
        
        # Annotate significant differences
        y_max = plot_df[metric].max()
        y_min = plot_df[metric].min()
        y_offset = (y_max - y_min) * 0.05
        
        for k, (a, b, pval) in enumerate(sig_pairs):  
            x1 = list(dunn.index).index(a)
            x2 = list(dunn.columns).index(b)
            y = y_max + (k+1)*2*y_offset
            ax.plot([x1, x1, x2, x2], [y, y + y_offset, y + y_offset, y], 
                   lw=1.5, c='black')
            
            if pval < 0.001:
                text = '***'
            elif pval < 0.01:
                text = '**'
            elif pval < 0.05:
                text = '*'
            else:
                text = 'ns'
            ax.text((x1+x2)/2, y + 1.2*y_offset, text, ha='center', 
                   va='bottom', fontsize=10)
    
    if suptitle:
        fig.suptitle(suptitle, fontsize=14, fontweight='bold', y=0.995)
    
    plt.tight_layout()
    plt.show()
    
    return dunn_results

##############################################################################
#differences in taxonomic metrics

# 1. Differences across ecoregion
order = ['North Central Appalachians','Central Appalachians','Erie Drift Plain', 
         'Western Allegheny Plateau','Northern Allegheny Plateau','Ridge and Valley']

dunn_eco = plot_kruskal_panel(
    df=bmi_gdf,
    group_col='US_L3NAME',
    metrics=['IBI', 'Richness', 'richEPT', 'SHANdivers'],
    metric_labels=['IBI', 'Richness', 'RichEPT', 'Shannon Diversity'],
    suptitle='Biodiversity Metrics by Ecoregion',
    palette=eco_palette,
    order=order,
    xlabel='Ecoregion',
    rotate_xticks=True,
    figsize=(16, 11)
)

# 2. Differences by AMD presence
plot_wilcoxon_panel(
    df=bmi_gdf,
    group_col='has_aml',
    metrics=['IBI', 'Richness', 'richEPT', 'SHANdivers'],
    metric_labels=['IBI', 'Richness', 'RichEPT', 'Shannon Diversity'],
    suptitle='Biodiversity Metrics by AMD Presence',
    palette=aml_palette,
    order=[0, 1],
    group_labels=['No AMD Present', 'AMD Present'],
    xlabel='AMD Status'
)

# 3. Differences by stream size
plot_wilcoxon_panel(
    df=bmi_gdf,
    group_col='Small',
    metrics=['IBI', 'Richness', 'richEPT', 'SHANdivers'],
    metric_labels=['IBI', 'Richness', 'RichEPT', 'Shannon Diversity'],
    suptitle='Biodiversity Metrics by Stream Size',
    palette='Set3',
    order=[0, 1],
    group_labels=['Semi-wadeable', 'Wadeable'],
    xlabel='Stream Size'
)

# 4. Differences by season
plot_wilcoxon_panel(
    df=bmi_gdf,
    group_col='season',
    metrics=['IBI', 'Richness', 'richEPT', 'SHANdivers'],
    metric_labels=['IBI', 'Richness', 'RichEPT', 'Shannon Diversity'],
    suptitle='Biodiversity Metrics by Sampling Season',
    palette=['#a6bddb', '#fa9fb5'],
    xlabel='Season'
)

plot_wilcoxon_panel(
    df=bmi_gdf,
    group_col='Small',
    metrics=['DRNAREA', 'DEV', 'unconventional_density', 'conventional_density'],
    metric_labels=['Drainage Area', '% DLC', 'UOGD density', 'COGD density'],
    suptitle='watershed attributes for wadeable and semi-wadeable streams',
    palette='Set1',
    order=[0, 1],
    group_labels=['Semi-wadeable', 'wadeable'],
    xlabel='Stream size'
)


#############################################################################
#differences in functional metrics

def plot_wilcoxon_panel_5(
    df,
    group_col,
    metrics=['CG', 'FC', 'PR', 'SC', 'SH'],
    metric_labels=None,
    suptitle=None,
    palette=None,
    order=None,
    group_labels=None,
    xlabel=None,
    figsize=(15, 9)
):
    """
    Create a 2x3 panel plot for 5 metrics with Wilcoxon tests
    """
    if metric_labels is None:
        metric_labels = metrics
    
    fig, axes = plt.subplots(2, 3, figsize=figsize)
    axes = axes.flatten()
    
    for idx, (metric, label) in enumerate(zip(metrics, metric_labels)):
        ax = axes[idx]
        
        # Filter data
        plot_df = df[[group_col, metric]].dropna()
        
        # Create boxplot
        sns.boxplot(x=group_col, y=metric, data=plot_df, palette=palette, 
                   order=order, fliersize=0.5, ax=ax)
        
        if group_labels and order:
            ax.set_xticklabels(group_labels)
        
        ax.set_xlabel(xlabel or group_col, fontsize=11)
        ax.set_ylabel(label, fontsize=11)
        ax.set_title(label, fontsize=12, fontweight='bold')
        
        # Perform Mann–Whitney U test
        groups = plot_df[group_col].unique()
        if len(groups) == 2:
            g1, g2 = groups
            data1 = plot_df.loc[plot_df[group_col] == g1, metric]
            data2 = plot_df.loc[plot_df[group_col] == g2, metric]
            
            stat, pval = mannwhitneyu(data1, data2, alternative='two-sided')
            print(f"{label} - Mann–Whitney U: U={stat:.3f}, p={pval:.3g}")
            
            # Annotate significance
            y_max = plot_df[metric].max()
            y_min = plot_df[metric].min()
            y_offset = (y_max - y_min) * 0.05
            y = y_max + y_offset
            
            x1, x2 = 0, 1
            ax.plot([x1, x1, x2, x2], [y, y + y_offset, y + y_offset, y], 
                   lw=1.5, c='black')
            
            if pval < 0.001:
                text = '***'
            elif pval < 0.01:
                text = '**'
            elif pval < 0.05:
                text = '*'
            else:
                text = 'ns'
            ax.text((x1+x2)/2, y + 1.2*y_offset, text, ha='center', 
                   va='bottom', fontsize=11)
    
    # Remove the extra subplot (6th position)
    fig.delaxes(axes[5])
    
    if suptitle:
        fig.suptitle(suptitle, fontsize=14, fontweight='bold', y=0.995)
    
    plt.tight_layout()
    plt.show()


def plot_kruskal_panel_5(
    df,
    group_col,
    metrics=['CG', 'FC', 'PR', 'SC', 'SH'],
    metric_labels=None,
    suptitle=None,
    palette=None,
    order=None,
    group_labels=None,
    xlabel=None,
    rotate_xticks=False,
    figsize=(16, 9)
):
    """
    Create a 2x3 panel plot for 5 metrics with Kruskal-Wallis and Dunn tests
    """
    if metric_labels is None:
        metric_labels = metrics
    
    fig, axes = plt.subplots(2, 3, figsize=figsize)
    axes = axes.flatten()
    
    dunn_results = {}
    
    for idx, (metric, label) in enumerate(zip(metrics, metric_labels)):
        ax = axes[idx]
        
        # Filter data
        plot_df = df[[group_col, metric]].dropna()
        
        # Create boxplot
        sns.boxplot(x=group_col, y=metric, data=plot_df, palette=palette, 
                   order=order, fliersize=0.5, ax=ax)
        
        if group_labels and order:
            ax.set_xticklabels(group_labels)
        
        ax.set_xlabel(xlabel or group_col, fontsize=11)
        ax.set_ylabel(label, fontsize=11)
        ax.set_title(label, fontsize=12, fontweight='bold')
        
        if rotate_xticks:
            ax.tick_params(axis='x', rotation=45)
            for tick in ax.get_xticklabels():
                tick.set_ha('right')
        
        # Kruskal-Wallis test
        groups = [group[metric].values for _, group in plot_df.groupby(group_col)]
        stat, p = kruskal(*groups)
        print(f"\n{label} - Kruskal–Wallis: H={stat:.3f}, p={p:.3g}")
        
        # Dunn's test
        dunn = sp.posthoc_dunn(plot_df, val_col=metric, group_col=group_col, 
                               p_adjust='bonferroni')
        dunn_results[metric] = dunn
        print(f"Dunn's test for {label}:")
        print(dunn)
        
        # Find significant pairs
        sig_pairs = np.argwhere(dunn.values < 0.05)
        sig_pairs = [(dunn.index[i], dunn.columns[j], dunn.values[i, j]) 
                     for i, j in sig_pairs if i < j]
        sig_pairs = sorted(sig_pairs, key=lambda x: x[2])
        
        # Annotate significant differences
        y_max = plot_df[metric].max()
        y_min = plot_df[metric].min()
        y_offset = (y_max - y_min) * 0.05
        
        for k, (a, b, pval) in enumerate(sig_pairs):
            x1 = list(dunn.index).index(a)
            x2 = list(dunn.columns).index(b)
            y = y_max + (k+1)*2*y_offset
            ax.plot([x1, x1, x2, x2], [y, y + y_offset, y + y_offset, y], 
                   lw=1.5, c='black')
            
            if pval < 0.001:
                text = '***'
            elif pval < 0.01:
                text = '**'
            elif pval < 0.05:
                text = '*'
            else:
                text = 'ns'
            ax.text((x1+x2)/2, y + 1.2*y_offset, text, ha='center', 
                   va='bottom', fontsize=10)
    
    # Remove the extra subplot (6th position)
    fig.delaxes(axes[5])
    
    if suptitle:
        fig.suptitle(suptitle, fontsize=14, fontweight='bold', y=0.995)
    
    plt.tight_layout()
    plt.show()
    
    return dunn_results


##############################################################################
# FUNCTIONAL FEEDING GROUPS PLOTTING
##############################################################################

# 1. Differences across ecoregion
order = ['North Central Appalachians','Central Appalachians','Erie Drift Plain', 
         'Western Allegheny Plateau','Northern Allegheny Plateau','Ridge and Valley']

dunn_eco_ffg = plot_kruskal_panel_5(
    df=bmi_gdf,
    group_col='US_L3NAME',
    metrics=['CG', 'FC', 'PR', 'SC', 'SH'],
    metric_labels=['Collector-Gatherers', 'Filter-Collectors', 'Predators', 'Scrapers', 'Shredders'],
    suptitle='Functional Feeding Groups by Ecoregion',
    palette=eco_palette,
    order=order,
    xlabel='Ecoregion',
    rotate_xticks=True,
    figsize=(17, 10)
)

# 2. Differences by AMD presence
plot_wilcoxon_panel_5(
    df=bmi_gdf,
    group_col='has_aml',
    metrics=['CG', 'FC', 'PR', 'SC', 'SH'],
    metric_labels=['Collector-Gatherers', 'Filter-Collectors', 'Predators', 'Scrapers', 'Shredders'],
    suptitle='Functional Feeding Groups by AMD Presence',
    palette=aml_palette,
    order=[0, 1],
    group_labels=['No AMD Present', 'AMD Present'],
    xlabel='AMD Status',
    figsize=(15, 9)
)

# 3. Differences by stream size
plot_wilcoxon_panel_5(
    df=bmi_gdf,
    group_col='Small',
    metrics=['CG', 'FC', 'PR', 'SC', 'SH'],
    metric_labels=['Collector-Gatherers', 'Filter-Collectors', 'Predators', 'Scrapers', 'Shredders'],
    suptitle='Functional Feeding Groups by Stream Size',
    palette="Set3",
    order=[0, 1],
    group_labels=['Semi-wadeable', 'Wadeable'],
    xlabel='Stream Size',
    figsize=(15, 9)
)

# 4. Differences by season
plot_wilcoxon_panel_5(
    df=bmi_gdf,
    group_col='season',
    metrics=['CG', 'FC', 'PR', 'SC', 'SH'],
    metric_labels=['Collector-Gatherers', 'Filter-Collectors', 'Predators', 'Scrapers', 'Shredders'],
    suptitle='Functional Feeding Groups by Sampling Season',
    palette=['#a6bddb', '#fa9fb5'],
    xlabel='Season',
    figsize=(15, 9)
)
##############################################################################
#regression analysis

rows = ['DEV', 'unconventional_density', 'conventional_density']
cols1 = ['Richness', 'richEPT', 'SHANdivers', 'IBI']
cols2 = ['CG', 'FC', 'PR', 'SC', 'SH']

def plot_custom_grid(rows, cols, df, title):
    nrows = len(rows)
    ncols = len(cols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(4*ncols, 4*nrows))
    
    if nrows == 1:
        axes = axes[np.newaxis, :]
    if ncols == 1:
        axes = axes[:, np.newaxis]
    
    for i, row_var in enumerate(rows):
        for j, col_var in enumerate(cols):
            ax = axes[i, j]
            
            # Scatter plot with regression line
            sns.regplot(x=col_var, y=row_var, data=df, ax=ax, scatter_kws={'s':20}, line_kws={'color':'red'})
            
            # Calculate Spearman correlation
            rho, p = stats.spearmanr(df[col_var], df[row_var])
            
            # Add box with rho and p-value
            ax.text(0.5, 0.9, f'ρ={rho:.2f}\np={p:.3f}', 
                    horizontalalignment='center', verticalalignment='center',
                    fontsize=12, fontweight='bold',
                    bbox=dict(boxstyle="round,pad=0.3", facecolor='lightgrey', edgecolor='black', alpha=0.8),
                    transform=ax.transAxes)
    
    # Adjust layout
    plt.suptitle(title, fontsize=20)
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.show()

# First grid
plot_custom_grid(rows, cols1, bmi_gdf, "predictors (DLC, UOGD and COGD) vs predictands (taxonomic metrics)")

# Second grid
plot_custom_grid(rows, cols2, bmi_gdf, "predictors (DLC, UOGD and COGD) vs predictands (functional metrics)")