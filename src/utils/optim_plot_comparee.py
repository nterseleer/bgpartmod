"""
Module pour comparer les paramètres optimisés entre plusieurs simulations d'optimisation.
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Dict, List, Optional, Any
from datetime import datetime
from sklearn.preprocessing import MinMaxScaler

from src.utils import optimization as optim
from src.utils import functions as fns
from src.config_model import varinfos
from src.config_system import path_config

# Constants
FIGURE_PATH = path_config.FIGURE_PATH

# Style configuration
plt.style.use('seaborn-v0_8-whitegrid')
sns.set_palette("husl")


def compare_optimizations_parameters(optimizations_dict: Dict[str, str],
                                   ncols: int = 3,
                                   figsize: Optional[tuple] = None,
                                   alpha: float = 0.7,
                                   savefig: bool = False,
                                   filename: Optional[str] = None,
                                   dateinname: bool = False) -> tuple:
    """
    Compare les paramètres optimisés entre plusieurs simulations d'optimisation.
    Génère deux graphiques : un par paramètre et un avec paramètres normalisés.
    
    Args:
        optimizations_dict: Dictionnaire avec pseudo-nom (clé) -> nom officiel OPT<number> (valeur)
        ncols: Nombre de colonnes pour l'arrangement des subplots
        figsize: Taille de la figure (largeur, hauteur). Si None, calculée automatiquement
        alpha: Transparence des boxplots
        savefig: Si True, sauvegarde la figure
        filename: Nom du fichier de sauvegarde (sans extension)
        dateinname: Si True, ajoute la date au nom du fichier
    
    Returns:
        Tuple de (figure1, figure2) avec les deux graphiques
    """
    
    # Charger les données
    all_best_params, parameters = _load_optimization_data(optimizations_dict)
    
    # Graphique 1: Boxplot par paramètre (amélioré)
    fig1 = _create_individual_parameter_plots(
        all_best_params, parameters, optimizations_dict,
        ncols, figsize, alpha, savefig, filename, dateinname
    )
    
    # Graphique 2: Boxplot avec paramètres normalisés
    fig2 = _create_normalized_comparison_plot(
        all_best_params, parameters, optimizations_dict,
        savefig, filename, dateinname
    )
    
    return fig1, fig2


def _load_optimization_data(optimizations_dict: Dict[str, str]) -> tuple:
    """Charge les données d'optimisation et extrait les paramètres."""
    print("Chargement des optimisations...")
    all_best_params = {}  # pseudo_name -> {param: value}
    all_parameters = set()
    
    for pseudo_name, opt_name in optimizations_dict.items():
        try:
            # Charger l'optimisation
            opt = optim.Optimization.load_existing(opt_name)
            
            # Traiter les résultats pour obtenir les meilleurs paramètres
            if not hasattr(opt, 'summary'):
                opt.process_results()
            
            # Stocker les meilleurs paramètres
            all_best_params[pseudo_name] = opt.summary['best_parameters']
            all_parameters.update(opt.summary['best_parameters'].keys())
            
            print(f"✓ {pseudo_name} ({opt_name}): {opt.summary['best_score']:.3f}")
            
        except Exception as e:
            print(f"✗ Erreur lors du chargement de {pseudo_name} ({opt_name}): {e}")
            continue
    
    if not all_best_params:
        raise ValueError("Aucune optimisation n'a pu être chargée")
    
    parameters = sorted(list(all_parameters))
    return all_best_params, parameters


def _create_individual_parameter_plots(all_best_params, parameters, optimizations_dict,
                                     ncols, figsize, alpha, savefig, filename, dateinname):
    """Crée le graphique amélioré avec un subplot par paramètre - version corrigée."""
    
    # Forcer une meilleure disposition : maximum 2 colonnes pour plus d'espace
    ncols = min(ncols, 3)
    nrows = int(np.ceil(len(parameters) / ncols))
    
    if figsize is None:
        figsize = (10 * ncols, 6 * nrows)  # Plus grand pour éviter le surcharge
    
    # Créer la figure avec style moderne
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, facecolor='white')
    fig.patch.set_facecolor('white')
    
    if nrows == 1 and ncols == 1:
        axes = np.array([axes])
    elif nrows == 1 or ncols == 1:
        axes = axes.flatten() if hasattr(axes, 'flatten') else np.array([axes])
    else:
        axes = axes.flatten()
    
    # Palette de couleurs distinctes et lisibles
    colors = sns.color_palette("tab10", len(optimizations_dict))  # Plus distinctes
    color_map = dict(zip(optimizations_dict.keys(), colors))
    
    # Créer une légende globale une seule fois
    legend_elements = []
    for opt_name, color in color_map.items():
        legend_elements.append(plt.Line2D([0], [0], marker='o', color='w', 
                                         markerfacecolor=color, markersize=12,
                                         label=opt_name, markeredgecolor='white', 
                                         markeredgewidth=2))
    
    # Créer un boxplot pour chaque paramètre
    for i, param in enumerate(parameters):
        if i < len(axes):
            ax = axes[i]
            
            # Collecter les valeurs et couleurs
            param_values = []
            labels = []
            point_colors = []
            
            for pseudo_name, best_params in all_best_params.items():
                if param in best_params:
                    param_values.append(best_params[param])
                    labels.append(pseudo_name)
                    point_colors.append(color_map[pseudo_name])
            
            if param_values:
                # Style du boxplot minimaliste
                bp = ax.boxplot([param_values], 
                              patch_artist=True,
                              widths=0.4,
                              boxprops=dict(facecolor='lightblue', alpha=0.4, linewidth=1.5),
                              medianprops=dict(color='darkred', linewidth=3),
                              whiskerprops=dict(linewidth=1.5),
                              capprops=dict(linewidth=1.5))
                
                # Points individuels SANS annotations individuelles
                x_positions = []
                for j, (color, label) in enumerate(zip(point_colors, labels)):
                    # Jitter horizontal pour séparer les points
                    jitter = (j - len(point_colors)/2 + 0.5) * 0.05
                    x_pos = 1 + jitter
                    x_positions.append(x_pos)
                
                # Afficher les points avec jitter
                for j, (x, y, color, label) in enumerate(zip(x_positions, param_values, point_colors, labels)):
                    ax.scatter(x, y, c=[color], s=150, alpha=0.9, 
                             edgecolors='white', linewidth=2.5, zorder=3)
                
                # Informations du paramètre (titre et statistiques simplifiés)
                _style_parameter_subplot_clean(ax, param, param_values)
                
        else:
            # Cacher les axes non utilisés
            axes[i].set_visible(False)
    
    # Légende globale à droite de la figure
    fig.legend(handles=legend_elements, bbox_to_anchor=(0.98, 0.5), 
              loc='center left', frameon=True, fancybox=True, shadow=True,
              fontsize=12, title='Optimisations', title_fontsize=14)
    
    # Titre principal élégant
    fig.suptitle(f'Comparaison des Paramètres Optimisés\n'
                f'{len(optimizations_dict)} optimisations • '
                f'{len(parameters)} paramètres', 
                fontsize=18, fontweight='bold', y=0.96)
    
    # Ajustements esthétiques avec espace pour la légende
    plt.tight_layout(rect=[0, 0, 0.85, 0.92])
    
    # Sauvegarder
    if savefig:
        _save_figure(fig, filename, dateinname, "_individual")
    
    return fig


def _create_normalized_comparison_plot(all_best_params, parameters, optimizations_dict,
                                     savefig, filename, dateinname):
    """Crée le graphique avec paramètres normalisés - version épurée."""
    
    # Préparer les données pour normalisation
    data_for_normalization = []
    for param in parameters:
        param_values = []
        for pseudo_name, best_params in all_best_params.items():
            if param in best_params:
                param_values.append(best_params[param])
        if param_values:
            data_for_normalization.append(param_values)
    
    # Normaliser par paramètre (Min-Max scaling)
    normalized_data = {}
    for i, param in enumerate(parameters):
        if i < len(data_for_normalization):
            values = np.array(data_for_normalization[i])
            if len(values) > 1 and np.std(values) > 0:
                # Min-Max normalization
                normalized = (values - np.min(values)) / (np.max(values) - np.min(values))
            else:
                normalized = np.zeros_like(values)
            
            # Associer aux optimisations
            j = 0
            for pseudo_name, best_params in all_best_params.items():
                if param in best_params:
                    if pseudo_name not in normalized_data:
                        normalized_data[pseudo_name] = {}
                    normalized_data[pseudo_name][param] = normalized[j]
                    j += 1
    
    # Créer le DataFrame pour le plot
    plot_data = []
    for pseudo_name, params_dict in normalized_data.items():
        for param, normalized_value in params_dict.items():
            plot_data.append({
                'Optimization': pseudo_name,
                'Parameter': param,
                'Normalized_Value': normalized_value
            })
    
    df = pd.DataFrame(plot_data)
    
    # Créer la figure avec plus d'espace pour la légende
    fig, ax = plt.subplots(figsize=(max(14, len(parameters) * 1.4), 10))
    fig.patch.set_facecolor('white')
    
    # Boxplot minimaliste avec seaborn
    sns.boxplot(data=df, x='Parameter', y='Normalized_Value', 
                color='lightgray', ax=ax, width=0.4)
    
    # Styliser les boxplots pour les rendre plus discrets
    for patch in ax.patches:
        patch.set_alpha(0.2)
        patch.set_linewidth(1)
    
    # Couleurs distinctes pour chaque optimisation
    colors = sns.color_palette("Set2", len(optimizations_dict))  # Palette plus douce
    color_map = dict(zip(optimizations_dict.keys(), colors))
    
    # Points plus grands et plus visibles avec jitter horizontal pour éviter superposition
    np.random.seed(42)  # Pour reproductibilité
    for pseudo_name in optimizations_dict.keys():
        subset = df[df['Optimization'] == pseudo_name]
        if not subset.empty:
            # Ajouter un léger jitter horizontal pour séparer les points
            x_positions = []
            for param in subset['Parameter']:
                param_idx = list(df['Parameter'].unique()).index(param)
                jitter = np.random.uniform(-0.15, 0.15)
                x_positions.append(param_idx + jitter)
            
            ax.scatter(x_positions, subset['Normalized_Value'],
                      color=color_map[pseudo_name], label=pseudo_name,
                      s=150, alpha=0.9, edgecolors='white', linewidth=2.5, zorder=4)
    
    # Style du graphique épuré
    ax.set_xlabel('Paramètres', fontsize=16, fontweight='bold', labelpad=15)
    ax.set_ylabel('Valeurs Normalisées [0-1]', fontsize=16, fontweight='bold', labelpad=15)
    ax.set_title('Variabilité des Paramètres Normalisés\nentre Optimisations', 
                 fontsize=18, fontweight='bold', pad=30)
    
    # Améliorer les labels des paramètres
    param_labels = []
    for param in df['Parameter'].unique():
        # Simplifier les noms de paramètres
        clean_param = param.replace('Phy+', '').replace('Det', 'D').replace('Bac', 'B')
        param_labels.append(clean_param)
    
    # Fixer les ticks avant de changer les labels
    ax.set_xticks(range(len(param_labels)))
    ax.set_xticklabels(param_labels, rotation=45, ha='right', fontsize=12)
    ax.tick_params(axis='y', labelsize=12)
    
    # Légende centralisée et élégante
    legend = ax.legend(bbox_to_anchor=(0.5, -0.15), loc='upper center', 
                      ncol=min(4, len(optimizations_dict)), frameon=True,
                      fancybox=True, shadow=True, fontsize=14,
                      title='Optimisations', title_fontsize=16)
    legend.get_title().set_fontweight('bold')
    
    # Améliorer le style des lignes de la légende (gérer différents types d'objets)
    for handle in legend.legend_handles:
        try:
            handle.set_markersize(12)
            handle.set_alpha(1.0)
        except AttributeError:
            # Pour les PathCollection objects
            handle.set_sizes([144])  # 12^2 = 144 pour la taille
            handle.set_alpha(1.0)
    
    # Grille très subtile
    ax.grid(True, alpha=0.2, linestyle=':', axis='y')
    ax.set_axisbelow(True)
    
    # Limites et style des axes
    ax.set_ylim(-0.05, 1.05)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(1.5)
    ax.spines['bottom'].set_linewidth(1.5)
    
    # Ajuster l'espacement pour la légende
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.2)  # Espace pour la légende
    
    # Ajouter une ligne horizontale à 0.5 pour référence
    ax.axhline(y=0.5, color='darkgray', linestyle='--', alpha=0.5, 
               linewidth=1, label='Milieu de la plage')
    
    # Sauvegarder
    if savefig:
        _save_figure(fig, filename, dateinname, "_normalized")
    
    return fig


def _style_parameter_subplot_clean(ax, param, param_values):
    """Style un subplot de paramètre individuel - version épurée."""
    try:
        param_info = varinfos.ref_values[param]
        complete_name = param_info['complete_name']
        symbol = param_info.get('symbol', '')
        reference_value = param_info['reference_value']
        units = param_info.get('units', '')
        
        # Titre avec le symbole ou le nom complet
        title = fns.cleantext(symbol) if symbol else complete_name
        ax.set_title(title, fontsize=14, fontweight='bold', pad=20)
        
        # Ligne de référence discrète
        ax.axhline(y=reference_value, color='gray', 
                 linestyle=':', alpha=0.6, linewidth=1.5)
        
        # Label Y avec unités
        ylabel = param.replace('+', ' ')
        if units:
            ylabel += f"\n[{fns.cleantext(units)}]"
        ax.set_ylabel(ylabel, fontsize=12, fontweight='bold')
        
    except (KeyError, AttributeError):
        # Si pas d'info dans varinfos
        clean_param = param.replace('+', ' ').replace('Phy', 'Phyto').replace('Det', 'Detritus')
        ax.set_title(clean_param, fontsize=14, fontweight='bold', pad=20)
        ax.set_ylabel(clean_param, fontsize=12, fontweight='bold')
    
    # Statistiques simplifiées et plus lisibles
    mean_val = np.mean(param_values)
    std_val = np.std(param_values)
    cv = std_val / mean_val * 100 if mean_val != 0 else 0
    
    # Positionnement intelligent des statistiques
    stats_text = f'μ={mean_val:.3g} | σ={std_val:.3g}\nCV={cv:.1f}%'
    ax.text(0.02, 0.98, stats_text, 
           transform=ax.transAxes, va='top', ha='left',
           bbox=dict(boxstyle='round,pad=0.5', 
                   facecolor='white', alpha=0.9,
                   edgecolor='gray', linewidth=1),
           fontsize=11, fontweight='bold')
    
    # Style moderne des axes
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(2)
    ax.spines['bottom'].set_linewidth(2)
    
    # Supprimer les ticks X
    ax.set_xticks([])
    ax.grid(True, alpha=0.2, axis='y', linestyle='-')
    ax.set_axisbelow(True)
    
    # Améliorer les ticks Y
    ax.tick_params(axis='y', labelsize=10, colors='black', width=1.5)


def _style_parameter_subplot(ax, param, param_values):
    """Style un subplot de paramètre individuel - version legacy."""
    # Appeler la nouvelle fonction pour maintenir la compatibilité
    _style_parameter_subplot_clean(ax, param, param_values)


def _save_figure(fig, filename, dateinname, suffix):
    """Sauvegarde une figure avec le bon nommage."""
    datestr = datetime.now().strftime('%Y%m%d_') if dateinname else ''
    
    if filename is None:
        filename = "parameters_comparison"
    
    fig_path = os.path.join(FIGURE_PATH, f'{datestr}{filename}{suffix}.png')
    fig.savefig(fig_path, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"Figure sauvegardée: {fig_path}")


def get_optimization_summary(optimizations_dict: Dict[str, str]) -> pd.DataFrame:
    """
    Génère un DataFrame résumé avec les informations principales de chaque optimisation.
    
    Args:
        optimizations_dict: Dictionnaire avec pseudo-nom -> nom officiel
        
    Returns:
        DataFrame avec les informations de résumé
    """
    
    summary_data = []
    
    for pseudo_name, opt_name in optimizations_dict.items():
        try:
            opt = optim.Optimization.load_existing(opt_name)
            
            if not hasattr(opt, 'summary'):
                opt.process_results()
            
            # Informations basiques
            row = {
                'Pseudo_Name': pseudo_name,
                'Official_Name': opt_name,
                'Best_Score': opt.summary['best_score'],
                'Generations_Completed': opt.summary['convergence']['generations_completed'],
                'Last_Improvement': opt.summary['convergence']['last_improvement'],
                'Num_Parameters': len(opt.summary['best_parameters'])
            }
            
            # Ajouter les meilleurs paramètres
            for param, value in opt.summary['best_parameters'].items():
                row[f'Best_{param}'] = value
            
            summary_data.append(row)
            
        except Exception as e:
            print(f"Erreur lors du traitement de {pseudo_name}: {e}")
            continue
    
    return pd.DataFrame(summary_data)


if __name__ == "__main__":
    # Exemple d'utilisation
    example_optimizations = {
        'Test1': 'OPT001',
        'Test2': 'OPT002',
        'Test3': 'OPT003'
    }
    
    print("Exemple d'utilisation:")
    print("optimizations = {'Test1': 'OPT001', 'Test2': 'OPT002', 'Test3': 'OPT003'}")
    print("fig = compare_optimizations_parameters(optimizations)")
    print("plt.show()")
    print("\n# Pour obtenir un résumé tabulaire:")
    print("summary = get_optimization_summary(optimizations)")
    print("print(summary)")