"""
[EN]
UPGMA draws the generated guide tree as a dendrogram.

What is a dendrogram?
A tree diagram showing hierarchical clustering.
The y-axis shows the merging height (distance).
Sections merging earlier (lower) are more similar.

[TR]
UPGMA'nın ürettiği kılavuz ağacı dendrogram olarak çizer.

Dendrogram nedir?
  Hiyerarşik kümelemeyi gösteren ağaç diyagramı.
  Y ekseni birleşme yüksekliğini (mesafeyi) gösterir.
  Daha erken (altta) birleşen diziler daha benzerdir.
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from dynamicMSA.core.distance import build_distance_matrix
from dynamicMSA.core.guide_tree import upgma
from scipy.cluster.hierarchy import dendrogram, linkage


def plot_guide_tree(
    sequences: list[str],
    names: list[str] = None,
    title: str = "UPGMA Kılavuz Ağacı",
    save_path: str = None,
    show: bool = True
) -> plt.Figure:
    """
    [EN]
    Draws a guide tree using Scipy's dendrogram function.

    Parameters
    ------------
    sequences : list[str] — raw sequences
    names : list[str] — sequence names
    title : str — graph title
    save_path : str — save path (optional)
    show : bool — call plt.show()

    Returns
    --------
    plt.Figure
    
    [TR]
    Scipy'nin dendrogram fonksiyonunu kullanarak kılavuz ağacı çizer.

    Parametreler
    ------------
    sequences : list[str]  — ham diziler
    names     : list[str]  — dizi isimleri
    title     : str        — grafik başlığı
    save_path : str        — kayıt yolu (opsiyonel)
    show      : bool       — plt.show() çağrılsın mı

    Döndürür
    --------
    plt.Figure
    """
    if names is None:
        names = [f"Seq{i}" for i in range(len(sequences))]

    # Mesafe matrisini scipy linkage formatına çevir
    dist_matrix = build_distance_matrix(sequences)
    n = len(sequences)

    # Üçgen matrisi düzleştir (condensed format)
    condensed = []
    for i in range(n):
        for j in range(i + 1, n):
            condensed.append(dist_matrix[i][j])

    # UPGMA = average linkage
    Z = linkage(condensed, method='average')

    fig, ax = plt.subplots(figsize=(10, 5))

    dendrogram(
        Z,
        labels=names,
        ax=ax,
        color_threshold=0.7 * max(Z[:, 2]) if len(Z) > 0 else 0,
        leaf_font_size=12,
        leaf_rotation=30,
        above_threshold_color='#4a4a8a',
    )

    ax.set_title(title, fontsize=14, fontweight='bold', pad=15)
    ax.set_xlabel("Diziler", fontsize=12)
    ax.set_ylabel("Mesafe (Evrimsel Uzaklık)", fontsize=12)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Arka plan rengi
    ax.set_facecolor('#f8f9fa')
    fig.patch.set_facecolor('white')

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight',
                    facecolor='white', edgecolor='none')

    if show:
        plt.show()

    return fig