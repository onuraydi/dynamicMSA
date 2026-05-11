"""
[EN]
It visualizes the Needleman-Wunsch DP matrix as a heat map. It also plots the traceback path on the matrix.

Why do we visualize it?
Seeing the DP matrix makes the algorithm work more concrete.
Seeing which cell we came from makes it easier to understand the traceback.

[TR]
Needleman-Wunsch'un DP matrisini ısı haritası olarak görselleştirir.
Traceback yolunu da matrisin üzerine çizer.

Neden görselleştiriyoruz?
  DP matrisini görmek, algoritmanın nasıl çalıştığını somut hale getirir.
  Hangi hücreden geldiğimizi görmek traceback'i anlamayı kolaylaştırır.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap


def plot_dp_matrix(
    matrix: np.ndarray,
    seq1: str,
    seq2: str,
    path: list = None,
    title: str = "Needleman-Wunsch DP Matrisi",
    save_path: str = None,
    show: bool = True
) -> plt.Figure:
    """
    [EN]
    Plots the DP matrix as a heat map.

    Parameters
    ------------
    matrix : np.ndarray — matrix populated from NW
    seq1, seq2: str — aligned arrays (for axis labels)
    path : list — traceback path [(i,j), ...] (optional)
    title : str — chart title
    save_path : str — path to save file (optional)
    show : bool — call plt.show()

    Returns
    --------
    plt.Figure

    [TR]
    DP matrisini ısı haritası olarak çizer.

    Parametreler
    ------------
    matrix    : np.ndarray  — NW'den gelen doldurulan matris
    seq1, seq2: str         — hizalanan diziler (eksen etiketleri için)
    path      : list        — traceback yolu [(i,j), ...] (opsiyonel)
    title     : str         — grafik başlığı
    save_path : str         — kaydedilecek dosya yolu (opsiyonel)
    show      : bool        — plt.show() çağrılsın mı

    Döndürür
    --------
    plt.Figure
    """
    fig, ax = plt.subplots(figsize=(max(8, len(seq2) * 0.7), max(6, len(seq1) * 0.7)))

    # Özel renk paleti: düşük skor = koyu mavi, yüksek skor = açık sarı
    cmap = LinearSegmentedColormap.from_list(
        "nw_cmap", ["#1a1a2e", "#16213e", "#0f3460", "#e94560", "#f5a623"]
    )

    im = ax.imshow(matrix, cmap=cmap, aspect='auto')
    plt.colorbar(im, ax=ax, label="Skor")

    # Eksen etiketleri
    x_labels = ['-'] + list(seq2)
    y_labels = ['-'] + list(seq1)
    ax.set_xticks(range(len(x_labels)))
    ax.set_yticks(range(len(y_labels)))
    ax.set_xticklabels(x_labels, fontsize=11, fontweight='bold')
    ax.set_yticklabels(y_labels, fontsize=11, fontweight='bold')
    ax.set_xlabel(f"Dizi 2: {seq2}", fontsize=12, labelpad=10)
    ax.set_ylabel(f"Dizi 1: {seq1}", fontsize=12, labelpad=10)

    # Hücre içine skor değerlerini yaz (küçük matrisler için)
    if matrix.shape[0] <= 15 and matrix.shape[1] <= 15:
        for i in range(matrix.shape[0]):
            for j in range(matrix.shape[1]):
                val = matrix[i, j]
                color = 'white' if val < matrix.max() * 0.6 else 'black'
                ax.text(j, i, f"{val:.0f}", ha='center', va='center',
                        fontsize=9, color=color, fontweight='bold')

    # Traceback yolunu çiz
    if path:
        path_i = [p[0] for p in path]
        path_j = [p[1] for p in path]
        ax.plot(path_j, path_i, 'o-', color='#00ff88', linewidth=2.5,
                markersize=8, markerfacecolor='white', markeredgecolor='#00ff88',
                markeredgewidth=2, label='Traceback yolu', zorder=5)

        # Başlangıç ve bitiş noktalarını işaretle
        ax.plot(path_j[0], path_i[0], 's', color='#00ff88', markersize=12,
                label='Başlangıç', zorder=6)
        ax.plot(path_j[-1], path_i[-1], '*', color='#ff6b6b', markersize=14,
                label='Bitiş', zorder=6)

        ax.legend(loc='upper left', fontsize=9, framealpha=0.8)

    ax.set_title(title, fontsize=14, fontweight='bold', pad=15)
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight',
                    facecolor='white', edgecolor='none')

    if show:
        plt.show()

    return fig