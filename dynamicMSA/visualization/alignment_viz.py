"""
[EN]
It visualizes the final alignment as colored blocks.

Each nucleotide/amino acid is shown with its own color:
For DNA: A=green, T=red, G=yellow, C=blue, -=gray
Gap columns are clearly visible. Matching positions in all sequences are highlighted.

[TR]
Final hizalamayı renkli bloklar halinde görselleştirir.

Her nükleotid/amino asit kendi rengiyle gösterilir:
  DNA için: A=yeşil, T=kırmızı, G=sarı, C=mavi, -=gri
  Gap sütunları açıkça belli olur.
  Tüm dizilerde eşleşen pozisyonlar vurgulanır.
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np


# DNA/RNA nükleotid renkleri (biyoinformatikte standart renkler)
NUCLEOTIDE_COLORS = {
    'A': '#2ecc71',   # yeşil
    'T': '#e74c3c',   # kırmızı
    'U': '#e74c3c',   # kırmızı (RNA)
    'G': '#f39c12',   # turuncu-sarı
    'C': '#3498db',   # mavi
    '-': '#bdc3c7',   # gri (gap)
    'N': '#9b59b6',   # mor (belirsiz)
}

DEFAULT_COLOR = '#95a5a6'


def plot_alignment(
    aligned: list[str],
    names: list[str] = None,
    title: str = "Çoklu Dizi Hizalaması",
    highlight_conserved: bool = True,
    save_path: str = None,
    show: bool = True
) -> plt.Figure:
    """
    [EN]
    It plots aligned arrays as a colored table.

    Parameters
    ------------
    aligned : list[str] — aligned arrays (of equal length)
    names : list[str] — array names
    title : str — chart title
    highlight_conserved: bool — highlight positions that are the same in all arrays
    save_path : str — save path (optional)
    show : bool

    Returns
    --------
    plt.Figure

    [TR]
    Hizalanmış dizileri renkli tablo olarak çizer.

    Parametreler
    ------------
    aligned            : list[str]  — hizalanmış diziler (eşit uzunlukta)
    names              : list[str]  — dizi isimleri
    title              : str        — grafik başlığı
    highlight_conserved: bool       — tüm dizilerde aynı olan pozisyonları vurgula
    save_path          : str        — kayıt yolu (opsiyonel)
    show               : bool

    Döndürür
    --------
    plt.Figure
    """
    if not aligned:
        raise ValueError("Hizalanmış dizi listesi boş.")

    n_seqs = len(aligned)
    length = len(aligned[0])

    if names is None:
        names = [f"Seq{i}" for i in range(n_seqs)]

    # Korunmuş pozisyonları belirle (tüm dizilerde aynı, gap değil)
    conserved = set()
    if highlight_conserved:
        for pos in range(length):
            col = [seq[pos] for seq in aligned]
            if len(set(col)) == 1 and col[0] != '-':
                conserved.add(pos)

    # Figür boyutu: uzunluğa ve dizi sayısına göre otomatik
    fig_width  = max(12, length * 0.45)
    fig_height = max(4, n_seqs * 0.8 + 2)

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    ax.set_xlim(0, length)
    ax.set_ylim(0, n_seqs)
    ax.axis('off')

    for seq_idx, (seq, name) in enumerate(zip(aligned, names)):
        y = n_seqs - seq_idx - 1  # Yukarıdan aşağıya sıralama

        # Dizi adını sol tarafa yaz
        ax.text(-0.5, y + 0.5, name, ha='right', va='center',
                fontsize=10, fontweight='bold', fontfamily='monospace')

        for pos, char in enumerate(seq):
            color = NUCLEOTIDE_COLORS.get(char.upper(), DEFAULT_COLOR)

            # Korunmuş pozisyonları çerçevele
            if pos in conserved:
                border_color = '#2c3e50'
                border_width = 2.0
                alpha = 1.0
            else:
                border_color = 'white'
                border_width = 0.5
                alpha = 0.85

            rect = mpatches.FancyBboxPatch(
                (pos + 0.05, y + 0.1),
                0.9, 0.8,
                boxstyle="round,pad=0.05",
                facecolor=color,
                edgecolor=border_color,
                linewidth=border_width,
                alpha=alpha
            )
            ax.add_patch(rect)

            # Karakter etiketi
            ax.text(pos + 0.5, y + 0.5, char,
                    ha='center', va='center',
                    fontsize=10, fontweight='bold',
                    fontfamily='monospace',
                    color='white' if char != '-' else '#7f8c8d')

    # Pozisyon numaraları (her 5'te bir)
    for pos in range(0, length, max(1, length // 10)):
        ax.text(pos + 0.5, n_seqs + 0.1, str(pos + 1),
                ha='center', va='bottom', fontsize=7, color='#7f8c8d')

    # Renk açıklamaları
    legend_handles = [
        mpatches.Patch(facecolor=c, label=nuc, edgecolor='white')
        for nuc, c in NUCLEOTIDE_COLORS.items() if nuc != 'N'
    ]
    if highlight_conserved and conserved:
        legend_handles.append(
            mpatches.Patch(facecolor='white', edgecolor='#2c3e50',
                           linewidth=2, label='Korunmuş pozisyon')
        )

    ax.legend(handles=legend_handles, loc='lower right',
              bbox_to_anchor=(1.0, -0.05), ncol=len(legend_handles),
              fontsize=9, framealpha=0.9)

    ax.set_title(title, fontsize=14, fontweight='bold', pad=20)
    fig.patch.set_facecolor('white')
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight',
                    facecolor='white', edgecolor='none')

    if show:
        plt.show()

    return fig