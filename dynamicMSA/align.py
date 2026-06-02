"""
[EN]
The main interface that the user calls directly.

Usage:
    from dynamicMSA import align

    result = align(
    
    sequences=["ACGT", "ACGGT", "ACGAT"],
    names=["Seq1", "Seq2", "Seq3"],
    visualize=True
    )
    
    print(result["aligned"])
    print(result["score"])

[TR]
Kullanıcının doğrudan çağırdığı ana arayüz.

Kullanım:
    from dynamicMSA import align

    result = align(
        sequences=["ACGT", "ACGGT", "ACGAT"],
        names=["Seq1", "Seq2", "Seq3"],
        visualize=True
    )

    print(result["aligned"])
    print(result["score"])
"""

from dynamicMSA.core.progressive import progressive_align
from dynamicMSA.utils.scoring import ScoringMatrix


def align(
    sequences: list[str],
    names: list[str] = None,
    match: int = 1,
    mismatch: int = -1,
    gap: int = -2,
    visualize: bool = True,
    save_figures: bool = False,
    output_dir: str = ".",
    verbose: bool = True # konsol çıktıları için
) -> dict:
    """
    [EN]
    It performs multi-array alignment and visualizes it optionally.

    Parameters
    ------------
    sequences : list[str] — sequences to align (at least 2)
    names : list[str] — sequence names (optional)
    match : int — match score (default: +1)
    mismatch : int — mismatch penalty (default: -1)
    gap : int — gap penalty (default: -2)
    visualize : bool — show matplotlib graphs (default: True)
    save_figures : bool — save graphs to a file (default: False)
    output_dir : str — save directory (default: ".")
    verbose : bool — progress messages (default: True)

    Returns
    --------
    dict with these keys:
    "aligned" : list[str] — aligned sequences
    "names" : list[str] — sequence names
    "score" : float — total alignment score
    "distance_matrix" : list[list[float]] — distances between array pairs
    "guide_tree" : dict — UPGMA tree
    
    [TR]
    Çoklu dizi hizalaması yapar ve isteğe bağlı olarak görselleştirir.

    Parametreler
    ------------
    sequences    : list[str]  — hizalanacak diziler (en az 2)
    names        : list[str]  — dizi isimleri (opsiyonel)
    match        : int        — eşleşme skoru (varsayılan: +1)
    mismatch     : int        — uyumsuzluk cezası (varsayılan: -1)
    gap          : int        — boşluk cezası (varsayılan: -2)
    visualize    : bool       — matplotlib grafiklerini göster (varsayılan: True)
    save_figures : bool       — grafikleri dosyaya kaydet (varsayılan: False)
    output_dir   : str        — kayıt dizini (varsayılan: ".")
    verbose      : bool       — ilerleme mesajları (varsayılan: True)

    Döndürür
    --------
    dict şu anahtarlarla:
      "aligned"         : list[str]         — hizalanmış diziler
      "names"           : list[str]         — dizi isimleri
      "score"           : float             — toplam hizalama skoru
      "distance_matrix" : list[list[float]] — dizi çiftleri arası mesafeler
      "guide_tree"      : dict              — UPGMA ağacı
    """
    # ---------------------------------------------------------------
    # Girdi doğrulama
    # ---------------------------------------------------------------
    if not sequences or len(sequences) < 2:
        raise ValueError("En az 2 dizi gereklidir.")

    for i, seq in enumerate(sequences):
        if not seq:
            raise ValueError(f"Dizi {i} boş olamaz.")
        if not all(c.upper() in 'ACGTUN-' for c in seq):
            # Protein dizisi de olabilir, sadece uyarı ver
            if verbose:
                print(f"⚠ Uyarı: Dizi {i} standart dışı karakterler içeriyor. "
                      f"Protein dizisi mi?")

    sequences = [seq.upper() for seq in sequences]

    if names is None:
        names = [f"Seq{i+1}" for i in range(len(sequences))]

    original_names = list(names)

    scoring = ScoringMatrix(match=match, mismatch=mismatch, gap=gap)

    # ---------------------------------------------------------------
    # ADIM 1: Hizalama
    # ---------------------------------------------------------------
    if verbose:
        print(f"🔬 dynamicMSA başlatılıyor...")
        print(f"   {len(sequences)} dizi | Skor: match={match}, "
              f"mismatch={mismatch}, gap={gap}")
        print("   Mesafe matrisi hesaplanıyor...")

    result = progressive_align(sequences, names, scoring)

    if verbose:
        _pretty_print(result, original_names)

    # ---------------------------------------------------------------
    # ADIM 2: Görselleştirme
    # ---------------------------------------------------------------
    if visualize:
        _visualize(sequences, names, result, save_figures, output_dir, verbose)

    return result


def _pretty_print(result: dict, original_names: list[str]) -> None:
    """Hizalama sonucunu düzenli formatta konsola yazdırır."""
    aligned = result["aligned"]
    names   = result["names"] # Hizalanmış (yeni) isim sırası
    dm      = result["distance_matrix"]
    history = result["guide_tree"].get("merge_history", [])

    # Sütun genişliklerini orijinal isimlere göre de kontrol et
    col_w = max(max(len(n) for n in names), max(len(n) for n in original_names)) + 2
    seq_w = len(aligned[0]) if aligned else 0
    inner_w = 2 + col_w + seq_w + 2

    print()
    print("╔" + "═" * inner_w + "╗")
    print("║" + " HIZALAMA SONUCU ".center(inner_w) + "║")
    print("╠" + "═" * inner_w + "╣")
    
    # Hizalama sonucunu yazdırırken YENİ sırayı kullanıyoruz
    for name, seq in zip(names, aligned):
        print(f"║  {name:<{col_w}}{seq}  ║")
    
    print("╠" + "═" * inner_w + "╣")
    score_text = f"  Toplam Skor : {result['score']}"
    print(f"║{score_text:<{inner_w}}║")
    print("╚" + "═" * inner_w + "╝")

    print()
    print("┌─ MESAFE MATRİSİ ")
    
    # DİKKAT: Matrisi yazdırırken ORİJİNAL isim sırasını kullanıyoruz!
    header = " " * (col_w + 2) + "".join(f"{n:^10}" for n in original_names)
    print("│ " + header)
    for i, name in enumerate(original_names):
        row = "".join(f"{dm[i][j]:^10.3f}" for j in range(len(original_names)))
        print(f"│  {name:<{col_w}}{row}")
    
    print("└" + "─" * (len(header) + 1))

    if history:
        print()
        print("┌─ BİRLEŞME SIRASI (UPGMA)")
        for step, h in enumerate(history, 1):
            a, b = h["merged"]
            print(f"│  Adım {step}: '{a}' + '{b}'  →  mesafe: {h['distance']:.4f}")
        print("└" + "─" * 40)


def _visualize(sequences, names, result, save_figures, output_dir, verbose):
    """Tüm görselleştirmeleri sırayla çağırır."""
    import os
    from dynamicMSA.core.needleman_wunsch import needleman_wunsch
    from dynamicMSA.utils.scoring import ScoringMatrix
    from dynamicMSA.visualization.dp_matrix import plot_dp_matrix
    from dynamicMSA.visualization.guide_tree_viz import plot_guide_tree
    from dynamicMSA.visualization.alignment_viz import plot_alignment

    scoring = ScoringMatrix()

    if verbose:
        print("\n📊 Görselleştirmeler oluşturuluyor...")

    # 1. DP Matrisi: tüm dizi çiftleri için göster
    for i in range(len(sequences)):
        for j in range(i + 1, len(sequences)):
            _, _, _, matrix, path = needleman_wunsch(sequences[i], sequences[j], scoring)
            plot_dp_matrix(
                matrix=matrix,
                seq1=sequences[i],
                seq2=sequences[j],
                path=path,
                title=f"DP Matrisi: {names[i]} vs {names[j]}",
                save_path=os.path.join(output_dir, f"dp_matrix_{i}_{j}.png") if save_figures else None,
                show=True
            )

    # 2. Kılavuz Ağaç
    plot_guide_tree(
        sequences=sequences,
        names=names,
        title="UPGMA Kılavuz Ağacı",
        save_path=os.path.join(output_dir, "guide_tree.png") if save_figures else None,
        show=True
    )

    # 3. Final Hizalama
    plot_alignment(
        aligned=result["aligned"],
        names=result["names"],
        title="Çoklu Dizi Hizalaması (Final)",
        highlight_conserved=True,
        save_path=os.path.join(output_dir, "alignment.png") if save_figures else None,
        show=True
    )

    if verbose and save_figures:
        print(f"   Grafikler '{output_dir}' dizinine kaydedildi.")