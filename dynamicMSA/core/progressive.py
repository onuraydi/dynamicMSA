"""
[EN]
Progressive alignment — the heart of the project.

Working logic:
1. Align the two closest arrays first according to the guide tree.
2. Treat the resulting alignment as a "profile".
3. Add the next array to this profile.
4. Continue until all arrays are added.

What is a profile?
A representation of multiple aligned arrays together.
For example, after two arrays are aligned:
AC-GT
ACGGT
These two rows are now a "profile". When adding the third array to it,
we move these two rows together.

[TR]
Aşamalı hizalama (Progressive Alignment) — projenin kalbi.
 
Çalışma mantığı:
  1. Kılavuz ağaca göre önce en yakın iki diziyi hizala
  2. Elde edilen hizalamayı bir "profil" olarak ele al
  3. Bir sonraki diziyi bu profile ekle
  4. Tüm diziler eklenene kadar devam et
 
Profil nedir?
  Birden fazla hizalanmış dizinin birlikte temsil edilmesi.
  Örneğin iki dizi hizalandıktan sonra:
    AC-GT
    ACGGT
  Bu iki satır artık bir "profil"dir. Üçüncü diziyi buna eklerken
  bu iki satırı birlikte hareket ettiririz.
"""

from dynamicMSA.core.needleman_wunsch import needleman_wunsch
from dynamicMSA.core.distance import build_distance_matrix
from dynamicMSA.core.guide_tree import upgma, get_merge_order
from dynamicMSA.utils.scoring import ScoringMatrix

def align_profile_to_sequence(
    profile: list[str],
    sequence: str,
    scoring: ScoringMatrix
) -> list[str]:
    """
    [EN]
    Aligns a profile (multiple aligned arrays) into a single array.

    Method:
    Extract the "consensus" array of the profile (most frequent character),
    align it with the new array NW,
    then apply gap positions to all profile rows.

    Parameters
    ------------
    profile : list[str] — an aligned array where each element is of equal length
    sequence : str — a new array to be added to the profile
    scoring : ScoringMatrix

    Returns
    --------
    list[str] : expanded profile (new array added)

    [TR]
    Bir profili (birden fazla hizalanmış dizi) tek bir diziye hizalar.
 
    Yöntem:
      Profilin "konsensüs" dizisini çıkar (en sık görülen karakter),
      bunu yeni dizi ile NW ile hizala,
      sonra gap pozisyonlarını tüm profil satırlarına uygula.
 
    Parametreler
    ------------
    profile  : list[str] — her elemanı eşit uzunlukta hizalanmış bir dizi
    sequence : str       — profile eklenecek yeni dizi
    scoring  : ScoringMatrix
 
    Döndürür
    --------
    list[str] : genişletilmiş profil (yeni dizi de eklendi)
    """

    consensus = _profile_consensus(profile)
    aligned_consensus, aligned_seq, _,_,_ = needleman_wunsch(consensus, sequence, scoring)

    # Konsensüs hizalamasındaki gap pozisyonlarını profile yansıt
    updated_profile = _apply_gaps_to_profile(profile, aligned_consensus, consensus)
    updated_profile.append(aligned_seq)

    return updated_profile

def _profile_consensus(profile: list[str]) -> str:
    """
    [EN]
    It selects the most frequently used character in each position of the profile. 
    The '-' character (or else 'gap') is not included in the character count.
    
    [TR]
    Profilin her pozisyonundaki en sık karakteri seçer.
    Gap '-' karakter sayımına dahil edilmez (yoksa hep gap kazanır).
    profillerden temsilci ve tek bir konsensüs dizisi oluşturur.
    Dizi1:  A - C T G
    Dizi2:  A - C G G
    Dizi3:  T - C G C
    temsilci: A - C G G  (her pozisyonda en sık karakter)
    """
    if not profile:
        return ""
    
    length = len(profile[0])
    consensus = []

    for pos in range(length):
        counts = {}
        for seq in profile:
            ch = seq[pos]
            if ch != '-':
                counts[ch] = counts.get(ch, 0) + 1

        if counts:
            consensus.append(max(counts, key=counts.get))
        else:
            consensus.append('-')

    return ''.join(consensus)


def _apply_gaps_to_profile(
    profile: list[str],
    aligned_consensus: str,
    original_consensus: str
) -> list[str]:
    """
    [EN]
    New gaps added to the consensus during NW alignment
    are added to all lines of the profile.

    Example:
    original_consensus: ACGT
    aligned_consensus: AC-GT ← Gap added to position 3
    → Add '-' to position 3 of all lines of the profile
    
    [TR]
    NW hizalaması sırasında konsensüse eklenen yeni gap'leri
    profilin tüm satırlarına ekler.
 
    Örnek:
      original_consensus : ACGT
      aligned_consensus  : AC-GT   ← 3. pozisyona gap eklendi
      → Profilin tüm satırlarının 3. pozisyonuna '-' ekle
    """

    # Orijinal konsensüste olmayan gap'lerin pozisyonlarını bul

    new_gaps = []
    orig_idx = 0

    for i, ch in enumerate(aligned_consensus):
        if ch == '-':
            new_gaps.append(i)
        else:
            orig_idx += 1

    # Her profil satırına yeni gapları ekle
    updated = []
    for seq in profile:
        seq_list = list(seq)
        for gap_pos in new_gaps:
            seq_list.insert(gap_pos, '-')
        updated.append(''.join(seq_list))

    return updated

def progressive_align(
    sequences: list[str],
    names: list[str] = None,
    scoring: ScoringMatrix = None
) -> dict:
    """
    [EN]
    Main function: Generates MSA with Progressive Alignment.

    Steps:
    1. Calculate distance matrix
    2. Create guide tree with UPGMA
    3. Merge profiles in the order shown by the tree

    Parameters
    ------------
    sequences : list[str] — raw sequences to align
    names : list[str] — sequence names (optional)
    scoring : ScoringMatrix

    Returns
    --------
    dict with the following keys:
    "aligned" : list[str] — aligned sequences
    "names" : list[str] — sequence names
    "distance_matrix" : list[list[float]] — distance matrix
    "guide_tree" : dict — UPGMA tree
    "score" : float — total alignment score

    [TR]
    Ana fonksiyon: Progressive Alignment ile MSA üretir.
 
    Adımlar:
      1. Mesafe matrisi hesapla
      2. UPGMA ile kılavuz ağaç oluştur
      3. Ağacın gösterdiği sırayla profilleri birleştir
 
    Parametreler
    ------------
    sequences : list[str]    — hizalanacak ham diziler
    names     : list[str]    — dizi isimleri (opsiyonel)
    scoring   : ScoringMatrix
 
    Döndürür
    --------
    dict şu anahtarlarla:
      "aligned"          : list[str]           — hizalanmış diziler
      "names"            : list[str]           — dizi isimleri
      "distance_matrix"  : list[list[float]]   — mesafe matrisi
      "guide_tree"       : dict                — UPGMA ağacı
      "score"            : float               — toplam hizalama skoru
    """
    
    if scoring is None:
        scoring = ScoringMatrix()

    if names is None:
        names = [f"seq{i}" for i in range(len(sequences))]

    n = len(sequences)

    if n == 1:
        return {
            "aligned": sequences[:],
            "names": names,
            "distance_matrix": [[0.0]],
            "guide_tree": {"index": 0, "name": names[0], "sequences": [0]},
            "score": 0.0
        }
    
    # ---------------------------------------------------------------
    # ADIM 1: Mesafe matrisi
    # ---------------------------------------------------------------
    
    dist_matrix = build_distance_matrix(sequences, scoring)
 
    # ---------------------------------------------------------------
    # ADIM 2: Kılavuz ağaç
    # ---------------------------------------------------------------
    
    tree = upgma(dist_matrix, names)

    # ---------------------------------------------------------------
    # ADIM 3: Aşamalı hizalama
    # profile_map: dizi indeksi → o ana kadar hizalanmış satırlar
    # ---------------------------------------------------------------

    # Her dizi kendi profilini oluşturur: {indeks: [dizi]}
    profile_map = {i: [sequences[i]] for i in range(n)}

    # UPGMA ağacını post-order gezerek birleşme sırasını al
    # Her adımda iki alt grubu birleştirip yeni bir profil oluşturuyoruz
    _merge_profiles(tree, profile_map, scoring)

    # Kök düğümün sequences listesi final sırayı verir
    final_order = tree["sequences"]
    root_key    = _get_root_key(tree)
    aligned     = profile_map[root_key]   # kök profili

    aligned_ordered = aligned          # zaten doğru sırada
    names_ordered   = [names[i] for i in final_order]

    # Toplam skor: tüm kolon çiftlerinin ortalama skoru
    total_score = _compute_alignment_score(aligned_ordered, scoring)

    return {
        "aligned"         : aligned_ordered,
        "names"           : names_ordered,
        "distance_matrix" : dist_matrix,
        "guide_tree"      : tree,
        "score"           : total_score
    }

def _merge_profiles(node: dict, profile_map: dict, scoring) -> None:
    """
    UPGMA ağacını post-order (sol → sağ → kök) gezerek profilleri birleştirir.

    Yaprak düğüm    : profile_map'te zaten var, bir şey yapma
    İç düğüm (birleşme) : sol ve sağ çocukları önce işle, sonra birleştir
    Sonuç         : profile_map[root_key] tüm dizileri içerir
    """
    # Yaprak düğüm: tek dizi, zaten profile_map'te
    if "index" in node:
        return

    # Önce sol ve sağ alt ağaçları işle (post-order)
    _merge_profiles(node["left"],  profile_map, scoring)
    _merge_profiles(node["right"], profile_map, scoring)

    # Sol ve sağ profillerin başlangıç anahtarını bul
    left_key  = _get_root_key(node["left"])
    right_key = _get_root_key(node["right"])

    left_profile  = profile_map[left_key]
    right_profile = profile_map[right_key]

    # Sol profili, sağ profildeki tüm dizilerle hizala
    combined = left_profile
    for seq in right_profile:
        combined = align_profile_to_sequence(combined, seq.replace('-', ''), scoring)

    # Birleşik profili sol alt ağacın ilk yaprağının anahtarıyla kaydet
    profile_map[left_key] = combined


def _get_root_key(node: dict) -> int:
    """
    Bir düğümün profile_map'teki anahtarını döndürür.
    Yaprak için: node["index"]
    İç düğüm için: sol alt ağacın ilk yaprağının indeksi
    """
    if "index" in node:
        return node["index"]
    return _get_root_key(node["left"])


def _compute_alignment_score(aligned: list[str], scoring: ScoringMatrix) -> float:
    """
    [EN]
    The total score of aligned sequences: the sum of all column pairs. All sequence pairs are evaluated for each column (sum-of-pairs score).
    
    [TR]
    Hizalanmış dizilerin toplam skoru: tüm kolon çiftlerinin toplamı.
    Her kolon için tüm dizi çiftleri değerlendirilir (sum-of-pairs skoru).
    """
    if not aligned or len(aligned[0]) == 0:
        return 0.0
    
    total = 0.0
    length =len(aligned[0])

    for pos in range(length):
        col = [seq[pos] for seq in aligned]
        for i in range(len(col)):
            for j in range(i + 1, len(col)):
                total += scoring.score(col[i], col[j])

    return total