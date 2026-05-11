"""
[EN]
Calculates the evolutionary distance between two sequences.

Why do we need distance?
Needleman-Wunsch gives us a "score" — high score = similar sequences.
But for the UPGMA algorithm, we need to convert the score to a "distance".
Distance = how DIFFERENT they are.

Method used: Kimura 2-Parameter (simplified)
p = mismatched position rate in alignment
distance = -ln(1 - p) [simplified version]

[TR]
İki dizi arasındaki evrimsel mesafeyi hesaplar.
 
Neden mesafeye ihtiyacımız var?
  Needleman-Wunsch bize bir "skor" verir — yüksek skor = benzer diziler.
  Ama UPGMA algoritması için skoru bir "mesafe"ye (uzaklık) dönüştürmemiz gerekir.
  Mesafe = ne kadar FARKLI oldukları.
 
Kullanılan yöntem: Kimura 2-Parameter (basitleştirilmiş)
  p = hizalamada eşleşmeyen pozisyon oranı
  distance = -ln(1 - p)   [basit versiyon]
"""

import math
from dynamicMSA.core.needleman_wunsch import needleman_wunsch
from dynamicMSA.utils.scoring import ScoringMatrix


def percent_identity(aligned1: str, aligned2: str) -> float:
    """
    [EN]
    Calculates the percentage of identity between two aligned arrays.

    Identity = (number of matching positions) / (total number of positions)

    Parameters
    ------------
    aligned1, aligned2 : str
    Aligned arrays of the same length, which may contain gaps.

    Returns
    --------
    float : a value between 0.0 (no match) and 1.0 (exact match)

    [TR]
    İki hizalanmış dizi arasındaki özdeşlik yüzdesini hesaplar.
 
    Özdeşlik = (eşleşen pozisyon sayısı) / (toplam pozisyon sayısı)
 
    Parametreler
    ------------
    aligned1, aligned2 : str
        Aynı uzunlukta, gap içerebilen hizalanmış diziler.
 
    Döndürür
    --------
    float : 0.0 (hiç eşleşme yok) ile 1.0 (tam eşleşme) arasında değer
    """
    if len(aligned1) != len(aligned2):
        raise ValueError("Hizalanmış diziler aynı uzunlukta olmalı.")
    
    matches = sum(a == b for a, b in zip(aligned1, aligned2) if a != '-' and b != '-')
    total   = sum(1 for a, b in zip(aligned1, aligned2) if a != '-' and b != '-')

    return matches / total if total > 0 else 0.0


def compute_distance(seq1: str, seq2: str, scoring: ScoringMatrix = None) -> float:
    """
    [EN]
    Calculates the evolutionary distance between two arrays.

    Processing steps:
    1. Align with NW
    2. Calculate the percentage of identity
    3. Convert to distance: d = -ln(identity)
    (if identity = 1, then d = 0, meaning the arrays are the same)

    Parameters
    ------------
    seq1, seq2 : str
    Raw (unaligned) arrays.
    scoring : ScoringMatrix

    Returns
    --------
    float : distance value (0 = identical, larger values ​​mean further apart)
    
    [TR]
    İki dizi arasındaki evrimsel mesafeyi hesaplar.
 
    İşlem adımları:
      1. NW ile hizala
      2. Özdeşlik yüzdesini hesapla
      3. Mesafeye çevir: d = -ln(identity)
         (identity = 1 ise d = 0, yani aynı diziler)
 
    Parametreler
    ------------
    seq1, seq2 : str
        Ham (hizalanmamış) diziler.
    scoring : ScoringMatrix
 
    Döndürür
    --------
    float : mesafe değeri (0 = özdeş, büyüdükçe daha uzak)
    """

    if scoring is None:
        scoring = ScoringMatrix()

    aligned1, aligned2,_, _,_ = needleman_wunsch(seq1, seq2, scoring)
    identity = percent_identity(aligned1, aligned2)

    # identity 0 ise log tanımsız → çok büyük bir mesafe döndür
    if identity <= 0:
        return float('inf')
    
    return -math.log(identity)

 
def build_distance_matrix(sequences: list[str], scoring: ScoringMatrix = None) -> list[list[float]]:
    """
    [EN]
    Creates an NxN distance matrix for N sequences.

    The matrix is ​​symmetric: distance(i,j) == distance(j,i)
    The diagonal is zero: distance(i,i) == 0

    Parameters
    ------------
    sequences : list[str]
    List of sequences to align.
    scoring : ScoringMatrix

    Returns
    --------
    list[list[float]] : NxN distance matrix
    
    [TR]
    N dizi için NxN mesafe matrisi oluşturur.
 
    Matris simetriktir: distance(i,j) == distance(j,i)
    Köşegen sıfırdır: distance(i,i) == 0
 
    Parametreler
    ------------
    sequences : list[str]
        Hizalanacak dizilerin listesi.
    scoring : ScoringMatrix
 
    Döndürür
    --------
    list[list[float]] : NxN mesafe matrisi
    """

    n = len(sequences)
    matrix = [[0.0] * n for _ in range(n)]

    for i in range(n):
        for j in range(i + 1, n):
            d = compute_distance(sequences[i], sequences[j], scoring)
            matrix[i][j] = d
            matrix[j][i] = d

    return matrix
  