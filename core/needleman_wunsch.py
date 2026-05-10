"""
[EN]
Needleman-Wunsch algorithm for global array alignment.

The algorithm has 3 stages:
1. INITIALIZATION: Fill the first row and column of the matrix with gap penalties.
2. FILLING: Fill each cell by looking at its three neighbors (top, left, diagonal).
3. TRACEBACK: Find the alignment by tracing from the bottom right of the matrix to the top left.

[TR]
Global dizi hizalaması için Needleman-Wunsch algoritması.
 
Algoritmanın 3 aşaması:
  1. İNİSYALİZASYON : Matrisin ilk satır ve sütununu gap cezalarıyla doldur
  2. DOLDURMA       : Her hücreyi üç komşusuna bakarak doldur (yukarı, sol, çapraz)
  3. TRACEBACK      : Matrisin sağ-altından başlayıp sola-yukarı izleyerek hizalamayı bul
"""

import numpy as np
from dynamicMSA.utils.scoring import ScoringMatrix

def needleman_wunsch(
        seq1: str,
        seq2: str,
        scoring: ScoringMatrix = None
) -> tuple[str, str, float, np.ndarray, list] :
    """
    [EN]
    Parameters
    ------------
    seq1, seq2 : str
    DNA/RNA/protein sequences to be aligned.
    scoring : ScoringMatrix
    Score matrix. Default values ​​are used if not provided.

    Returns
    --------
    aligned1 : str — first aligned sequence (gaps are indicated by '-')
    aligned2 : str — second aligned sequence
    score : float — total alignment score
    matrix : np.ndarray — populated DP matrix (for visualization)
    path : list — traceback path [(i,j), ...] (for visualization)
        

    [TR]
    Parametreler
    ------------
    seq1, seq2 : str
        Hizalanacak DNA/RNA/protein dizileri.
    scoring : ScoringMatrix
        Skor matrisi. Verilmezse varsayılan değerler kullanılır.
 
    Döndürür
    --------
    aligned1 : str   — hizalanmış birinci dizi (gap'ler '-' ile gösterilir)
    aligned2 : str   — hizalanmış ikinci dizi
    score    : float — toplam hizalama skoru
    matrix   : np.ndarray — doldurulan DP matrisi (görselleştirme için)
    path     : list  — traceback yolu [(i,j), ...] (görselleştirme için)
    """

    if scoring is None:
        scoring = ScoringMatrix()

    n = len(seq1)
    m = len(seq2)

    # ---------------------------------------------------------------
    # AŞAMA 1: Matris oluşturma ve inisyalizasyon
    # Satırlar seq1'e, sütunlar seq2'ye karşılık gelir.
    # İlk satır ve sütun kümülatif gap cezalarıyla başlar.
    # ---------------------------------------------------------------


    matrix = np.zeros((n+1, m+1), dtype = float)

    for i in range(n+1):
        matrix[i][0] = i * scoring.gap # ilk satırdaki 0'lar
    for j in range(m+1):
        matrix[0][j] = j* scoring.gap # ilk sütundaki 0'lar


    # ---------------------------------------------------------------
    # AŞAMA 2: Matrisi doldurma
    # Her hücre için üç seçenek değerlendirilir:
    #   - Çapraz (diagonal) : iki karakteri eşleştir
    #   - Yukarıdan         : seq1'e gap ekle
    #   - Soldan            : seq2'ye gap ekle
    # En yüksek değer seçilir.
    # ---------------------------------------------------------------

    for i in range(1, n+1):
        for j in range(1, m+1):
            diagonal = matrix[i-1][j-1] + scoring.score(seq1[i-1], seq2[j-1])
            up = matrix[i-1][j] + scoring.gap
            left = matrix[i][j-1] + scoring.gap

            matrix[i][j] = max(diagonal, up, left)

    # ---------------------------------------------------------------
    # AŞAMA 3: Traceback
    # Sağ-alt köşeden başlayıp hangi hücre seçildiyse oraya giderek
    # hizalamayı geriye doğru inşa ederiz.
    # ---------------------------------------------------------------
    
    aligned1, aligned2 = [], []
    path = []
    i,j = n,m

    while i > 0 or j > 0:
        path.append((i,j))
        current = matrix[i][j]

        if i > 0 or j > 0:
            diag_val = matrix[i-1][j-1] + scoring.score(seq1[i-1],seq2[j-1])
        else:
            diag_val = float('-inf')

        if i > 0 and j > 0 and current == diag_val:
            # çapraz: her iki diziden de karakter al
            aligned1.append(seq1[i-1])
            aligned2.append(seq2[j-1])
            i -= 1
            j -= 1
        elif i > 0 and current == matrix[i-1][j] + scoring.gap:
            # Yukarıdan: seq1'den karakter, seq2'ye gap
            aligned1.append(seq1[i-1])
            aligned2.append('-')
            i -= 1
        else:
            # Soldan: seq1'e gap, seq2'den karakter
            aligned1.append('-')
            aligned2.append(seq2[j-1])
            j -= 1
    
    path.append((0,0))

    # Geriye doğru inşa ettiğimiz için ters çevir
    aligned1 = ''.join(reversed(aligned1))
    aligned2 = ''.join(reversed(aligned2))
    path = list(reversed(path))

    return aligned1, aligned2, matrix[n][m], matrix, path