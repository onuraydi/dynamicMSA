"""
[EN]
Defines the score values ​​to be used during alignment.

Key concepts:
- MATCH: Points earned if two nucleotides/amino acids are the same.
- MISMATCH: Penalty (negative) if they are different.
- GAP: Penalty (negative) when a gap (-) is added to a sequence.

These values ​​directly affect how the Needleman-Wunsch matrix is ​​filled.

[TR]
Hizalama sırasında kullanılacak skor değerlerini tanımlar.
 
Temel kavramlar:
- MATCH     : İki nükleotid/amino asit aynıysa kazanılan puan
- MISMATCH  : Farklılarsa verilen ceza (negatif)
- GAP       : Bir diziye boşluk (-) eklendiğinde verilen ceza (negatif)
 
Bu değerler Needleman-Wunsch matrisinin nasıl dolacağını doğrudan etkiler.
"""

class ScoringMatrix:
    """
    [EN]
    A simple score matrix.

    Parameters
    ------------
    match : int
    Points awarded when the same characters match. (default: +1)
    mismatch : int
    Penalty awarded when different characters match. (default: -1)
    gap : int
    Penalty awarded when a space (indel) is created. (default: -2)

    [TR]
    Basit bir skor matrisi.
 
    Parametreler
    ------------
    match : int
        Aynı karakterler eşleştiğinde verilecek puan. (varsayılan: +1)
    mismatch : int
        Farklı karakterler eşleştiğinde verilecek ceza. (varsayılan: -1)
    gap : int
        Boşluk (indel) açıldığında verilecek ceza. (varsayılan: -2)
    """

    def __init__(self, match: int = 1, mismatch: int = -1, gap: int = -2):
        self.match = match
        self.mismatch = mismatch
        self.gap = gap

    def score(self, a: str, b: str) -> int:
        """
        [EN]
        Compares two characters and returns a score.

        Parameters
        ------------
        a, b : str
        Single characters to compare. '-' means gap.

        Returns
        --------
        int : match, mismatch, or gap score
        
        [TR]
        İki karakteri karşılaştırır ve skor döndürür.
 
        Parametreler
        ------------
        a, b : str
            Karşılaştırılacak tek karakterler. '-' gap anlamına gelir.
 
        Döndürür
        --------
        int : match, mismatch veya gap skoru
        """

        if a == '-' or b == '-':
            return self.gap
        return self.match if a == b else self.mismatch
    
    def __repr__(self):
        return (f"ScoringMatrix(match={self.match}, "
                f"mismatch={self.mismatch}, gap={self.gap})")