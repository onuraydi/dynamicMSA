"""
[EN]
UPGMA (Unweighted Pair Group Method with Arithmetic Mean) algorithm
creates a guide tree.

Why do we need a guide tree?
In multiple alignments, we need to know which sequences to align first.
UPGMA applies the "merge the most similar pair first" rule.

How does UPGMA work?
1. Each sequence starts as its own set: {A}, {B}, {C}, {D}
2. Find the two sets with the smallest distance → merge
3. Update the distance of the new set relative to the others (take the average)
4. Repeat until only one set remains

Output format:
Each node is a dict:
{"left": ..., "right": ..., "height": float, "sequences": [int, ...]}
Leaf nodes:
{"index": int, "name": str, "sequences": [int]}

[TR]
UPGMA (Unweighted Pair Group Method with Arithmetic Mean) algoritması ile
kılavuz ağaç (guide tree) oluşturur.
 
Neden kılavuz ağaca ihtiyacımız var?
  Çoklu hizalamada hangi dizileri önce hizalayacağımızı bilmemiz gerekir.
  UPGMA bize "en benzer ikili önce birleştir" kuralını uygular.
 
UPGMA nasıl çalışır?
  1. Her dizi kendi kümesi olarak başlar: {A}, {B}, {C}, {D}
  2. En küçük mesafeli iki kümeyi bul → birleştir
  3. Yeni kümenin mesafesini diğerlerine göre güncelle (ortalama al)
  4. Tek küme kalana kadar tekrarla
 
Çıktı formatı:
  Her düğüm bir dict'tir:
    {"left": ..., "right": ..., "height": float, "sequences": [int, ...]}
  Yaprak düğümler:
    {"index": int, "name": str, "sequences": [int]}
"""

import copy

def upgma(distance_matrix: list[list[float]], names: list[str] = None) -> dict:
    """
    [EN]
    Performs hierarchical clustering with the UPGMA algorithm.

    Parameters
    ------------
    distance_matrix : list[list[float]]
    NxN matrix generated with build_distance_matrix().
    names : list[str]
    Array names (optional, for visualization).

    Returns
    --------
    dict : Root node of the tree (in nested dict structure)

    [TR]
    UPGMA algoritması ile hiyerarşik kümeleme yapar.
 
    Parametreler
    ------------
    distance_matrix : list[list[float]]
        build_distance_matrix() ile üretilen NxN matris.
    names : list[str]
        Dizi isimleri (opsiyonel, görselleştirme için).
 
    Döndürür
    --------
    dict : Ağacın kök düğümü (nested dict yapısında)
    """
    n = len(distance_matrix)

    if names is None:
        names = [f"seq{i}" for i in range(n)]

        # Çalışma kopyası: algoritma matrisi değiştirir
        dist = copy.deepcopy(distance_matrix)

        # Her küme başlangıçta tek bir yaprak düğüm
        clusters = [
            {"index": i, "name": names[i], "sequences": [i], "height": 0.0}
            for i in range(n)
        ]

        # Birleştirme adımlarını kaydet (görselleştirme için)
        merge_history = []

        while len(clusters) > 1:
            # En Küçük mesafeli çifti bul
            min_dist = float('inf')
            min_i, min_j = 0,1

            for i in range(len(clusters)):
                for j in range(i + 1, len(clusters)):
                    ci = clusters[i]["sequences"]
                    cj = clusters[j]["sequences"]

                    # küme ortalama mesafesi
                    avg = _average_distance(dist, ci,cj)

                    if avg < min_dist:
                        min_dist = avg
                        min_i, min_j = i, j

            # yeni iç düğüm oluştur
            # height = birleşme noktasının dallanma yüksekliği
            left == clusters[min_i]
            right = clusters[min_j]
            height = min_dist / 2.0

            new_node = {
                "left": left,
                "right": right,
                "height": height,
                "sequences": left["sequences"] + right["sequences"]
            }

            merge_history.append({
                "merged": (left.get("name","node"), right.get("name","node")),
                "distance": min_dist,
                "height": height
            })


            # Kümeleri güncelle: min_i ve min_j'yi kaldır, yenisini ekle
            # (önce büyük indeksi kaldır, küçük indeks kaymasın)
            clusters.pop(max(min_i, min_j))
            clusters.pop(min(min_i, min_j))
            clusters.append(new_node)

        root = clusters[0]
        root["merge_history"] = merge_history
        return root
    
def _average_distance(
    dist: list[list[float]],
    group_a: list[int],
    group_b: list[int]
) -> float:
    """
    [EN]
    It calculates the average distance between two sets.
    This is the "unweighted" part of UPGMA: each pair has an equal weight.

    [TR]
    İki küme arasındaki ortalama mesafeyi hesaplar.
    UPGMA'nın "unweighted" kısmı budur: her çiftin ağırlığı eşittir.
    """
    total = 0.0
    count = 0

    for i in group_a:
        for j in group_b:
            total += dist[i][j]
            count += 1

    return total / count if count > 0 else 0.0

def get_merge_order(root: dict) -> list[list[int]]:
    """
    [EN]
    Extracts the alignment order starting from the root node.

    Returns
    --------
    list[list[int]] : List of array indices to be joined at each step
    Example: [[0,1], [2], [3]] → first align 0 and 1, then add 2, then 3

    [TR]
    Kök düğümden başlayarak hizalama sırasını çıkarır.
 
    Döndürür
    --------
    list[list[int]] : Her adımda birleştirilecek dizi indekslerinin listesi
      Örnek: [[0,1], [2], [3]] → önce 0 ve 1 hizalanır, sonra 2 eklenir, sonra 3
    """
    order = []
    _traverse(root, order)
    return order

def _traverse(node: dict, order: list):
    """
    [EN]
    It travels through the tree in post-order, collecting the order in which the trees merge.

    [TR]
    Ağacı post-order gezerek birleşme sırasını toplar.
    """
    if "index" in node:
        # yaprak düğüm
        order.append([node["index"]])
        return
    
    _traverse(node["left"], order)
    _traverse(node["right"], order)
    order.append(node["sequences"])