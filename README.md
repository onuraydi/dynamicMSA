# dynamicMSA

Dinamik Programlama (Needleman-Wunsch + UPGMA) tabanlı **Çoklu Dizi Hizalaması** kütüphanesi.

## Kurulum

```bash
pip install dynamicMSA
```

## Hızlı Kullanım

```python
from dynamicMSA import align

sequences = ["ACGTACGT", "ACGTTCGT", "ACGT-CGT"]
names     = ["Dizi1", "Dizi2", "Dizi3"]

result = align(sequences, names=names, visualize=True)

print(result["aligned"])   # hizalanmış diziler
print(result["score"])     # toplam hizalama skoru
```

## Özellikler

- **Needleman-Wunsch** ile global ikili hizalama
- **UPGMA** ile kılavuz ağaç oluşturma
- **Progressive Alignment** ile çoklu dizi hizalaması
- **Matplotlib** ile 3 farklı görselleştirme:
  - DP matrisi ısı haritası (traceback yolu ile)
  - UPGMA kılavuz ağacı (dendrogram)
  - Renkli hizalama tablosu (korunmuş pozisyonlar vurgulanır)

## Parametreler

| Parametre | Açıklama | Varsayılan |
|-----------|----------|-----------|
| `match` | Eşleşme skoru | `+1` |
| `mismatch` | Uyumsuzluk cezası | `-1` |
| `gap` | Boşluk cezası | `-2` |
| `visualize` | Grafikleri göster | `True` |
| `save_figures` | Grafikleri kaydet | `False` |

## Algoritma

```
1. Her dizi çifti için Needleman-Wunsch → mesafe matrisi
2. UPGMA algoritması → kılavuz ağaç
3. Progressive Alignment → final MSA
```

## Lisans

MIT