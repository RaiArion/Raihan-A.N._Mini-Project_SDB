#  Pipeline Analisis Sekuens DNA 16S rRNA: Karakterisasi Genomik Bakteri Agen Bioremediasi

Repository ini dibuat untuk memenuhi tugas **Mini Project Mata Kuliah Struktur Data Bioinformatika (BIF1223)** - IPB University.

**Anggota Kelompok/Individu:**
* Nama: Raihan Arionanda Nasution
* NIM: G0401241033
* Departemen: Bioinformatika

---

##  Latar Belakang & Topik Analisis
Projek ini melakukan analisis komparatif terhadap sekuens **16S rRNA gene (partial sequence)** dari **5 bakteri agen bioremediasi** yang diunduh langsung dari database NCBI untuk mengevaluasi stabilitas termal genom (*GC Content*) serta dimensi fisik materi genetiknya (*Berat Molekul*):
1. **Ideonella sakaiensis** (Kemampuan mendegradasi plastik PET)
2. **Pseudomonas putida** (Kemostasis polutan organik/hidrokarbon)
3. **Alcanivorax borkumensis** (Pendegradasi utama tumpahan minyak di laut)
4. **Deinococcus radiodurans** (Resisten radiasi & penyerap limbah radioaktif)
5. **Cupriavidus metallidurans** (Pendetoks dan akumulator logam berat)

---

## 🛠️ Implementasi Struktur Data (Pipeline Flow)
Sesuai ketentuan, program Python (`main.py`) dibangun menggunakan struktur data dasar mandiri:
* **List:** Pipeline utama untuk menampung seluruh objek sekuens dari kelima berkas FASTA.
* **Dictionary:** Memproses hitung frekuensi individual nukleotida (A, T, G, C) secara efisien.
* **Dataframe & Sorting:** Mengurutkan data berdasarkan *GC Content* secara menurun (*descending*) untuk menyaring **Top 3 sekuens terbaik**.

---

## Analisis Tambahan: Estimasi Berat Molekul (*Molecular Weight*)
Selain parameter *GC Content*, projek ini menambahkan fitur kalkulasi **Berat Molekul Untai Tunggal (Single-Stranded DNA Molecular Weight)** menggunakan rumus biokimia empiris:
$$\text{Berat Molekul} = (A \times 313.2) + (T \times 304.2) + (G \times 329.2) + (C \times 289.2)$$

---

## Hasil Analisis Lengkap (Urutan GC Content)
Berdasarkan hasil eksekusi program pada data sekuens nyata, berikut adalah profil lengkap kelima bakteri:

| Peringkat | Spesies Bakteri | Frekuensi (A, T, G, C) | GC Content (%) | Berat Molekul (g/mol) |
| :---: | :--- | :---: | :---: | :---: |
| **1** | *Ideonella sakaiensis* | (367, 284, 487, 355) | **56.40%** | 464,323.6 |
| **2** | *Deinococcus radiodurans* | (362, 289, 474, 341) | **55.59%** | 455,950.2 |
| **3** | *Alcanivorax borkumensis* | (374, 303, 476, 341) | **54.69%** | 464,625.8 |
| 4 | *Cupriavidus metallidurans* | (378, 299, 477, 338) | 54.59% | 464,123.4 |
| 5 | *Pseudomonas putida* | (386, 321, 476, 344) | 53.70% | 474,727.4 |

---

## Visualisasi Grafik
Hasil visualisasi hubungan grafik batang kontras tinggi antara persentase *GC Content* dan grafik *Berat Molekul* telah disimpan otomatis di dalam repository ini dengan nama berkas:
 `grafik_gc_dan_mw.png`

---

##  Cara Menjalankan Program
1. Pastikan dependensi telah terinstal:
   ```bash
   pip install pandas matplotlib
