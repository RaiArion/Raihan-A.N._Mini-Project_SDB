import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.gridspec import GridSpec

files_bakteri = {
    'Ideonella sakaiensis': 'bakteri1.fasta',
    'Pseudomonas putida': 'bakteri2.fasta',
    'Alcanivorax borkumensis': 'bakteri3.fasta',
    'Deinococcus radiodurans': 'bakteri4.fasta',
    'Cupriavidus metallidurans': 'bakteri5.fasta'
}
data_pipeline = []

def baca_fasta(file_path):
    """Membaca sekuens DNA dari file FASTA, mengabaikan baris header."""
    sekuens = ""
    if os.path.exists(file_path):
        with open(file_path, 'r') as f:
            for line in f:
                if not line.startswith('>'):
                    sekuens += line.strip().upper()
    else:
        print(f"Peringatan: File {file_path} tidak ditemukan!")
    return sekuens

for nama_bakteri, file_name in files_bakteri.items():
    seq = baca_fasta(file_name)
    if seq:
        data_pipeline.append({'bakteri': nama_bakteri, 'sekuens': seq})

hasil_analisis = []
for item in data_pipeline:
    nama = item['bakteri']
    seq = item['sekuens']
    total_basa = len(seq)
    frekuensi = {
        'A': seq.count('A'),
        'T': seq.count('T'),
        'G': seq.count('G'),
        'C': seq.count('C')
    }
    gc_content = ((frekuensi['G'] + frekuensi['C']) / total_basa) * 100 if total_basa > 0 else 0
    # Estimasi berat molekul DNA untai tunggal (g/mol)
    # Rumus: (A * 313.2) + (T * 304.2) + (G * 329.2) + (C * 289.2)
    berat_molekul = (
        frekuensi['A'] * 313.2 +
        frekuensi['T'] * 304.2 +
        frekuensi['G'] * 329.2 +
        frekuensi['C'] * 289.2
    )
    hasil_analisis.append({
        'Bakteri': nama,
        'Frekuensi_A': frekuensi['A'],
        'Frekuensi_T': frekuensi['T'],
        'Frekuensi_G': frekuensi['G'],
        'Frekuensi_C': frekuensi['C'],
        'GC_Content(%)': round(gc_content, 2),
        'Berat_Molekul(g/mol)': round(berat_molekul, 2)
    })

df = pd.DataFrame(hasil_analisis)
df_sorted = df.sort_values(by='GC_Content(%)', ascending=False).reset_index(drop=True)

top_3 = df_sorted.head(3)
print("\n" + "="*45)
print("     === TOP 3 BAKTERI BERDASARKAN GC CONTENT ===")
print("="*45)
print(top_3[['Bakteri', 'GC_Content(%)', 'Berat_Molekul(g/mol)']].to_string(index=False))
print("="*45 + "\n")

df_sorted.to_csv('hasil_analisis_bioremediasi.csv', index=False)
print("Sukses: Hasil analisis lengkap berhasil disimpan ke 'hasil_analisis_bioremediasi.csv'!")

print("Sedang membuat grafik visualisasi...")

COLOR_PRIMARY = '#2a78d6'
COLOR_SECONDARY = '#1baf7a'
COLOR_DIM = '#d3d1c7'
COLOR_TEXT = '#2c2c2a'
COLOR_MUTED = '#888780'
COLOR_GRID = '#e1e0d9'
COLOR_BG = '#fafafa'

labels = [name.replace(' ', '\n') for name in df_sorted['Bakteri']]
n = len(df_sorted)

def make_bar_colors(base, dim, n):
    """Warna penuh untuk batang tertinggi, redup untuk sisanya."""
    return [base if i == 0 else dim for i in range(n)]

plt.rcParams.update({
    'font.family': 'sans-serif',
    'axes.spines.top': False,
    'axes.spines.right': False,
    'axes.spines.left': False,
    'axes.spines.bottom': False,
})

fig = plt.figure(figsize=(11, 8), facecolor=COLOR_BG)
fig.suptitle(
    'Analisis Bakteri Agen Bioremediasi',
    fontsize=13, fontweight='medium', color=COLOR_TEXT,
    x=0.06, ha='left', y=0.97
)
fig.text(
    0.06, 0.93,
    'GC Content dan Berat Molekul DNA Untai Tunggal',
    fontsize=10, color=COLOR_MUTED, ha='left'
)

gs = GridSpec(2, 1, figure=fig, hspace=0.55, top=0.88, bottom=0.08, left=0.06, right=0.97)

# Panel atas: GC Content
ax1 = fig.add_subplot(gs[0])
ax1.set_facecolor(COLOR_BG)

bar_colors_gc = make_bar_colors(COLOR_PRIMARY, '#b5d4f4', n)
bars = ax1.bar(
    range(n),
    df_sorted['GC_Content(%)'],
    color=bar_colors_gc,
    width=0.5,
    zorder=3
)

for bar, val in zip(bars, df_sorted['GC_Content(%)']):
    ax1.text(
        bar.get_x() + bar.get_width() / 2,
        bar.get_height() + 0.4,
        f'{val:.1f}%',
        ha='center', va='bottom',
        fontsize=9, color=COLOR_MUTED
    )

ax1.set_xticks(range(n))
ax1.set_xticklabels(labels, fontsize=9, color=COLOR_TEXT, linespacing=1.3)
ax1.tick_params(axis='x', length=0, pad=6)
ax1.tick_params(axis='y', labelsize=9, colors=COLOR_MUTED, length=0)
ax1.set_ylabel('GC Content (%)', fontsize=9, color=COLOR_MUTED, labelpad=8)
ax1.yaxis.set_major_formatter(mticker.FormatStrFormatter('%.0f%%'))
ax1.set_ylim(0, df_sorted['GC_Content(%)'].max() * 1.18)
ax1.grid(axis='y', color=COLOR_GRID, linewidth=0.8, zorder=0)
ax1.set_axisbelow(True)

# Panel bawah: Berat Molekul
ax2 = fig.add_subplot(gs[1])
ax2.set_facecolor(COLOR_BG)

bar_colors_mw = make_bar_colors(COLOR_SECONDARY, '#9fe1cb', n)
bars2 = ax2.bar(
    range(n),
    df_sorted['Berat_Molekul(g/mol)'],
    color=bar_colors_mw,
    width=0.5,
    zorder=3
)

for bar, val in zip(bars2, df_sorted['Berat_Molekul(g/mol)']):
    ax2.text(
        bar.get_x() + bar.get_width() / 2,
        bar.get_height() + df_sorted['Berat_Molekul(g/mol)'].max() * 0.01,
        f'{val/1e6:.2f}M',
        ha='center', va='bottom',
        fontsize=9, color=COLOR_MUTED
    )

ax2.set_xticks(range(n))
ax2.set_xticklabels(labels, fontsize=9, color=COLOR_TEXT, linespacing=1.3)
ax2.tick_params(axis='x', length=0, pad=6)
ax2.tick_params(axis='y', labelsize=9, colors=COLOR_MUTED, length=0)
ax2.set_ylabel('Berat Molekul (g/mol)', fontsize=9, color=COLOR_MUTED, labelpad=8)
ax2.yaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: f'{x/1e6:.1f}M'))
ax2.set_ylim(0, df_sorted['Berat_Molekul(g/mol)'].max() * 1.18)
ax2.grid(axis='y', color=COLOR_GRID, linewidth=0.8, zorder=0)
ax2.set_axisbelow(True)

# Label panel
for ax, label, color in [(ax1, 'GC Content (%)', COLOR_PRIMARY), (ax2, 'Berat Molekul (g/mol)', COLOR_SECONDARY)]:
    ax.set_title(label, fontsize=10, color=color, fontweight='medium', loc='left', pad=8)

plt.savefig('grafik_gc_dan_mw.png', dpi=300, bbox_inches='tight', facecolor=COLOR_BG)
print("Sukses: Gambar grafik visualisasi berhasil disimpan dengan nama 'grafik_gc_dan_mw.png'!")
plt.show()