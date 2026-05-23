PCR Primer Automation System

A comprehensive bioinformatics tool for designing high-quality primer pairs for PCR amplification. This pipeline automatically scans DNA sequences, evaluates primer candidates using multiple quality metrics, and generates publication-ready reports and visualizations.

---

## 📑 Table of Contents

- [Overview](#-overview)
- [Features](#-features)
- [Installation](#-installation)
- [Usage](#-usage)
- [Output Files](#-output-files)
- [Graph Explanations](#-graph-explanations)
- [Algorithm Details](#-algorithm-details)
- [Output Columns](#-output-columns)
- [Examples](#-examples)
- [License](#-license)

---

## 🎯 Overview

This pipeline is designed for molecular biologists and bioinformatics researchers who need to design reliable primer pairs for PCR amplification. It processes DNA sequences from FASTA or plain text files and returns ranked primer pairs with comprehensive quality metrics.

### Key Capabilities:
- **High-throughput scanning** of thousands of primer pair candidates
- **Multi-criteria filtering** based on length, GC content, and melting temperature
- **Quality scoring system** (0-10 scale) evaluating homopolymers, poly-N regions, and G/C clamps
- **Self-complementarity analysis** using primer3 thermodynamic calculations
- **Cross-complementarity evaluation** between forward and reverse primers
- **Automated visualization** with 4 publication-ready graphs

---

## ✨ Features

| Feature | Description |
|---------|-------------|
| 🔍 **Smart Scanning** | Grid-based scanning with configurable step size and product ranges |
| 🌡️ **Tm Calculation** | Wallace rule (2×AT + 4×GC) for accurate melting temperature |
| 📊 **GC Content** | Precise GC percentage calculation for each primer |
| 🏆 **Quality Scoring** | Comprehensive 0-10 quality score with penalty system |
| 🔗 **Dimer Analysis** | Homodimer and heterodimer stability via primer3 |
| 📈 **4 Visualizations** | Histograms, scatter plots, KDE curves, and heatmaps |
| 📋 **CSV Export** | Full candidate list with all metrics in primer_candidates.csv |
| 📄 **Text Report** | Top 5 primers summary in top5_report.txt |

---

## 🔧 Installation

### Prerequisites

```bash
# Python 3.8 or higher required
python --version
```

### Required Packages

```bash
# Install all dependencies
pip install biopython primer3-py pandas numpy matplotlib seaborn
```

Or use the requirements file:

```bash
pip install -r requirements.txt
```

### Dependencies Detail

| Package | Version | Purpose |
|---------|---------|---------|
| **BioPython** | ≥1.80 | Sequence handling, GC calculation |
| **primer3-py** | ≥2.0 | Thermodynamic calculations for dimers |
| **pandas** | ≥1.5 | Data manipulation and CSV export |
| **numpy** | ≥1.23 | Numerical operations |
| **matplotlib** | ≥3.6 | Graph generation |
| **seaborn** | ≥0.12 | Statistical visualizations |

---

## 🚀 Usage

### Basic Usage

```python
# Simply run the script - it will automatically look for dna.txt
python primer_design_pipeline.py
```

### Input File

Place your DNA sequence in a file named `dna.txt` in the `/content` directory.

**Supported formats:**
- **FASTA format**: Standard FASTA with sequence header
- **Plain text**: Raw DNA sequence (A, T, G, C, N characters)

### Example dna.txt:

```
> Wheat_BAX_Inhibitor
ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA
TCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA
...
```

### Output Location

All output files are saved to `/content/MINI_PROJECT/`:
- `primer_candidates.csv` - Complete candidate list
- `top5_report.txt` - Top 5 primers summary
- `graph1_penalty.png` - Penalty score distribution
- `graph2_tm_scatter.png` - Tm scatter plot
- `graph3_gc_dist.png` - GC content distribution
- `graph4_heatmap.png` - Feature correlation heatmap

---

## 📁 Output Files

### 1. primer_candidates.csv

The main output containing all primer pair candidates with 23 columns of metrics, sorted by penalty score (best candidates first).

### 2. top5_report.txt

Human-readable summary of the top 5 primer pairs with key metrics:

```
🏆 TOP 5 PRIMER PAIRS
Target: Triticum_BAX_Inhibitor
Sequence Length: 1,234 bp
==================================================

1. Penalty: 0.156
   Forward:  ATGCGATCGATGCGATCGAT
   Reverse:  CGATCGATCGCGATCGTAG
   Tm: 60.2/59.8°C  |  GC: 50.0/48.0%
   Product: 200 bp
   Tm Diff: 0.4°C | Cross 3': 2.3
...
```

### 3. Graph Visualizations (PNG)

## 📊 Graph Explanations

### 📈 Graph 1: Penalty Score Distribution

**File:** `graph1_penalty.png`

![Penalty Score Distribution](public/images/graph1_penalty.png)

**Description:**
This histogram shows the distribution of penalty scores across all analyzed primer pairs.

| Aspect | Explanation |
|--------|-------------|
| **X-axis** | Penalty score (lower is better, range 0-5 typically) |
| **Y-axis** | Count of primer pairs with that score |
| **Color** | Blue bars with KDE (Kernel Density Estimate) overlay |
| **Insight** | Peaks near 0 indicate many high-quality candidates |

**How to Interpret:**
- **Left tail (near 0):** Best primer candidates with minimal penalties
- **Right tail (>2):** Lower quality primers with multiple issues
- **Wide distribution:** Diverse candidate pool
- **Narrow peak:** Consistent primer quality

---

### 🎯 Graph 2: Forward vs Reverse Tm Scatter

**File:** `graph2_tm_scatter.png`

![Tm Scatter Plot](public/images/graph2_tm_scatter.png)

**Description:**
Bivariate scatter plot comparing melting temperatures of forward and reverse primers.

| Aspect | Explanation |
|--------|-------------|
| **X-axis** | Forward primer Tm (°C) |
| **Y-axis** | Reverse primer Tm (°C) |
| **Color** | Penalty score (viridis colormap: purple=low, yellow=high) |
| **Ideal Zone** | Both Tm around 60°C (diagonal area) |
| **Insight** | Points near the diagonal have balanced Tm |

**How to Interpret:**
- **Diagonal line (y=x):** Primers with equal Tm
- **Cluster around 60°C:** Optimal annealing temperature
- **Color gradient:** Yellow points have higher penalty scores
- **Scatter spread:** Tm variation in candidates

**Ideal Scenario:** Tight cluster at coordinates (60, 60)

---

### 🧬 Graph 3: GC Content Distribution

**File:** `graph3_gc_dist.png`

![GC Content KDE](public/images/graph3_gc_dist.png)

**Description:**
Kernel Density Estimation (KDE) plots showing GC content distributions for forward and reverse primers.

| Aspect | Explanation |
|--------|-------------|
| **Red line** | Forward primer GC% distribution |
| **Green line** | Reverse primer GC% distribution |
| **X-axis** | GC percentage (0-100%) |
| **Y-axis** | Density probability |
| **Shaded areas** | 40-65% acceptable range |

**How to Interpret:**
- **Peaks at 50-60%:** Ideal GC content for stable annealing
- **Bimodal distribution:** Two populations of primers
- **Overlap area:** Both primers have similar GC content
- **Tails outside 40-65%:** Potentially problematic primers

**Why GC Matters:**
- **Too low (<40%):** Weak binding, unstable duplex
- **Too high (>65%):** Non-specific binding, secondary structures
- **Optimal:** 50-60% for most applications

---

### 🔥 Graph 4: Feature Correlation Heatmap

**File:** `graph4_heatmap.png`

![Correlation Heatmap](public/images/graph4_heatmap.png)

**Description:**
Correlation matrix heatmap showing relationships between all numeric features.

| Aspect | Explanation |
|--------|-------------|
| **Color scale** | Red = positive correlation, Blue = negative correlation |
| **Values** | Pearson correlation coefficient (-1 to +1) |
| **Diagonal** | Always 1.0 (self-correlation) |

**Key Correlations to Look For:**

| Correlation | Interpretation |
|-------------|----------------|
| **fwd_tm ↔ rev_tm** | High positive = balanced Tm design |
| **penalty_score ↔ cross_3p** | Positive = cross stability affects quality |
| **fwd_gc ↔ fwd_tm** | Positive = GC influences Tm (expected) |
| **tm_diff ↔ penalty** | Positive = larger Tm diff = higher penalty |

**How to Read:**
- **+1.0 (dark red):** Perfect positive correlation
- **0.0 (white):** No correlation
- **-1.0 (dark blue):** Perfect negative correlation

---

## 🔬 Algorithm Details

### Melting Temperature (Wallace Rule)

```
Tm = 2 × (A + T) + 4 × (G + C)
```

This simplified formula is widely used for primers under 25 bp and assumes:
- 50mM Na⁺ concentration
- 50nM primer concentration

### Quality Score Calculation

The quality score (0-10) penalizes problematic features:

| Feature | Penalty | Reason |
|---------|---------|--------|
| **Homopolymer (5+ repeats)** | -4.0 | Replication slippage |
| **Homopolymer (4 repeats)** | -2.0 | Potential instability |
| **G/C clamp at 3' end** | +0.5 | Stable annealing |
| **Poly-N regions** | -3.0 | Ambiguous binding |
| **Dinucleotide repeats at 3'** | -1.5 | Mispriming risk |

### Penalty Score Formula

```
penalty = (|fwd_tm - 60| × 0.05) + (|rev_tm - 60| × 0.05)
        + (|fwd_gc - 50| × 0.03) + (|rev_gc - 50| × 0.03)
        + ((10 - fwd_quality) × 0.1) + ((10 - rev_quality) × 0.1)
```

**Lower penalty = Better primer pair**

### Constraints

| Parameter | Minimum | Maximum | Optimal Target |
|-----------|---------|---------|----------------|
| **Primer Length** | 18 bp | 25 bp | 20 bp |
| **GC Content** | 40% | 65% | 50% |
| **Melting Temp** | 55°C | 65°C | 60°C |
| **Product Size** | 150 bp | 400 bp | 200-300 bp |

---

## 📋 Output Columns

| Column | Type | Description |
|--------|------|-------------|
| `forward` | string | Forward primer sequence (5'→3') |
| `reverse` | string | Reverse primer sequence (5'→3') |
| `Tm` | string | Forward/Reverse Tm (e.g., "60.2/59.8") |
| `GC` | string | Forward/Reverse GC% (e.g., "50.0/48.0") |
| `fwd_tm` | float | Forward primer Tm in °C |
| `rev_tm` | float | Reverse primer Tm in °C |
| `fwd_gc` | float | Forward primer GC% |
| `rev_gc` | float | Reverse primer GC% |
| `product_size` | int | Amplicon length in bp |
| `penalty_score` | float | Combined penalty (lower=better) |
| `tm_diff` | float | Absolute Tm difference between primers |
| `fwd_start` | int | Start position of forward primer |
| `fwd_end` | int | End position of forward primer |
| `rev_start` | int | Start position of reverse primer |
| `rev_end` | int | End position of reverse primer |
| `fwd_self_any` | float | Forward self-complementarity (any) |
| `fwd_self_3p` | float | Forward self-complementarity (3' end) |
| `rev_self_any` | float | Reverse self-complementarity (any) |
| `rev_self_3p` | float | Reverse self-complementarity (3' end) |
| `cross_any` | float | Cross-complementarity (any) |
| `cross_3p` | float | Cross-complementarity (3' end) |
| `fwd_quality` | float | Forward primer quality score (0-10) |
| `rev_quality` | float | Reverse primer quality score (0-10) |

---

## 💡 Examples

### Running the Pipeline

```bash
# Clone the repository
git clone https://github.com/yourusername/primer-design-pipeline.git
cd primer-design-pipeline

# Install dependencies
pip install -r requirements.txt

# Prepare your sequence file
cp your_sequence.fasta dna.txt

# Run the pipeline
python primer_design_pipeline.py
```

### Sample Output

```
🚀 PRIMER DESIGN PIPELINE v2.0
=============================================

🔍 BACKGROUND SEARCH: Scanning for dna.txt...
✅ LOADING: dna.txt
📏 SEQUENCE: Triticum_BAX_Inhibitor (2,456 bp)

🔬 SCANNING HIGH-QUALITY PRIMER PAIRS...
✅ Generated 347 candidate pairs

💾 primer_candidates.csv → /content/MINI_PROJECT/primer_candidates.csv
📊 top5_report.txt → /content/MINI_PROJECT/top5_report.txt
📈 GENERATING GRAPHS...
✅ All graphs generated!

🎯 BEST PRIMER:
FWD: ATGCGATCGATGCGATCGAT
REV: CGATCGATCGCGATCGTAG
Tm: 60.2/59.8°C | GC: 50.0/48.0%
Size: 200 bp | Score: 0.156
Tm Diff: 0.4°C | Cross Any: 3.2

✅ PIPELINE COMPLETED SUCCESSFULLY!
```

### Loading Results in Python

```python
import pandas as pd

# Load the results
df = pd.read_csv('/content/MINI_PROJECT/primer_candidates.csv')

# Get top 10 primers
top10 = df.head(10)

# Filter by GC content
high_gc = df[(df['fwd_gc'] > 50) & (df['rev_gc'] > 50)]

# Find primers with balanced Tm
balanced = df[df['tm_diff'] < 1.0]
```

---

## 🧪 Workflow Diagram

```
┌─────────────────────────────────────────────────────────────┐
│                    INPUT SEQUENCE                           │
│                   (dna.txt file)                            │
└─────────────────────────┬───────────────────────────────────┘
                          │
                          ▼
┌─────────────────────────────────────────────────────────────┐
│              BACKGROUND SEARCH & PARSING                    │
│         (FASTA or plain text format support)                │
└─────────────────────────┬───────────────────────────────────┘
                          │
                          ▼
┌─────────────────────────────────────────────────────────────┐
│               PRIMER PAIR SCANNING                          │
│   • Length: 18-25 bp                                        │
│   • GC: 40-65%                                              │
│   • Tm: 55-65°C                                             │
│   • Quality ≥ 6                                             │
└─────────────────────────┬───────────────────────────────────┘
                          │
                          ▼
┌─────────────────────────────────────────────────────────────┐
│              QUALITY & PENALTY SCORING                      │
│   • Melting temperature analysis                            │
│   • GC content calculation                                  │
│   • Homopolymer/repeat detection                           │
│   • Self-complementarity (primer3)                         │
│   • Cross-complementarity analysis                          │
└─────────────────────────┬───────────────────────────────────┘
                          │
                          ▼
┌─────────────────────────────────────────────────────────────┐
│              SORTING & RANKING                             │
│          (By penalty score, ascending)                     │
└─────────────────────────┬───────────────────────────────────┘
                          │
                          ▼
┌─────────────────────────────────────────────────────────────┐
│                   OUTPUT GENERATION                        │
├─────────────────────────────────────────────────────────────┤
│  📄 primer_candidates.csv    →  500 best candidates        │
│  📋 top5_report.txt          →  Top 5 detailed summary     │
│  📊 graph1_penalty.png       →  Penalty distribution       │
│  📈 graph2_tm_scatter.png    →  Tm correlation scatter     │
│  📉 graph3_gc_dist.png       →  GC KDE curves              │
│  🔥 graph4_heatmap.png       →  Feature correlations       │
└─────────────────────────────────────────────────────────────┘
```

---

## ⚙️ Customization

### Modify Scanning Parameters

Edit the `constraints` dictionary in `scan_primer_pairs()`:

```python
constraints = {
    'length': (18, 25),      # Primer length range
    'gc': (40, 65),          # GC content range (%)
    'tm': (55, 65),          # Melting temperature range (°C)
    'product': (150, 400)    # Amplicon size range (bp)
}
```

### Adjust Step Size

For finer scanning (more candidates, slower):

```python
step_size = 10  # Instead of 30
```

### Change Output Limit

```python
self.candidates = candidates[:1000]  # Top 1000 instead of 500
```

---

## 📄 License

This project is licensed under the MIT License.

---

## 🙏 Acknowledgments

- **BioPython** - For sequence handling utilities
- **primer3-py** - For thermodynamic calculations
- **seaborn** - For statistical visualizations

Author
Mahnoor - m96837320@gmail.com
