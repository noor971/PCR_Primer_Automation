import os
import sys
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
import primer3
from pathlib import Path
import re
import warnings
warnings.filterwarnings('ignore')

# FIXED PATHS
INPUT_DIR = Path(r"/content")
OUTPUT_DIR = Path(r"/content/MINI_PROJECT")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Graph imports
import matplotlib.pyplot as plt
import seaborn as sns


class CompletePrimerPipeline:
    """Complete primer design pipeline with quality scoring"""

    def __init__(self):
        self.sequence = ""
        self.seq_name = "Triticum_BAX_Inhibitor"
        self.candidates = []

    def background_search_sequence(self):
        """Load sequence from dna.txt file"""
        print("🔍 BACKGROUND SEARCH: Scanning for dna.txt...")

        fasta_file_path = INPUT_DIR / "dna.txt" # Create a Path object
        fasta_file_str = str(fasta_file_path) # Convert to string for os.path.exists and file operations

        if not os.path.exists(fasta_file_str): # Use os.path.exists for robust check
            print(f"❌ dna.txt not found in {INPUT_DIR}")
            print("💡 Please ensure 'dna.txt' is uploaded to the /content directory.")
            return False

        print(f"✅ LOADING: {fasta_file_path.name}")

        # Use try-except for robust loading as a plain sequence if FASTA parsing fails
        try:
            records = list(SeqIO.parse(fasta_file_str, "fasta")) # Use string path for SeqIO.parse
            if records:
                record = records[0]
                self.sequence = str(record.seq).upper()
                self.seq_name = fasta_file_path.stem # Use Path object for .stem
                print(f"📏 SEQUENCE: {self.seq_name} ({len(self.sequence):,} bp)")
                return True
            else:
                # Fallback for plain text DNA files
                with open(fasta_file_str, 'r') as f: # Use string path for open
                    self.sequence = f.read().strip().upper()
                self.seq_name = fasta_file_path.stem # Use Path object for .stem
                print(f"📏 SEQUENCE (Plain Text): {self.seq_name} ({len(self.sequence):,} bp)")
                return True
        except Exception as e:
            print(f"Error loading sequence: {e}")
            return False

    def calculate_tm_basic(self, seq):
        """Wallace rule: Tm = 2(A+T) + 4(G+C)"""
        return 2 * (seq.count('A') + seq.count('T')) + 4 * (seq.count('G') + seq.count('C'))

    def gc_content(self, seq):
        """GC percentage"""
        return round(gc_fraction(seq) * 100, 1)

    def quality_score(self, seq):
        """Comprehensive quality scoring (0-10, higher=better)"""
        score = 10.0

        # Homopolymer penalty
        if re.search(r'([ATGC])\1{4,}', seq):
            score -= 4
        elif re.search(r'([ATGC])\1{3}', seq):
            score -= 2

        # 3' end G/C clamp (good)
        if seq[-1] in 'GC':
            score += 0.5

        # Poly-N penalty
        if re.search(r'N{2,}', seq):
            score -= 3

        # Repeat penalty
        if re.search(r'([ATGC])\1{2}$', seq):
            score -= 1.5

        return max(0, round(score, 1))

    def calculate_self_complementarity(self, seq):
        """Calculate self-complementarity using primer3"""
        try:
            homodimer_result = primer3.calc_homodimer(seq, mv_conc=50, dv_conc=1.5, dntp_conc=0.25, dna_conc=500)
            self_any_dg = homodimer_result.dg if hasattr(homodimer_result, 'dg') else 0.0
        except:
            self_any_dg = 0.0

        try:
            self_3p_result = primer3.calc_end_stability(seq, seq, mv_conc=50, dv_conc=1.5, dntp_conc=0.25, dna_conc=500)
            self_3p_dg = self_3p_result.dg if hasattr(self_3p_result, 'dg') else 0.0
        except:
            self_3p_dg = 0.0

        return {
            'self_any': round(abs(self_any_dg), 1),
            'self_3p': round(abs(self_3p_dg), 1)
        }

    def calculate_cross_complementarity(self, fwd_seq, rev_seq):
        """Calculate cross-complementarity between forward and reverse primers"""
        try:
            heterodimer_result = primer3.calc_heterodimer(fwd_seq, rev_seq, mv_conc=50, dv_conc=1.5, dntp_conc=0.25, dna_conc=500)
            cross_any_dg = heterodimer_result.dg if hasattr(heterodimer_result, 'dg') else 0.0
        except:
            cross_any_dg = 0.0

        try:
            cross_3p_result = primer3.calc_end_stability(fwd_seq, rev_seq, mv_conc=50, dv_conc=1.5, dntp_conc=0.25, dna_conc=500)
            cross_3p_dg = cross_3p_result.dg if hasattr(cross_3p_result, 'dg') else 0.0
        except:
            cross_3p_dg = 0.0

        return {
            'cross_any': round(abs(cross_any_dg), 1),
            'cross_3p': round(abs(cross_3p_dg), 1)
        }

    def penalty_score(self, fwd_seq, rev_seq):
        """Combined penalty (lower = better)"""
        fwd_tm = self.calculate_tm_basic(fwd_seq)
        rev_tm = self.calculate_tm_basic(rev_seq)
        fwd_gc = self.gc_content(fwd_seq)
        rev_gc = self.gc_content(rev_seq)

        tm_target = 60
        gc_target = 50

        tm_penalty = abs(fwd_tm - tm_target) * 0.05 + abs(rev_tm - tm_target) * 0.05
        gc_penalty = abs(fwd_gc - gc_target) * 0.03 + abs(rev_gc - gc_target) * 0.03
        qual_penalty = (10 - self.quality_score(fwd_seq)) * 0.1 + (10 - self.quality_score(rev_seq)) * 0.1

        return round(tm_penalty + gc_penalty + qual_penalty, 3)

    def scan_primer_pairs(self):
        """Complete primer scanning pipeline"""
        print("🔬 SCANNING HIGH-QUALITY PRIMER PAIRS...")

        candidates = []
        constraints = {
            'length': (18, 25),
            'gc': (40, 65),
            'tm': (55, 65),
            'product': (150, 400)
        }

        step_size = 30
        for start in range(50, len(self.sequence) - 450, step_size):
            for prod_size in range(150, 401, 50):
                end = start + prod_size
                if end >= len(self.sequence) - 50:
                    continue

                for length in range(18, 26):
                    fwd = self.sequence[start:start+length]
                    rev_seq_on_original_strand = self.sequence[end-length:end]
                    rev = str(Seq(rev_seq_on_original_strand).reverse_complement())

                    fwd_gc = self.gc_content(fwd)
                    rev_gc = self.gc_content(rev)
                    fwd_tm = self.calculate_tm_basic(fwd)
                    rev_tm = self.calculate_tm_basic(rev)

                    if (
                        constraints['gc'][0] <= fwd_gc <= constraints['gc'][1] and
                        constraints['gc'][0] <= rev_gc <= constraints['gc'][1] and
                        constraints['tm'][0] <= fwd_tm <= constraints['tm'][1] and
                        constraints['tm'][0] <= rev_tm <= constraints['tm'][1] and
                        self.quality_score(fwd) >= 6 and self.quality_score(rev) >= 6
                    ):

                        penalty = self.penalty_score(fwd, rev)

                        fwd_self = self.calculate_self_complementarity(fwd)
                        rev_self = self.calculate_self_complementarity(rev)
                        cross = self.calculate_cross_complementarity(fwd, rev)

                        candidates.append({
                            'forward': fwd,
                            'reverse': rev,
                            'Tm': f"{fwd_tm:.1f}/{rev_tm:.1f}",
                            'GC': f"{fwd_gc:.1f}/{rev_gc:.1f}",
                            'fwd_tm': fwd_tm,
                            'rev_tm': rev_tm,
                            'fwd_gc': fwd_gc,
                            'rev_gc': rev_gc,
                            'product_size': prod_size,
                            'penalty_score': penalty,
                            'tm_diff': abs(fwd_tm - rev_tm),
                            'fwd_start': start,
                            'fwd_end': start + length,
                            'rev_start': end - length,
                            'rev_end': end,
                            'fwd_self_any': fwd_self['self_any'],
                            'fwd_self_3p': fwd_self['self_3p'],
                            'rev_self_any': rev_self['self_any'],
                            'rev_self_3p': rev_self['self_3p'],
                            'cross_any': cross['cross_any'],
                            'cross_3p': cross['cross_3p'],
                            'fwd_quality': self.quality_score(fwd),
                            'rev_quality': self.quality_score(rev)
                        })

        candidates.sort(key=lambda x: x['penalty_score'])
        print(f"✅ Generated {len(candidates)} candidate pairs")
        self.candidates = candidates[:500]
        return self.candidates

    def generate_outputs(self, candidates):
        """Generate output files with requested columns"""

        # FIXED: Added fwd_tm, rev_tm, and fwd_gc to ensure plotting works
        columns_order = [
            'forward', 'reverse', 'Tm', 'GC', 'fwd_tm', 'rev_tm', 'fwd_gc', 'rev_gc',
            'product_size', 'penalty_score',
            'tm_diff', 'fwd_start', 'fwd_end', 'rev_start', 'rev_end',
            'fwd_self_any', 'fwd_self_3p', 'rev_self_any', 'rev_self_3p',
            'cross_any', 'cross_3p', 'fwd_quality', 'rev_quality'
        ]

        df = pd.DataFrame(candidates)
        df = df[columns_order]

        csv_path = OUTPUT_DIR / "primer_candidates.csv"
        df.to_csv(csv_path, index=False)
        print(f"💾 primer_candidates.csv → {csv_path}")

        # Generate text report
        report_path = OUTPUT_DIR / "top5_report.txt"
        with open(report_path, 'w', encoding='utf-8') as f:
            f.write("🏆 TOP 5 PRIMER PAIRS\n")
            f.write(f"Target: {self.seq_name}\n")
            f.write(f"Sequence Length: {len(self.sequence):,} bp\n")
            f.write("="*50 + "\n\n")

            for i, pair in enumerate(candidates[:5], 1):
                f.write(f"{i}. Penalty: {pair['penalty_score']:.3f}\n")
                f.write(f"   Forward:  {pair['forward']}\n")
                f.write(f"   Reverse:  {pair['reverse']}\n")
                f.write(f"   Tm: {pair['Tm']}°C  |  GC: {pair['GC']}%\n")
                f.write(f"   Product: {pair['product_size']} bp\n")
                f.write(f"   Tm Diff: {pair['tm_diff']}°C | Cross 3': {pair['cross_3p']}\n\n")

        print(f"📊 top5_report.txt → {report_path}")

        # Best primer info
        if candidates:
            best = candidates[0]
            print("\n🎯 BEST PRIMER:")
            print(f"FWD: {best['forward']}")
            print(f"REV: {best['reverse']}")
            print(f"Tm: {best['Tm']}°C | GC: {best['GC']}% ")
            print(f"Size: {best['product_size']} bp | Score: {best['penalty_score']}")
            print(f"Tm Diff: {best['tm_diff']}°C | Cross Any: {best['cross_any']}")

        # Generate graphs
        self.generate_all_graphs(df)

        return df

    def generate_all_graphs(self, df):
        """Generate all visualization graphs"""
        print("📈 GENERATING GRAPHS...")

        plt.style.use('seaborn-v0_8-whitegrid')

        # Graph 1: Penalty Score Distribution
        plt.figure(figsize=(8, 4))
        sns.histplot(df['penalty_score'], color='#3498db', kde=True)
        plt.title('Penalty Score Distribution', fontsize=12)
        plt.xlabel('Penalty Score')
        plt.ylabel('Count')
        plt.savefig(OUTPUT_DIR / 'graph1_penalty.png', dpi=150, bbox_inches='tight')
        plt.close()

        # Graph 2: Forward vs Reverse Tm Scatter
        plt.figure(figsize=(8, 4))
        scatter = plt.scatter(df['fwd_tm'], df['rev_tm'], c=df['penalty_score'],
                            cmap='viridis', alpha=0.7, s=30)
        plt.colorbar(scatter, label='Penalty Score')
        plt.title('Forward vs Reverse Tm', fontsize=12)
        plt.xlabel('Forward Tm (°C)')
        plt.ylabel('Reverse Tm (°C)')
        plt.savefig(OUTPUT_DIR / 'graph2_tm_scatter.png', dpi=150, bbox_inches='tight')
        plt.close()

        # Graph 3: GC Content Distribution
        plt.figure(figsize=(8, 4))
        sns.kdeplot(df['fwd_gc'], label='Forward GC', color='#e74c3c', linewidth=2)
        sns.kdeplot(df['rev_gc'], label='Reverse GC', color='#2ecc71', linewidth=2)
        plt.title('GC Content Distribution', fontsize=12)
        plt.xlabel('GC %')
        plt.ylabel('Density')
        plt.legend()
        plt.savefig(OUTPUT_DIR / 'graph3_gc_dist.png', dpi=150, bbox_inches='tight')
        plt.close()

        # Graph 4: Feature Correlation Heatmap
        plt.figure(figsize=(10, 8))
        corr_cols = ['fwd_tm', 'rev_tm', 'fwd_gc', 'rev_gc', 'penalty_score', 'cross_3p']
        sns.heatmap(df[corr_cols].corr(), annot=True, cmap='coolwarm',
                   center=0, fmt='.2f', square=True)
        plt.title('Feature Correlation Heatmap', fontsize=12)
        plt.savefig(OUTPUT_DIR / 'graph4_heatmap.png', dpi=150, bbox_inches='tight')
        plt.close()

        print("✅ All graphs generated!")


def main():
    """Main execution function"""
    print("🚀 PRIMER DESIGN PIPELINE v2.0")
    print("=" * 45)

    pipeline = CompletePrimerPipeline()

    if not pipeline.background_search_sequence():
        return

    candidates = pipeline.scan_primer_pairs()

    if candidates:
        pipeline.generate_outputs(candidates)
        print("\n✅ PIPELINE COMPLETED SUCCESSFULLY!")
    else:
        print("\n❌ No candidate primers found matching criteria.")


if __name__ == "__main__":
    main()
