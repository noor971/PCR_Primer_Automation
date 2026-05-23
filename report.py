#!/usr/bin/env python3
"""
Professional HTML Report Generator
Generates beautiful HTML reports for primer design results
"""

import os
import pandas as pd
from pathlib import Path
from datetime import datetime

# Paths
INPUT_DIR = Path(r"/content")
OUTPUT_DIR = Path(r"/content/MINI_PROJECT")
CSV_FILE = OUTPUT_DIR / "primer_candidates.csv"
HTML_OUTPUT = OUTPUT_DIR / "report.html"


def read_csv_data():
    """Read primer candidates CSV"""
    if not CSV_FILE.exists():
        print(f"❌ CSV file not found: {CSV_FILE}")
        return None
    
    df = pd.read_csv(CSV_FILE)
    print(f"✅ Loaded {len(df)} primer candidates")
    return df


def generate_html_table(df, top_n=20):
    """Generate HTML table for top primers"""
    html = '<table class=" primers-table">\n'
    html += '  <thead>\n    <tr>\n'
    headers = ['#', 'Forward', 'Reverse', 'Tm', 'GC', 'Product', 'Penalty', 'Quality']
    for h in headers:
        html += f'      <th>{h}</th>\n'
    html += '    </tr>\n  </thead>\n'
    
    html += '  <tbody>\n'
    for i, row in df.head(top_n).iterrows():
        html += '    <tr>\n'
        html += f'      <td>{i+1}</td>\n'
        html += f'      <td class="seq"><code>{row["forward"]}</code></td>\n'
        html += f'      <td class="seq"><code>{row["reverse"]}</code></td>\n'
        html += f'      <td>{row["Tm"]}°C</td>\n'
        html += f'      <td>{row["GC"]}%</td>\n'
        html += f'      <td>{row["product_size"]} bp</td>\n'
        html += f'      <td class="penalty">{row["penalty_score"]:.3f}</td>\n'
        fwd_q = row.get('fwd_quality', 0)
        rev_q = row.get('rev_quality', 0)
        avg_q = (fwd_q + rev_q) / 2
        html += f'      <td><span class="score">{avg_q:.1f}</span></td>\n'
        html += '    </tr>\n'
    html += '  </tbody>\n</table>'
    
    return html


def generate_statistics(df):
    """Generate statistics section"""
    stats = {
        'Total Primers': len(df),
        'Avg Penalty': f"{df['penalty_score'].mean():.3f}",
        'Best Penalty': f"{df['penalty_score'].min():.3f}",
        'Avg Tm (Fwd)': f"{df['fwd_tm'].mean():.1f}°C",
        'Avg Tm (Rev)': f"{df['rev_tm'].mean():.1f}°C",
        'Avg GC (Fwd)': f"{df['fwd_gc'].mean():.1f}%",
        'Avg GC (Rev)': f"{df['rev_gc'].mean():.1f}%",
        'Avg Product': f"{df['product_size'].mean():.0f} bp",
    }
    
    html = '<div class="stats-grid">\n'
    for key, value in stats.items():
        html += f'  <div class="stat-card">\n'
        html += f'    <div class="stat-value">{value}</div>\n'
        html += f'    <div class="stat-label">{key}</div>\n'
        html += '  </div>\n'
    html += '</div>'
    
    return html


def generate_charts_section():
    """Generate charts section HTML"""
    graphs = [
        ('graph1_penalty.png', 'Penalty Score Distribution'),
        ('graph2_tm_scatter.png', 'Forward vs Reverse Tm'),
        ('graph3_gc_dist.png', 'GC Content Distribution'),
        ('graph4_heatmap.png', 'Feature Correlation Heatmap'),
    ]
    
    html = '<div class="charts-grid">\n'
    for graph, title in graphs:
        graph_path = OUTPUT_DIR / graph
        html += f'  <div class="chart-card">\n'
        html += f'    <h3>{title}</h3>\n'
        if graph_path.exists():
            html += f'    <img src="{graph}" alt="{title}" class="chart-img">\n'
        else:
            html += f'    <div class="chart-placeholder">Chart Not Generated</div>\n'
        html += '  </div>\n'
    html += '</div>'
    
    return html


def generate_html_report(df):
    """Generate complete HTML report"""
    
    top5 = df.head(5)
    
    html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Primer Design Report | Mahnoor</title>
    <style>
        * {{
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }}
        
        body {{
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            background: linear-gradient(135deg, #1a1a2e 0%, #16213e 100%);
            min-height: 100vh;
            padding: 20px;
            color: #eee;
        }}
        
        .container {{
            max-width: 1400px;
            margin: 0 auto;
        }}
        
        /* Header */
        .header {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            padding: 40px;
            border-radius: 20px;
            margin-bottom: 30px;
            text-align: center;
            box-shadow: 0 20px 60px rgba(0,0,0,0.3);
        }}
        
        .header h1 {{
            font-size: 2.5em;
            margin-bottom: 10px;
            text-shadow: 2px 2px 4px rgba(0,0,0,0.3);
        }}
        
        .header .subtitle {{
            font-size: 1.2em;
            opacity: 0.9;
        }}
        
        .header .meta {{
            margin-top: 20px;
            font-size: 0.9em;
            opacity: 0.8;
        }}
        
        /* Sections */
        .section {{
            background: rgba(255,255,255,0.05);
            backdrop-filter: blur(10px);
            border-radius: 20px;
            padding: 30px;
            margin-bottom: 30px;
            border: 1px solid rgba(255,255,255,0.1);
        }}
        
        .section h2 {{
            font-size: 1.8em;
            margin-bottom: 20px;
            color: #667eea;
            border-bottom: 2px solid #667eea;
            padding-bottom: 10px;
        }}
        
        /* Stats Grid */
        .stats-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 20px;
            margin-top: 20px;
        }}
        
        .stat-card {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            padding: 25px;
            border-radius: 15px;
            text-align: center;
            transition: transform 0.3s;
        }}
        
        .stat-card:hover {{
            transform: translateY(-5px);
        }}
        
        .stat-value {{
            font-size: 2em;
            font-weight: bold;
            margin-bottom: 5px;
        }}
        
        .stat-label {{
            font-size: 0.9em;
            opacity: 0.8;
            text-transform: uppercase;
            letter-spacing: 1px;
        }}
        
        /* Table */
        .primers-table {{
            width: 100%;
            border-collapse: collapse;
            margin-top: 20px;
            font-size: 0.9em;
        }}
        
        .primers-table thead {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        }}
        
        .primers-table th {{
            padding: 15px;
            text-align: left;
            font-weight: 600;
            text-transform: uppercase;
            letter-spacing: 1px;
        }}
        
        .primers-table td {{
            padding: 12px 15px;
            border-bottom: 1px solid rgba(255,255,255,0.1);
        }}
        
        .primers-table tbody tr:hover {{
            background: rgba(255,255,255,0.05);
        }}
        
        .primers-table code {{
            background: rgba(102, 126, 234, 0.3);
            padding: 3px 8px;
            border-radius: 4px;
            font-family: 'Courier New', monospace;
            color: #a8d8ff;
        }}
        
        .primers-table .seq {{
            font-family: 'Courier New', monospace;
            font-size: 0.85em;
            max-width: 150px;
            overflow: hidden;
            text-overflow: ellipsis;
        }}
        
        .primers-table .penalty {{
            color: #00ff88;
            font-weight: bold;
        }}
        
        .primers-table .score {{
            background: linear-gradient(135deg, #00ff88 0%, #667eea 100%);
            padding: 4px 10px;
            border-radius: 20px;
            font-weight: bold;
        }}
        
        /* Charts */
        .charts-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(400px, 1fr));
            gap: 25px;
            margin-top: 20px;
        }}
        
        .chart-card {{
            background: rgba(255,255,255,0.03);
            padding: 20px;
            border-radius: 15px;
            border: 1px solid rgba(255,255,255,0.1);
        }}
        
        .chart-card h3 {{
            margin-bottom: 15px;
            color: #a8d8ff;
        }}
        
        .chart-img {{
            width: 100%;
            border-radius: 10px;
            box-shadow: 0 10px 30px rgba(0,0,0,0.3);
        }}
        
        .chart-placeholder {{
            height: 250px;
            display: flex;
            align-items: center;
            justify-content: center;
            background: rgba(255,255,255,0.05);
            border-radius: 10px;
            color: #666;
        }}
        
        /* Best Primer Section */
        .best-primer {{
            background: linear-gradient(135deg, #00ff88 0%, #667eea 100%);
            padding: 30px;
            border-radius: 20px;
            margin-bottom: 30px;
        }}
        
        .best-primer h2 {{
            color: #fff;
            border-bottom-color: rgba(255,255,255,0.3);
        }}
        
        .primer-pair {{
            display: grid;
            grid-template-columns: 1fr 1fr;
            gap: 20px;
            margin-top: 20px;
        }}
        
        .primer-box {{
            background: rgba(0,0,0,0.3);
            padding: 20px;
            border-radius: 15px;
        }}
        
        .primer-box h4 {{
            color: #00ff88;
            margin-bottom: 10px;
            font-size: 1.1em;
        }}
        
        .primer-box .sequence {{
            font-family: 'Courier New', monospace;
            font-size: 1.3em;
            word-break: break-all;
            color: #a8d8ff;
            background: rgba(255,255,255,0.1);
            padding: 15px;
            border-radius: 10px;
        }}
        
        .primer-box .details {{
            margin-top: 15px;
            display: grid;
            grid-template-columns: repeat(3, 1fr);
            gap: 10px;
        }}
        
        .detail-item {{
            text-align: center;
        }}
        
        .detail-item .label {{
            font-size: 0.8em;
            opacity: 0.7;
            text-transform: uppercase;
        }}
        
        .detail-item .value {{
            font-size: 1.2em;
            font-weight: bold;
        }}
        
        /* Footer */
        .footer {{
            text-align: center;
            padding: 30px;
            color: #666;
            font-size: 0.9em;
        }}
        
        .footer a {{
            color: #667eea;
            text-decoration: none;
        }}
        
        /* Responsive */
        @media (max-width: 768px) {{
            .header h1 {{
                font-size: 1.8em;
            }}
            
            .charts-grid {{
                grid-template-columns: 1fr;
            }}
            
            .primer-pair {{
                grid-template-columns: 1fr;
            }}
            
            .primers-table {{
                font-size: 0.8em;
            }}
            
            .primers-table .seq {{
                max-width: 100px;
            }}
        }}
    </style>
</head>
<body>
    <div class="container">
        <!-- Header -->
        <div class="header">
            <h1>🧬 Primer Design Report</h1>
            <div class="subtitle">Complete Analysis Results</div>
            <div class="meta">
                Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')} | 
                Project: Primer Helper Design | 
                Author: Mahnoor
            </div>
        </div>
        
        <!-- Statistics Section -->
        <div class="section">
            <h2>📊 Summary Statistics</h2>
            {generate_statistics(df)}
        </div>
        
        <!-- Best Primer Section -->
        <div class="section best-primer">
            <h2>🏆 Best Primer Pair</h2>
            <div class="primer-pair">
                <div class="primer-box">
                    <h4>Forward Primer</h4>
                    <div class="sequence">{top5.iloc[0]['forward']}</div>
                    <div class="details">
                        <div class="detail-item">
                            <div class="label">Tm</div>
                            <div class="value">{top5.iloc[0]['fwd_tm']}°C</div>
                        </div>
                        <div class="detail-item">
                            <div class="label">GC</div>
                            <div class="value">{top5.iloc[0]['fwd_gc']}%</div>
                        </div>
                        <div class="detail-item">
                            <div class="label">Quality</div>
                            <div class="value">{top5.iloc[0].get('fwd_quality', 0):.1f}</div>
                        </div>
                    </div>
                </div>
                <div class="primer-box">
                    <h4>Reverse Primer</h4>
                    <div class="sequence">{top5.iloc[0]['reverse']}</div>
                    <div class="details">
                        <div class="detail-item">
                            <div class="label">Tm</div>
                            <div class="value">{top5.iloc[0]['rev_tm']}°C</div>
                        </div>
                        <div class="detail-item">
                            <div class="label">GC</div>
                            <div class="value">{top5.iloc[0]['rev_gc']}%</div>
                        </div>
                        <div class="detail-item">
                            <div class="label">Quality</div>
                            <div class="value">{top5.iloc[0].get('rev_quality', 0):.1f}</div>
                        </div>
                    </div>
                </div>
            </div>
        </div>
        
        <!-- Top Primers Table -->
        <div class="section">
            <h2>📋 Top 20 Primer Pairs</h2>
            {generate_html_table(df, top_n=20)}
        </div>
        
        <!-- Charts Section -->
        <div class="section">
            <h2>📈 Visualization Analysis</h2>
            {generate_charts_section()}
        </div>
        
        <!-- Footer -->
        <div class="footer">
            <p>🔗 <a href="https://github.com/noor971/primer_helper_design">GitHub Repository</a></p>
            <p>Made with ❤️ by Mahnoor</p>
        </div>
    </div>
</body>
</html>"""
    
    return html_content


def main():
    """Main function to generate report"""
    print("🎨 GENERATING HTML REPORT...")
    print("=" * 45)
    
    # Read data
    df = read_csv_data()
    if df is None:
        return
    
    # Generate HTML
    html_content = generate_html_report(df)
    
    # Save to file
    with open(HTML_OUTPUT, 'w', encoding='utf-8') as f:
        f.write(html_content)
    
    print(f"✅ Report generated: {HTML_OUTPUT}")
    print(f"📄 File size: {HTML_OUTPUT.stat().st_size / 1024:.1f} KB")


if __name__ == "__main__":
    main()
