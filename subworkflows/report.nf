// ─────────────────────────────────────────────────────────────────────────────
// subworkflows/report.nf
//
// GENERATE_REPORT  –  Read combined_results.parquet and produce a self-contained
// HTML report in the same dark-theme style as demo_report.html.
// No LLM interpretation — data, charts, and tables only.
// ─────────────────────────────────────────────────────────────────────────────

process GENERATE_REPORT {
    publishDir "${params.outdir}", mode: 'copy'
    container 'python:3.11-slim'
    conda     'conda-forge::python=3.11 conda-forge::pandas conda-forge::pyarrow'
    when: params.run_report && params.run_merge

    input:
    path parquet

    output:
    path 'annotation_report.html', emit: html

    script:
    """
    pip install -q pandas pyarrow

    python3 << 'PYEOF'
import pandas as pd
import json
import sys
import html as H
from collections import Counter
from datetime import datetime

PARQUET  = '${parquet}'
PALETTE  = ['#5b8dee','#e05c8a','#4ecdc4','#f7b731','#af52de','#ff6b6b','#a8e063','#fd9644',
            '#5b8dee','#e05c8a','#4ecdc4','#f7b731','#af52de','#ff6b6b','#a8e063']

# ── Load ──────────────────────────────────────────────────────────────────────
try:
    df = pd.read_parquet(PARQUET)
except Exception as e:
    print(f"Error reading parquet: {e}", file=sys.stderr)
    df = pd.DataFrame({'source_tool': pd.Series([], dtype='str')})

# ── Helpers ───────────────────────────────────────────────────────────────────
def tdf(tool):
    sub = df[df['source_tool'] == tool].copy()
    return sub.dropna(axis=1, how='all') if not sub.empty else sub

def gcol(d, *names):
    """Return first matching column as filled Series, or empty Series."""
    for n in names:
        if n in d.columns:
            return d[n].fillna('')
    return pd.Series([''] * len(d), index=d.index)

def top_counts(series, n=10):
    """Return (labels_list, values_list) for top-n non-empty values."""
    c = Counter(str(v) for v in series
                if str(v).strip() and str(v) not in ('', 'nan', 'None', 'NaN'))
    top = c.most_common(n)
    return [x[0] for x in top], [x[1] for x in top]

def esc(v):
    return H.escape(str(v)) if pd.notna(v) and str(v) not in ('', 'nan') else ''

def make_table(d, cols=None, max_rows=500):
    if d is None or d.empty:
        return '<p class="muted" style="padding:16px 0">No data available</p>'
    show = [c for c in (cols or list(d.columns)) if c in d.columns and c != 'source_tool']
    if not show:
        show = [c for c in d.columns if c != 'source_tool']
    d2 = d[show].head(max_rows)
    rows = ['<div class="table-wrap"><table><thead><tr>']
    for c in d2.columns:
        rows.append('<th>' + esc(c) + '</th>')
    rows.append('</tr></thead><tbody>')
    for _, row in d2.iterrows():
        rows.append('<tr>')
        for v in row:
            rows.append('<td>' + esc(v) + '</td>')
        rows.append('</tr>')
    rows.append('</tbody></table></div>')
    if len(d) > max_rows:
        rows.append(f'<p class="muted" style="font-size:11px;padding:4px 0 16px">'
                    f'Showing first {max_rows:,} of {len(d):,} rows</p>')
    return ''.join(rows)

def stat_card(val, label, detail='', cls=''):
    inner = (f'<div class="val {cls}">{val}</div>'
             f'<div class="label">{label}</div>'
             + (f'<div class="detail">{detail}</div>' if detail else ''))
    return f'<div class="card">{inner}</div>'

def section_open(sid, active=False):
    cls = 'section active' if active else 'section'
    return f'<section id="{sid}" class="{cls}">'

def two_col(*boxes):
    return '<div class="two-col">' + ''.join(boxes) + '</div>'

def chart_box(cid, title='', height=300):
    return (f'<div class="chart-box">'
            + (f'<h3 class="chart-title">{title}</h3>' if title else '')
            + f'<canvas id="{cid}" height="{height}"></canvas></div>')

def full_chart(cid, title='', height=340):
    return (f'<div class="full-chart">'
            + (f'<h3 class="chart-title">{title}</h3>' if title else '')
            + f'<canvas id="{cid}" height="{height}"></canvas></div>')

def jd(v):
    """Serialize to JSON string."""
    return json.dumps(v)

# ── Per-tool frames ───────────────────────────────────────────────────────────
amr  = tdf('amrfinder')
rgi  = tdf('card_rgi')
kof  = tdf('kofam')
pf   = tdf('plasmidfinder')
intf = tdf('integronfinder')
mob  = tdf('mob_suite')
ise  = tdf('isescan')

# ── Overview ──────────────────────────────────────────────────────────────────
OV_LABELS = ['AMRFinder','CARD RGI','KofamScan','PlasmidFinder','IntegronFinder','mob_suite','ISEScan']
OV_KEYS   = ['amrfinder','card_rgi','kofam','plasmidfinder','integronfinder','mob_suite','isescan']
OV_COUNTS = [int(df[df['source_tool'] == k].shape[0]) for k in OV_KEYS]

# ── AMRFinder ─────────────────────────────────────────────────────────────────
amr_class_l, amr_class_v = top_counts(gcol(amr, 'Class', 'class'), 12)
amr_sub_l,   amr_sub_v   = top_counts(gcol(amr, 'Subclass', 'subclass'), 12)
amr_gene_l,  amr_gene_v  = top_counts(gcol(amr, 'Gene symbol', 'gene_symbol'), 15)
amr_uniq_genes = len(set(str(v) for v in gcol(amr, 'Gene symbol', 'gene_symbol') if str(v).strip()))
amr_uniq_class = len(set(str(v) for v in gcol(amr, 'Class', 'class') if str(v).strip()))

# ── CARD RGI ──────────────────────────────────────────────────────────────────
rgi_mech_l, rgi_mech_v = top_counts(gcol(rgi, 'Resistance Mechanism', 'resistance_mechanism'), 8)
rgi_drug_l, rgi_drug_v = top_counts(gcol(rgi, 'Drug Class', 'drug_class'), 15)
rgi_aro_l,  rgi_aro_v  = top_counts(gcol(rgi, 'Best_Hit_ARO', 'best_hit_aro'), 15)
rgi_uniq_aro   = len(set(str(v) for v in gcol(rgi, 'Best_Hit_ARO') if str(v).strip()))
rgi_uniq_drug  = len(set(str(v) for v in gcol(rgi, 'Drug Class') if str(v).strip()))
rgi_uniq_mech  = len(set(str(v) for v in gcol(rgi, 'Resistance Mechanism') if str(v).strip()))

# ── KofamScan ─────────────────────────────────────────────────────────────────
kof_def_l, kof_def_v = top_counts(gcol(kof, 'KO_definition'), 15)
kof_ko_l,  kof_ko_v  = top_counts(gcol(kof, 'KO', 'ko'), 15)
try:
    kof_sig = int((pd.to_numeric(gcol(kof, 'E_value', 'e_value'), errors='coerce') < 1e-5).sum())
except Exception:
    kof_sig = 0
kof_uniq_ko = len(set(str(v) for v in gcol(kof, 'KO', 'ko') if str(v).strip()))

# ── ISEScan ───────────────────────────────────────────────────────────────────
ise_fam_l,  ise_fam_v  = top_counts(gcol(ise, 'family', 'Family'), 12)
ise_complete_n = int((gcol(ise, 'isComplete') == 'y').sum()) if not ise.empty else 0
ise_partial_n  = max(0, len(ise) - ise_complete_n)
ise_uniq_fam   = len(set(str(v) for v in gcol(ise, 'family', 'Family') if str(v).strip()))

# ── PlasmidFinder ─────────────────────────────────────────────────────────────
pf_match_l, pf_match_v = top_counts(gcol(pf, 'match_name', 'Name'), 12)

# ── IntegronFinder ────────────────────────────────────────────────────────────
intf_type_l, intf_type_v = top_counts(gcol(intf, 'annotation_type', 'type', 'Type', 'element_type'), 10)

# ── mob_suite ─────────────────────────────────────────────────────────────────
mob_rep_l,  mob_rep_v  = top_counts(gcol(mob, 'rep_type', 'rep_type(s)'), 10)
mob_mob_l,  mob_mob_v  = top_counts(gcol(mob, 'mob_type', 'MOB_type', 'relaxase_type(s)'), 10)

# ── Build HTML ────────────────────────────────────────────────────────────────
now = datetime.utcnow().strftime('%Y-%m-%d %H:%M UTC')
total_hits = len(df)
total_genomes_pf = len(set(gcol(pf, 'genome'))) if not pf.empty else '—'

CSS = '''
  :root {
    --bg: #0f1117; --surface: #1a1d27; --surface2: #22263a; --border: #2e3350;
    --accent: #5b8dee; --accent2: #e05c8a; --accent3: #4ecdc4; --accent4: #f7b731;
    --text: #e4e8f7; --muted: #7a82a8;
    --font: 'Segoe UI', system-ui, -apple-system, sans-serif;
  }
  * { box-sizing: border-box; margin: 0; padding: 0; }
  body { background: var(--bg); color: var(--text); font-family: var(--font);
         font-size: 14px; line-height: 1.6; }
  header {
    background: linear-gradient(135deg, #12172a 0%, #1a2040 60%, #0f1a30 100%);
    border-bottom: 1px solid var(--border); padding: 32px 40px 24px;
  }
  header h1 { font-size: 26px; font-weight: 700; color: #fff; letter-spacing: -0.3px; }
  header h1 span { color: var(--accent); }
  header p.sub { color: var(--muted); margin-top: 6px; font-size: 13px; }
  .tag-row { display: flex; gap: 8px; flex-wrap: wrap; margin-top: 14px; }
  .tag { background: var(--surface2); border: 1px solid var(--border); border-radius: 20px;
         padding: 3px 12px; font-size: 12px; color: var(--muted); }
  .tag.hi { border-color: var(--accent); color: var(--accent); }
  nav { background: var(--surface); border-bottom: 1px solid var(--border);
        display: flex; gap: 2px; padding: 0 40px; overflow-x: auto; }
  nav button {
    background: none; border: none; color: var(--muted); cursor: pointer;
    font-family: var(--font); font-size: 13px; padding: 14px 18px;
    border-bottom: 2px solid transparent; white-space: nowrap; transition: color .15s;
  }
  nav button:hover { color: var(--text); }
  nav button.active { color: var(--accent); border-bottom-color: var(--accent); }
  main { max-width: 1280px; margin: 0 auto; padding: 32px 40px; }
  .section { display: none; }
  .section.active { display: block; }
  h2 { font-size: 20px; font-weight: 600; margin-bottom: 6px; }
  .section-desc { color: var(--muted); font-size: 13px; margin-bottom: 24px; }
  h3.chart-title { font-size: 13px; font-weight: 600; color: var(--muted);
                   text-transform: uppercase; letter-spacing: .5px; margin-bottom: 14px; }
  h3.sub-heading { font-size: 15px; font-weight: 600; margin: 28px 0 12px; color: var(--text); }
  .card-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(160px,1fr));
               gap: 16px; margin-bottom: 28px; }
  .card { background: var(--surface); border: 1px solid var(--border); border-radius: 10px;
          padding: 20px; }
  .val { font-size: 28px; font-weight: 700; color: var(--accent); }
  .val.warn { color: var(--accent4); }
  .val.ok   { color: var(--accent3); }
  .val.pink { color: var(--accent2); }
  .label  { font-size: 12px; color: var(--muted); margin-top: 4px; }
  .detail { font-size: 12px; color: var(--text); margin-top: 2px; }
  .two-col { display: grid; grid-template-columns: 1fr 1fr; gap: 24px; margin-bottom: 28px; }
  .chart-box { background: var(--surface); border: 1px solid var(--border); border-radius: 10px;
               padding: 24px; }
  .chart-box canvas { max-height: 300px; }
  .full-chart { background: var(--surface); border: 1px solid var(--border); border-radius: 10px;
                padding: 24px; margin-bottom: 28px; }
  .full-chart canvas { max-height: 360px; }
  .table-wrap { overflow-x: auto; border-radius: 10px; border: 1px solid var(--border);
                margin-bottom: 28px; }
  table { width: 100%; border-collapse: collapse; font-size: 13px; }
  thead tr { background: var(--surface2); }
  th { padding: 11px 14px; text-align: left; font-weight: 600; color: var(--muted);
       white-space: nowrap; border-bottom: 1px solid var(--border); }
  td { padding: 9px 14px; border-bottom: 1px solid var(--border); color: var(--text); }
  tbody tr:last-child td { border-bottom: none; }
  tbody tr:hover { background: var(--surface2); }
  .muted { color: var(--muted); }
  .badge { display: inline-block; padding: 2px 9px; border-radius: 12px;
           font-size: 11px; font-weight: 600; }
  .badge-blue   { background: rgba(91,141,238,.15); color: #7daaff; }
  .badge-pink   { background: rgba(224,92,138,.15);  color: #f07aa8; }
  .badge-teal   { background: rgba(78,205,196,.15);  color: #5dd8d0; }
  .badge-amber  { background: rgba(247,183,49,.15);  color: #f7c94a; }
  .badge-purple { background: rgba(175,82,222,.15);  color: #c47ef5; }
  .divider { border: none; border-top: 1px solid var(--border); margin: 28px 0; }
  @media (max-width: 760px) {
    header, nav, main { padding-left: 16px; padding-right: 16px; }
    .two-col { grid-template-columns: 1fr; }
    .card-grid { grid-template-columns: repeat(2,1fr); }
  }
'''

# ─── Section builders ──────────────────────────────────────────────────────────

def build_overview():
    cards = (
        stat_card(total_hits, 'Total Annotations', 'across all tools')
      + stat_card(OV_COUNTS[0], 'AMRFinder Hits', '')
      + stat_card(OV_COUNTS[1], 'CARD RGI Hits', '')
      + stat_card(OV_COUNTS[2], 'KEGG Annotations', 'KofamScan')
      + stat_card(OV_COUNTS[3], 'Plasmid Replicons', 'PlasmidFinder', 'ok' if OV_COUNTS[3] == 0 else 'warn')
      + stat_card(OV_COUNTS[4], 'Integrons', 'IntegronFinder', 'ok' if OV_COUNTS[4] == 0 else 'warn')
      + stat_card(OV_COUNTS[5], 'Plasmid Contigs', 'mob_suite')
      + stat_card(OV_COUNTS[6], 'IS Elements', 'ISEScan', 'warn' if OV_COUNTS[6] > 0 else 'ok')
    )
    # Summary table of tool results
    rows = ''.join(
        f'<tr><td>{l}</td><td>{c:,}</td></tr>'
        for l, c in zip(OV_LABELS, OV_COUNTS)
    )
    summary_table = (
        '<div class="table-wrap"><table><thead><tr><th>Tool</th><th>Annotation Rows</th></tr></thead>'
        f'<tbody>{rows}</tbody></table></div>'
    )
    return (
        section_open('overview', active=True)
        + '<h2>Overview</h2>'
        + '<p class="section-desc">Summary of all annotation results across the pipeline.</p>'
        + '<div class="card-grid">' + cards + '</div>'
        + full_chart('ovBar', 'Annotation counts per tool', height=260)
        + '<h3 class="sub-heading">Tool summary</h3>'
        + summary_table
        + '</section>'
    )

def build_amr():
    n_strict = int((gcol(amr, 'Method', 'method') != '').sum()) if not amr.empty else 0
    cards = (
        stat_card(len(amr), 'Total Hits', '')
      + stat_card(amr_uniq_genes, 'Unique Genes', '')
      + stat_card(amr_uniq_class, 'Drug Classes', '')
    )
    amr_cols = ['Gene symbol', 'Sequence name', 'Class', 'Subclass', 'Method',
                '% Coverage of reference sequence', '% Identity to reference sequence']
    return (
        section_open('amr')
        + '<h2>AMRFinder</h2>'
        + '<p class="section-desc">AMR gene detection via NCBI AMRFinderPlus.</p>'
        + '<div class="card-grid">' + cards + '</div>'
        + two_col(chart_box('amrClassDonut', 'Drug class distribution'),
                  chart_box('amrSubBar',     'Top antibiotic subclasses'))
        + '<h3 class="sub-heading">Full results</h3>'
        + make_table(amr, amr_cols)
        + '</section>'
    )

def build_card():
    cards = (
        stat_card(len(rgi), 'Total Hits', '')
      + stat_card(rgi_uniq_aro, 'Unique AROs', '')
      + stat_card(rgi_uniq_drug, 'Drug Classes', '')
      + stat_card(rgi_uniq_mech, 'Mechanisms', '')
    )
    rgi_cols = ['Best_Hit_ARO', 'Drug Class', 'Resistance Mechanism', 'AMR Gene Family',
                'Best_Identities', 'Cut_Off', 'Contig']
    return (
        section_open('card')
        + '<h2>CARD RGI</h2>'
        + '<p class="section-desc">Comprehensive AMR annotation via the CARD Resistance Gene Identifier.</p>'
        + '<div class="card-grid">' + cards + '</div>'
        + two_col(chart_box('rgiMechDonut', 'Resistance mechanism'),
                  chart_box('rgiDrugBar',   'Top drug classes'))
        + '<h3 class="sub-heading">Full results</h3>'
        + make_table(rgi, rgi_cols)
        + '</section>'
    )

def build_kegg():
    cards = (
        stat_card(len(kof), 'Total Annotations', '')
      + stat_card(kof_uniq_ko, 'Unique KO Terms', '')
      + stat_card(kof_sig, 'E-value < 1e-5', '', 'ok')
    )
    kof_cols = ['gene_name', 'KO', 'score', 'E_value', 'KO_definition']
    return (
        section_open('kegg')
        + '<h2>KEGG / KofamScan</h2>'
        + '<p class="section-desc">KEGG Orthology annotation via KofamScan.</p>'
        + '<div class="card-grid">' + cards + '</div>'
        + full_chart('kofBar', 'Top KO definitions (by frequency)', height=340)
        + '<h3 class="sub-heading">Full results (sorted by score)</h3>'
        + make_table(
            kof.sort_values('score', ascending=False) if 'score' in kof.columns else kof,
            kof_cols)
        + '</section>'
    )

def build_plasmids():
    cards_pf = (
        stat_card(len(pf), 'PlasmidFinder Hits', '')
      + stat_card(len(set(str(v) for v in gcol(pf,'match_name') if str(v).strip())),
                  'Unique Replicons', '', 'ok' if len(pf) == 0 else '')
    )
    cards_mob = (
        stat_card(len(mob), 'mob_suite Contigs', '')
      + stat_card(len(set(str(v) for v in gcol(mob,'mob_type','MOB_type') if str(v).strip())),
                  'MOB Types', '')
    )
    pf_cols  = ['genome', 'contig', 'match_name', 'match_id', 'coverage', 'identity', 'hit_length']
    mob_cols = ['sample_id', 'contig_id', 'rep_type', 'mob_type', 'num_contigs',
                'total_length', 'gc', 'primary_cluster_id']
    return (
        section_open('plasmids')
        + '<h2>Plasmid Detection</h2>'
        + '<p class="section-desc">Plasmid replicons (PlasmidFinder) and plasmid reconstruction / MOB typing (mob_suite).</p>'
        + '<h3 class="sub-heading">PlasmidFinder</h3>'
        + '<div class="card-grid">' + cards_pf + '</div>'
        + two_col(chart_box('pfDonut', 'Replicon type distribution'),
                  chart_box('mobRepBar', 'mob_suite — replicon types'))
        + make_table(pf, pf_cols)
        + '<hr class="divider">'
        + '<h3 class="sub-heading">mob_suite</h3>'
        + '<div class="card-grid">' + cards_mob + '</div>'
        + make_table(mob, mob_cols)
        + '</section>'
    )

def build_integrons():
    cards = (
        stat_card(len(intf), 'Total Elements', '')
      + stat_card(len(intf_type_l), 'Element Types', '')
    )
    return (
        section_open('integrons')
        + '<h2>IntegronFinder</h2>'
        + '<p class="section-desc">Integron and gene cassette detection via IntegronFinder.</p>'
        + '<div class="card-grid">' + cards + '</div>'
        + full_chart('intfBar', 'Element type distribution', height=280)
        + '<h3 class="sub-heading">Full results</h3>'
        + make_table(intf)
        + '</section>'
    )

def build_is():
    cards = (
        stat_card(len(ise), 'Total IS Elements', '', 'warn' if len(ise) > 0 else 'ok')
      + stat_card(ise_uniq_fam, 'Families', '')
      + stat_card(ise_complete_n, 'Complete', '', 'ok')
      + stat_card(ise_partial_n, 'Partial', '', 'warn')
    )
    ise_cols = ['seqID', 'family', 'cluster', 'isBegin', 'isEnd', 'isLen',
                'ncopy4is', 'typeIS', 'isComplete', 'score', 'e-value']
    return (
        section_open('is')
        + '<h2>IS Elements — ISEScan</h2>'
        + '<p class="section-desc">Insertion sequence element detection and annotation via ISEScan.</p>'
        + '<div class="card-grid">' + cards + '</div>'
        + two_col(chart_box('iseFamDonut', 'IS family distribution'),
                  chart_box('iseTypeDonut', 'Complete vs partial'))
        + '<h3 class="sub-heading">Full results</h3>'
        + make_table(ise, ise_cols)
        + '</section>'
    )

# ── Chart.js initialisation script ────────────────────────────────────────────
CHART_JS = (
    '<script src="https://cdn.jsdelivr.net/npm/chart.js@4.4.0/dist/chart.umd.min.js"></script>'
    '<script>'
    'const PAL=' + jd(PALETTE) + ';'
    'Chart.defaults.color="#7a82a8";'
    'Chart.defaults.borderColor="#2e3350";'
    'Chart.defaults.font.family="Segoe UI,system-ui,sans-serif";'
    'Chart.defaults.font.size=12;'

    # data
    + 'const ovLabels=' + jd(OV_LABELS) + ';'
    + 'const ovCounts=' + jd(OV_COUNTS) + ';'
    + 'const amrClassL=' + jd(amr_class_l) + ';'
    + 'const amrClassV=' + jd(amr_class_v) + ';'
    + 'const amrSubL='   + jd(amr_sub_l)   + ';'
    + 'const amrSubV='   + jd(amr_sub_v)   + ';'
    + 'const rgiMechL='  + jd(rgi_mech_l)  + ';'
    + 'const rgiMechV='  + jd(rgi_mech_v)  + ';'
    + 'const rgiDrugL='  + jd(rgi_drug_l)  + ';'
    + 'const rgiDrugV='  + jd(rgi_drug_v)  + ';'
    + 'const kofDefL='   + jd(kof_def_l)   + ';'
    + 'const kofDefV='   + jd(kof_def_v)   + ';'
    + 'const pfMatchL='  + jd(pf_match_l)  + ';'
    + 'const pfMatchV='  + jd(pf_match_v)  + ';'
    + 'const mobRepL='   + jd(mob_rep_l)   + ';'
    + 'const mobRepV='   + jd(mob_rep_v)   + ';'
    + 'const intfTypeL=' + jd(intf_type_l) + ';'
    + 'const intfTypeV=' + jd(intf_type_v) + ';'
    + 'const iseFamL='   + jd(ise_fam_l)   + ';'
    + 'const iseFamV='   + jd(ise_fam_v)   + ';'
    + 'const iseComplN=' + jd(ise_complete_n) + ';'
    + 'const isePartN='  + jd(ise_partial_n)  + ';'

    # navigation
    + '''
const sections = document.querySelectorAll('.section');
const navBtns  = document.querySelectorAll('nav button');
const shown    = {};
function show(id, btn) {
  sections.forEach(s => { s.classList.remove('active'); });
  navBtns.forEach(b => b.classList.remove('active'));
  const sec = document.getElementById(id);
  if (sec) sec.classList.add('active');
  if (btn) btn.classList.add('active');
  if (!shown[id]) { shown[id]=true; initCharts(id); }
}
show('overview', document.querySelector('nav button'));
'''

    # chart helpers
    + '''
function hBar(id, labels, data, pal) {
  new Chart(document.getElementById(id), {
    type:'bar',
    data:{ labels:labels, datasets:[{ data:data, backgroundColor:pal||PAL,
           borderColor:'transparent', borderRadius:4 }] },
    options:{ responsive:true, maintainAspectRatio:false, indexAxis:'y',
      plugins:{ legend:{ display:false } },
      scales:{ x:{ grid:{ color:'#2e3350' } }, y:{ grid:{ display:false } } } }
  });
}
function vBar(id, labels, data, pal) {
  new Chart(document.getElementById(id), {
    type:'bar',
    data:{ labels:labels, datasets:[{ data:data, backgroundColor:pal||PAL,
           borderColor:'transparent', borderRadius:4 }] },
    options:{ responsive:true, maintainAspectRatio:false,
      plugins:{ legend:{ display:false } },
      scales:{ x:{ grid:{ display:false } }, y:{ grid:{ color:'#2e3350' } } } }
  });
}
function donut(id, labels, data, pal) {
  new Chart(document.getElementById(id), {
    type:'doughnut',
    data:{ labels:labels, datasets:[{ data:data, backgroundColor:pal||PAL,
           borderColor:'#1a1d27', borderWidth:2 }] },
    options:{ responsive:true, maintainAspectRatio:false,
      plugins:{ legend:{ position:'right', labels:{ boxWidth:12, padding:14 } } },
      cutout:'62%' }
  });
}
'''

    # lazy chart init per section
    + '''
function initCharts(id) {
  if (id==='overview') {
    hBar('ovBar', ovLabels, ovCounts);
  }
  if (id==='amr') {
    donut('amrClassDonut', amrClassL, amrClassV);
    hBar ('amrSubBar',     amrSubL,   amrSubV);
  }
  if (id==='card') {
    donut('rgiMechDonut', rgiMechL, rgiMechV);
    hBar ('rgiDrugBar',   rgiDrugL, rgiDrugV);
  }
  if (id==='kegg') {
    hBar('kofBar', kofDefL, kofDefV);
  }
  if (id==='plasmids') {
    donut('pfDonut',   pfMatchL, pfMatchV);
    hBar ('mobRepBar', mobRepL,  mobRepV);
  }
  if (id==='integrons') {
    vBar('intfBar', intfTypeL, intfTypeV);
  }
  if (id==='is') {
    donut('iseFamDonut',  iseFamL, iseFamV);
    donut('iseTypeDonut', ["Complete","Partial"], [iseComplN, isePartN],
          ["#4ecdc4","#f7b731"]);
  }
}
'''
    + '</script>'
)

# ── Assemble ──────────────────────────────────────────────────────────────────
NAV_ITEMS = [
    ('overview',  '&#127774; Overview'),
    ('amr',       '&#128138; AMRFinder'),
    ('card',      '&#127919; CARD RGI'),
    ('kegg',      '&#128202; KofamScan'),
    ('plasmids',  '&#128260; Plasmids'),
    ('integrons', '&#129522; Integrons'),
    ('is',        '&#128302; IS Elements'),
]

nav_html = '<nav>' + ''.join(
    f'<button onclick="show(\'{sid}\',this)">{label}</button>'
    for sid, label in NAV_ITEMS
) + '</nav>'

tag_items = [
    (f'{total_hits:,} total annotations', False),
    (f'{OV_COUNTS[0]} AMRFinder hits', False),
    (f'{OV_COUNTS[1]} CARD RGI hits', False),
    (f'{OV_COUNTS[6]} IS elements', OV_COUNTS[6] > 0),
    (f'{OV_COUNTS[3]} plasmid replicons', OV_COUNTS[3] > 0),
    (f'Generated {now}', False),
]
tags_html = '<div class="tag-row">' + ''.join(
    f'<span class="tag{\" hi\" if hi else \"\"}">{t}</span>'
    for t, hi in tag_items
) + '</div>'

html_parts = [
    '<!DOCTYPE html>',
    '<html lang="en">',
    '<head>',
    '<meta charset="UTF-8">',
    '<meta name="viewport" content="width=device-width, initial-scale=1.0">',
    '<title>Genome Annotation Report</title>',
    f'<style>{CSS}</style>',
    '</head>',
    '<body>',
    '<header>',
    '<h1>Genome Annotation Report — <span>nf-annotation</span></h1>',
    '<p class="sub">Automated multi-tool annotation pipeline</p>',
    tags_html,
    '</header>',
    nav_html,
    '<main>',
    build_overview(),
    build_amr(),
    build_card(),
    build_kegg(),
    build_plasmids(),
    build_integrons(),
    build_is(),
    '</main>',
    CHART_JS,
    '</body>',
    '</html>',
]

with open('annotation_report.html', 'w') as fh:
    fh.write('\\n'.join(html_parts))
print("Report written to annotation_report.html", file=sys.stderr)
PYEOF
    """
}

// ─────────────────────────────────────────────────────────────────────────────
// REPORT workflow
// ─────────────────────────────────────────────────────────────────────────────
workflow REPORT {

    take:
    parquet

    main:
    GENERATE_REPORT(parquet)

    emit:
    html = GENERATE_REPORT.out.html
}
