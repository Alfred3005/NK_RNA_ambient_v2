"""
Fase 3: GSEA Preranked por Subtipo NK — Diseño Tricéfalo
=========================================================
Tres análisis independientes, todos usando la MATRIZ COMPLETA de genes
(no solo los significativos) como ranked list preranked:

  A. GLOBAL (nk_cell_general) — visión de conjunto del envejecimiento NK
  B. NK CD56dim               — subtipo maduro/citotóxico
  C. NK CD56bright             — subtipo inmunomodulador/IFN

Para cada análisis se aplican:
  1. GSEA Prerank (Wald stat) — Gold standard DESeq2
  2. GSEA Prerank (Métrica combinada: sign(LFC) × -log10(padj)) — captura
     cambios sutiles y coordinados (recomendado en el manual de referencia)
  3. ORA con g:Profiler + background personalizado (todos los genes analizados)

Gene sets: MSigDB Hallmark (H), KEGG 2021, Reactome 2022, GO:BP 2023

Referencia metodológica:
  "Análisis de Enriquecimiento de Vías Genéticas" — scAR_python_validation_v2_clean
  Script previo: scAR_python_validation_v4_clean/scripts/20_functional_enrichment.py
"""

import os
import numpy as np
import pandas as pd
import gseapy as gp
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import seaborn as sns
import warnings
warnings.filterwarnings("ignore")

# ─────────────────────────────────────────────────────────────────────────────
# CONFIGURACIÓN GLOBAL
# ─────────────────────────────────────────────────────────────────────────────
BASE_RESULTS = '/mnt/c/Users/PREDATOR/Documents/Antigravity_workspaces/NK_pipeline_RNA_ambient/scAR_python_validation_v4_clean_subtypes_abundance/results/subtypes'
GSEA_DIR     = os.path.join(BASE_RESULTS, 'gsea')
os.makedirs(GSEA_DIR, exist_ok=True)

# Tres análisis: (clave_archivo, nombre_display, subdirectorio_output, ruta_absoluta_opcional)
# Para el análisis Global se usa el dataset completo V4-clean (todos los NK, sin split por subtipo)
GLOBAL_DESEQ_PATH = '/mnt/c/Users/PREDATOR/Documents/Antigravity_workspaces/NK_pipeline_RNA_ambient/scAR_python_validation_v4_clean/results/pydeseq2/deseq2_results_v4_final.csv'

ANALYSES = [
    ('global',       'NK Global (Todos los NK · V4-Clean)', 'global',     GLOBAL_DESEQ_PATH),
    ('nk_cd56dim',   'NK CD56dim',                          'cd56dim',    None),
    ('nk_cd56bright','NK CD56bright',                        'cd56bright', None),
]

# Gene sets a usar en GSEA (disponibles offline vía gseapy)
GENE_SETS_GSEA = [
    'MSigDB_Hallmark_2020',
    'KEGG_2021_Human',
    'Reactome_2022',
    'GO_Biological_Process_2023',
]

# Gene sets a usar en ORA vía g:Profiler
ORA_SOURCES = ['GO:BP', 'KEGG', 'REAC', 'HP']

# Parámetros GSEA
MIN_SIZE     = 15
MAX_SIZE     = 500
PERMUTATIONS = 1000
FDR_CUT      = 0.25    # Umbral estándar GSEA
SEED         = 42

# ─────────────────────────────────────────────────────────────────────────────
# FUNCIONES DE RANKING
# ─────────────────────────────────────────────────────────────────────────────
def build_rankings(df: pd.DataFrame) -> dict:
    """
    Construye dos rankings a partir del DataFrame completo de DESeq2:
      - wald_stat:       estadístico de Wald (gold standard DESeq2)
      - combined_metric: sign(LFC) × -log10(padj) [Manual de referencia, pág. 15]
    Retorna un dict con ambas Series ordenadas descendentemente.
    """
    df = df.copy()
    df = df[~df.index.isna() & ~df.index.duplicated(keep='first')]

    # 1. Wald stat
    wald = df['stat'].dropna().replace([np.inf, -np.inf], np.nan).dropna()
    wald = wald.sort_values(ascending=False)

    # 2. Métrica combinada
    eps = 1e-300
    df['padj_clean'] = df['padj'].fillna(1.0).clip(lower=eps)
    df['lfc_clean']  = df['log2FoldChange'].fillna(0.0)
    df['combined']   = np.sign(df['lfc_clean']) * (-np.log10(df['padj_clean']))
    combined = df['combined'].dropna().replace([np.inf, -np.inf], np.nan).dropna()
    combined = combined.sort_values(ascending=False)

    return {'wald_stat': wald, 'combined_metric': combined}

# ─────────────────────────────────────────────────────────────────────────────
# FUNCIONES GSEA PRERANK
# ─────────────────────────────────────────────────────────────────────────────
def _extract_col(df: pd.DataFrame, patterns: list[str]) -> pd.Series | None:
    """Devuelve la primera columna cuyo nombre (lowercase) coincida con algún patrón."""
    for col in df.columns:
        lc = col.lower().strip()
        for pat in patterns:
            if pat in lc:
                val = df[col]
                if isinstance(val, pd.DataFrame):
                    val = val.iloc[:, 0]
                return val.copy()
    return None


def run_gsea_prerank(ranked: pd.Series, gene_set: str, out_dir: str,
                     metric_label: str) -> pd.DataFrame | None:
    """Ejecuta gseapy.prerank y retorna un DataFrame limpio con columnas estándar."""
    os.makedirs(out_dir, exist_ok=True)
    try:
        res = gp.prerank(
            rnk=ranked,
            gene_sets=gene_set,
            outdir=out_dir,
            min_size=MIN_SIZE,
            max_size=MAX_SIZE,
            permutation_num=PERMUTATIONS,
            ascending=False,
            no_plot=True,
            verbose=False,
            seed=SEED,
            threads=4,
        )
        if res is None or res.res2d is None or len(res.res2d) == 0:
            return None
        raw = res.res2d.copy()

        # In gseapy 1.x, Term is typically a column. We prioritize the column, then fallback to index.
        term_series = _extract_col(raw, ['term'])
        if term_series is not None:
            term_vals = term_series.astype(str).values
        else:
            name_series = _extract_col(raw, ['name'])
            if name_series is not None and not name_series.astype(str).str.lower().str.strip().isin(['prerank', '']).all():
                term_vals = name_series.astype(str).values
            else:
                term_vals = raw.index.astype(str).values

        nes_col  = _extract_col(raw, ['nes'])
        fdr_col  = _extract_col(raw, ['fdr', 'q val', 'qval'])
        pval_col = _extract_col(raw, ['nom p', 'nom_p', 'p-val'])

        clean = pd.DataFrame({
            'Term': term_vals,
            'NES':  pd.to_numeric(nes_col,  errors='coerce').values if nes_col  is not None else np.nan,
            'FDR':  pd.to_numeric(fdr_col,  errors='coerce').values if fdr_col  is not None else np.nan,
            'pval': pd.to_numeric(pval_col, errors='coerce').values if pval_col is not None else np.nan,
        })
        clean['metric']      = metric_label
        clean['gene_set_db'] = gene_set
        return clean

    except Exception as e:
        print(f"    ❌ Error GSEA [{gene_set}|{metric_label}]: {e}")
        return None




def gsea_for_analysis(ranked_dict: dict, analysis_key: str, display_name: str):
    """
    Corre GSEA Prerank con ambas métricas para todas las colecciones.
    Si ya existe un CSV de resultados para esa combinación, lo carga en lugar de re-ejecutar.
    Retorna dict: gene_set -> df_results_merged
    """
    out_base = os.path.join(GSEA_DIR, analysis_key)
    all_res: dict[str, list] = {}

    for gs in GENE_SETS_GSEA:
        # ── Ruta del CSV de resultados ya consolidado ──────────────────────────
        merged_csv = os.path.join(out_base, f'gsea_{gs}.csv')
        if os.path.exists(merged_csv):
            print(f"    ♻️  Cargando resultados previos: {gs}")
            merged = pd.read_csv(merged_csv)
            # Asegurar columna Term como string Serie
            if 'Term' in merged.columns:
                if isinstance(merged['Term'], pd.DataFrame):
                    merged['Term'] = merged['Term'].iloc[:, 0]
                merged['Term'] = merged['Term'].astype(str)
            all_res[gs] = merged
            if 'FDR' in merged.columns:
                sig = merged[merged['FDR'].astype(float) < FDR_CUT]
                print(f"      ✅ Significativos (FDR<{FDR_CUT}): {len(sig)}")
            continue

        frames = []
        for metric_label, ranked in ranked_dict.items():
            gs_dir = os.path.join(out_base, gs, metric_label)
            print(f"    🔬 {gs} [{metric_label}]...")
            df_res = run_gsea_prerank(ranked, gs, gs_dir, metric_label)
            if df_res is not None:
                frames.append(df_res)

        if frames:
            merged = pd.concat(frames, ignore_index=True)
            # Asegurar Term como string Serie antes de guardar
            if 'Term' in merged.columns:
                if isinstance(merged['Term'], pd.DataFrame):
                    merged['Term'] = merged['Term'].iloc[:, 0]
                merged['Term'] = merged['Term'].astype(str)
            all_res[gs] = merged
            merged.to_csv(merged_csv, index=False)

            if 'FDR' in merged.columns:
                sig = merged[merged['FDR'].astype(float) < FDR_CUT]
                sig.to_csv(os.path.join(out_base, f'gsea_sig_{gs}.csv'), index=False)
                print(f"      ✅ Significativos (FDR<{FDR_CUT}): {len(sig)}")

    return all_res

# ─────────────────────────────────────────────────────────────────────────────
# FUNCIÓN ORA con g:Profiler
# ─────────────────────────────────────────────────────────────────────────────
def run_ora(df_full: pd.DataFrame, analysis_key: str, display_name: str):
    """ORA con g:Profiler usando background = todos los genes analizados."""
    try:
        from gprofiler import GProfiler
    except ImportError:
        print("    ⚠️  gprofiler-official no instalado — saltando ORA.")
        return None

    out_dir = os.path.join(GSEA_DIR, analysis_key, 'ora')
    os.makedirs(out_dir, exist_ok=True)

    background_genes = df_full.dropna(subset=['padj']).index.tolist()
    sig_genes = df_full[(df_full['padj'] < 0.05)].index.tolist()

    print(f"    📋 ORA g:Profiler — Query: {len(sig_genes)} DEGs | Background: {len(background_genes)} genes")

    if len(sig_genes) < 5:
        print("    ⚠️  Muy pocos genes significativos para ORA.")
        return None

    gpc = GProfiler(return_dataframe=True)
    try:
        ora_res = gpc.profile(
            organism='hsapiens',
            query=sig_genes,
            no_evidences=False,
            background=background_genes,
            sources=ORA_SOURCES,
        )
    except Exception as e:
        print(f"    ❌ Error ORA g:Profiler: {e}")
        return None

    if ora_res is None or ora_res.empty:
        print("    ⚠️  Sin términos ORA significativos.")
        return None

    ora_res.to_csv(os.path.join(out_dir, 'ora_enrichment_results.csv'), index=False)
    sig_ora = ora_res[ora_res['p_value'] < 0.05]
    sig_ora.to_csv(os.path.join(out_dir, 'ora_sig_results.csv'), index=False)
    print(f"      ✅ ORA completado. {len(sig_ora)} términos significativos (p<0.05).")
    return ora_res

# ─────────────────────────────────────────────────────────────────────────────
# VISUALIZACIONES
# ─────────────────────────────────────────────────────────────────────────────
def plot_gsea_dotplot(df_res: pd.DataFrame, display_name: str, gs_name: str,
                      output_path: str, top_n: int = 25):
    """
    Dotplot combinado de ambas métricas (Wald + Combined):
      - Eje Y: términos (top N por |NES| del Wald stat)
      - Color: NES divergente
      - Tamaño: -log10(FDR)
      - Forma: metric (Wald = círculo, Combined = diamante)
    """
    if 'FDR' not in df_res.columns or 'NES' not in df_res.columns:
        return

    df_plot = df_res.copy()
    df_plot['NES'] = df_plot['NES'].astype(float)
    df_plot['FDR'] = df_plot['FDR'].astype(float).clip(lower=1e-10)
    df_plot['neg_log_fdr'] = -np.log10(df_plot['FDR'])
    df_plot['abs_NES'] = df_plot['NES'].abs()

    if 'Term' not in df_plot.columns:
        df_plot['Term'] = df_plot.index.astype(str)
    else:
        # Garantizar que Term sea una Serie de strings (no un DataFrame por columnas dup)
        if isinstance(df_plot['Term'], pd.DataFrame):
            df_plot['Term'] = df_plot['Term'].iloc[:, 0]
        df_plot['Term'] = df_plot['Term'].astype(str)

    df_plot['Term'] = df_plot['Term'].str.replace('_', ' ').str.title().str[:65]

    # Seleccionar top terms por Wald stat
    wald_df = df_plot[df_plot['metric'] == 'wald_stat']
    if wald_df.empty:
        wald_df = df_plot

    top_terms = wald_df.sort_values('abs_NES', ascending=False).head(top_n)['Term'].tolist()
    # Deduplicar manteniendo orden (un pathway puede aparecer en ambas métricas)
    seen = set()
    top_terms_unique = []
    for t in top_terms:
        if t not in seen:
            seen.add(t)
            top_terms_unique.append(t)
    top_terms = top_terms_unique

    df_plot = df_plot[df_plot['Term'].isin(top_terms)]

    if df_plot.empty:
        return

    # Orden por NES (Wald), deduplicado
    order_raw = (wald_df[wald_df['Term'].isin(top_terms)]
                 .sort_values('NES', ascending=True)['Term'].tolist())
    seen2 = set()
    order = []
    for t in order_raw:
        if t not in seen2:
            seen2.add(t)
            order.append(t)

    if not order:
        return

    df_plot = df_plot.copy()
    df_plot['Term'] = pd.Categorical(df_plot['Term'], categories=order, ordered=True)
    df_plot = df_plot.sort_values('Term')


    n_terms = len(order)
    fig_h = max(8, n_terms * 0.45 + 3)
    fig, ax = plt.subplots(figsize=(13, fig_h))

    # Mapa term → posición Y
    term_pos = {t: i for i, t in enumerate(order)}

    markers = {'wald_stat': 'o', 'combined_metric': 'D'}
    vmax = max(abs(df_plot['NES'].max()), abs(df_plot['NES'].min())) + 0.1
    norm = plt.Normalize(vmin=-vmax, vmax=vmax)

    for metric, mdf in df_plot.groupby('metric', observed=True):
        # Operar fila a fila — evita problemas de .loc con índices duplicados
        mdf = mdf.copy().reset_index(drop=True)
        mdf['y_pos'] = mdf['Term'].map(term_pos)
        mdf = mdf.dropna(subset=['y_pos'])
        if mdf.empty:
            continue
        mdf['y_pos'] = mdf['y_pos'].astype(int)

        ax.scatter(
            mdf['NES'].values,
            mdf['y_pos'].values,
            c=mdf['NES'].values,
            cmap='RdBu_r', norm=norm,
            s=(mdf['neg_log_fdr'].values * 45).clip(30, 700),
            marker=markers.get(metric, 'o'),
            alpha=0.82, edgecolors='white', linewidths=0.6,
            label='Wald stat' if metric == 'wald_stat' else 'Sign × -log₁₀(padj)',
            zorder=3,
        )

    ax.axvline(0, color='#94a3b8', ls='--', lw=1, alpha=0.7, zorder=1)
    ax.set_yticks(range(n_terms))
    ax.set_yticklabels(order, fontsize=9)
    ax.set_xlabel('NES (Normalized Enrichment Score)', fontsize=11)
    ax.set_title(
        f'GSEA Preranked — {display_name}\n{gs_name.replace("_", " ")} · FDR<{FDR_CUT} · Top {n_terms}',
        fontsize=12, weight='bold', pad=15
    )
    sm = plt.cm.ScalarMappable(cmap='RdBu_r', norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, pad=0.01, shrink=0.5)
    cbar.set_label('NES', fontsize=10)
    ax.legend(title='Métrica', fontsize=9, loc='lower right', framealpha=0.85)
    ax.set_facecolor('#f8fafc')
    ax.grid(axis='x', color='white', lw=1.5, zorder=2)
    ax.spines[['top','right']].set_visible(False)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"    📊 Dotplot guardado: {os.path.basename(output_path)}")


def plot_ora_dotplot(ora_res: pd.DataFrame, display_name: str, output_path: str, top_n: int = 15):
    """Dotplot ORA estilo g:Profiler."""
    if ora_res is None or ora_res.empty:
        return

    sig = ora_res[ora_res['p_value'] < 0.05].copy()
    if sig.empty:
        return

    sig['gene_ratio'] = sig['intersection_size'] / sig['term_size'].clip(lower=1)
    sig['neg_log_p']  = -np.log10(sig['p_value'].clip(lower=1e-300))
    sig['short_name'] = sig['name'].apply(lambda x: (x[:55] + '…') if len(x) > 55 else x)

    plot_df = sig.sort_values('p_value').head(top_n)

    fig, ax = plt.subplots(figsize=(11, max(5, len(plot_df) * 0.55 + 2)))
    scatter = ax.scatter(
        plot_df['gene_ratio'],
        range(len(plot_df)),
        s=plot_df['intersection_size'] * 25,
        c=plot_df['neg_log_p'],
        cmap='viridis',
        alpha=0.85,
        edgecolors='black',
        linewidths=0.5,
    )
    ax.set_yticks(range(len(plot_df)))
    ax.set_yticklabels(plot_df['short_name'], fontsize=9)
    ax.set_xlabel('Gene Ratio (Solapamiento / Tamaño Vía)', fontsize=11)
    ax.set_title(
        f'ORA g:Profiler — {display_name}\nGO:BP / KEGG / Reactome · p<0.05 · Top {len(plot_df)}',
        fontsize=12, weight='bold', pad=15
    )
    cbar = fig.colorbar(scatter, ax=ax, shrink=0.5)
    cbar.set_label('-log₁₀(p-value)', fontsize=10)

    sizes_leg = [5, 10, 20, 50]
    for sz in sizes_leg:
        ax.scatter([], [], s=sz*25, c='gray', alpha=0.5, label=str(sz), edgecolors='black')
    ax.legend(title='Genes encontrados', fontsize=8, loc='lower right', framealpha=0.85)
    ax.set_facecolor('#f8fafc')
    ax.spines[['top','right']].set_visible(False)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"    📊 Dotplot ORA guardado: {os.path.basename(output_path)}")


def plot_comparative_barplot(all_sig: dict, output_path: str):
    """
    Barplot comparativo: número de pathways significativos
    por análisis y colección de gene sets.
    """
    records = []
    for (key, display, _), gs_dict in all_sig.items():
        for gs, df_sig in gs_dict.items():
            if df_sig is not None and not df_sig.empty:
                if 'FDR' in df_sig.columns:
                    n_up   = len(df_sig[(df_sig['FDR'].astype(float) < FDR_CUT) & (df_sig['NES'].astype(float) > 0)])
                    n_down = len(df_sig[(df_sig['FDR'].astype(float) < FDR_CUT) & (df_sig['NES'].astype(float) < 0)])
                    records.append({'Análisis': display, 'Gene Set': gs.replace('_', ' '), 'Up': n_up, 'Down': -n_down})

    if not records:
        return

    rec_df = pd.DataFrame(records)
    analyses = rec_df['Análisis'].unique()
    gss = rec_df['Gene Set'].unique()

    fig, axes = plt.subplots(1, len(analyses), figsize=(5 * len(analyses), 6), sharey=False)
    if len(analyses) == 1:
        axes = [axes]

    colors_up   = '#ef4444'
    colors_down = '#3b82f6'

    for ax, analysis in zip(axes, analyses):
        sub = rec_df[rec_df['Análisis'] == analysis]
        y   = range(len(sub))
        ax.barh([g for g in sub['Gene Set']], sub['Up'],   color=colors_up,   alpha=0.85, label='Activados')
        ax.barh([g for g in sub['Gene Set']], sub['Down'], color=colors_down, alpha=0.85, label='Reprimidos')
        ax.axvline(0, color='#374151', lw=1)
        ax.set_title(analysis, fontsize=11, weight='bold')
        ax.set_xlabel('Pathways (FDR<0.25)', fontsize=10)
        ax.spines[['top','right']].set_visible(False)
        if ax == axes[0]:
            ax.legend(fontsize=9)

    fig.suptitle('Resumen GSEA: Pathways Significativos por Subtipo NK y Base de Datos',
                 fontsize=13, weight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  📊 Barplot comparativo guardado: {os.path.basename(output_path)}")


# ─────────────────────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────────────────────
def main():
    print("\n" + "="*70)
    print("FASE 3: GSEA TRICÉFALO — GLOBAL / CD56dim / CD56bright")
    print("Referencia: Manual de Enriquecimiento de Vías + 20_functional_enrichment.py")
    print("="*70)

    all_sig = {}   # Para el barplot comparativo final

    for (file_key, display_name, subdir, custom_path) in ANALYSES:
        print(f"\n{'─'*65}")
        print(f"▶ ANÁLISIS: {display_name}")
        print(f"{'─'*65}")

        # Ruta del archivo DESeq2
        if custom_path:
            deseq_path = custom_path
        else:
            deseq_path = os.path.join(BASE_RESULTS, f'deseq2_results_{file_key}.csv')

        if not os.path.exists(deseq_path):
            print(f"  ⚠️  Archivo no encontrado: {deseq_path}. Omitiendo.")
            continue

        df = pd.read_csv(deseq_path, index_col=0)
        df = df[~df.index.isna() & ~df.index.duplicated(keep='first')]
        print(f"  Genes en tabla completa: {len(df)} | Significativos (padj<0.05): {len(df[df['padj']<0.05]) if 'padj' in df.columns else 'N/A'}")

        # Guardar ranked lists de referencia
        ranked_dict = build_rankings(df)
        for metric_label, ranked in ranked_dict.items():
            rnk_path = os.path.join(GSEA_DIR, subdir, f'ranked_{metric_label}.rnk')
            os.makedirs(os.path.dirname(rnk_path), exist_ok=True)
            ranked.reset_index().to_csv(rnk_path, sep='\t', index=False, header=False)
        print(f"  Ranked lists: Wald={len(ranked_dict['wald_stat'])} genes | Combined={len(ranked_dict['combined_metric'])} genes")

        # ── GSEA Prerank ──────────────────────────────────────────────────────
        print(f"\n  [GSEA Prerank]")
        out_key  = os.path.join(GSEA_DIR, subdir)
        gsea_res = gsea_for_analysis(ranked_dict, subdir, display_name)

        # Generar dotplots por colección
        sig_by_gs = {}
        for gs, df_res in gsea_res.items():
            if 'FDR' not in df_res.columns:
                continue
            df_sig = df_res[df_res['FDR'].astype(float) < FDR_CUT]
            sig_by_gs[gs] = df_sig
            if len(df_sig) > 0:
                dotplot_path = os.path.join(out_key, f'dotplot_{gs}.png')
                plot_gsea_dotplot(df_res, display_name, gs, dotplot_path, top_n=25)

        all_sig[file_key, display_name, subdir] = sig_by_gs

        # ── ORA g:Profiler ────────────────────────────────────────────────────
        print(f"\n  [ORA g:Profiler]")
        ora_res = run_ora(df, subdir, display_name)
        if ora_res is not None:
            ora_plot_path = os.path.join(out_key, 'ora_dotplot.png')
            plot_ora_dotplot(ora_res, display_name, ora_plot_path, top_n=20)

    # ── Barplot comparativo final ─────────────────────────────────────────────
    print(f"\n{'─'*65}")
    print("▶ RESUMEN COMPARATIVO")
    print(f"{'─'*65}")
    comp_path = os.path.join(GSEA_DIR, 'comparative_summary_barplot.png')
    plot_comparative_barplot(all_sig, comp_path)

    print(f"\n{'='*70}")
    print(f"✅ FASE 3 COMPLETADA — Resultados en: {GSEA_DIR}")
    print(f"{'='*70}\n")


if __name__ == "__main__":
    main()
