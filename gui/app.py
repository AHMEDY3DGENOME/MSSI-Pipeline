import streamlit as st
import tempfile
import os
from pathlib import Path
from MSSI_Pipeline.loader     import MSSILoader
from MSSI_Pipeline.calculator import MSSICalculator
from MSSI_Pipeline.classifier import MSSIClassifier
from MSSI_Pipeline.visualizer import MSSIVisualizer
from MSSI_Pipeline.reporter   import MSSIReporter


STRESS_COLORS = {
    "very_high": "#E24B4A",
    "high":      "#EF9F27",
    "moderate":  "#378ADD",
    "low":       "#1D9E75",
}


def setup_page():
    st.set_page_config(
        page_title="MSSI Pipeline",
        page_icon="🧬",
        layout="wide",
    )
    st.title("🧬 MSSI Pipeline")
    st.caption("Multi-HSP Stress Signature Index — Spodoptera frugiperda")
    st.divider()


def sidebar_options(genes: list) -> dict:
    st.sidebar.header("Settings")
    st.sidebar.subheader("Gene Weights")
    weights = {}
    for gene in genes:
        weights[gene] = st.sidebar.slider(
            gene, min_value=0.1, max_value=3.0,
            value=1.0, step=0.1
        )
    st.sidebar.divider()
    show_plots  = st.sidebar.checkbox("Show Plots",   value=True)
    export_pdf  = st.sidebar.checkbox("Export PDF",   value=True)
    export_json = st.sidebar.checkbox("Export JSON",  value=True)
    return {
        "weights":     weights,
        "show_plots":  show_plots,
        "export_pdf":  export_pdf,
        "export_json": export_json,
    }


def show_metrics(results):
    st.subheader("MSSI Scores")
    cols = st.columns(len(results))
    for col, r in zip(cols, results):
        color = STRESS_COLORS[r.stress_level]
        col.markdown(
            f"""
            <div style='text-align:center; padding:16px;
                        border-radius:10px; border: 2px solid {color};'>
                <h3 style='color:{color}; margin:0'>{r.mssi_score:.3f}</h3>
                <p style='margin:4px 0; font-weight:600'>{r.location}</p>
                <p style='margin:0; font-size:12px; color:gray'>
                    {r.stress_label}</p>
                <p style='margin:0; font-size:12px'>
                    Confidence: {int(r.confidence*100)}%</p>
            </div>
            """,
            unsafe_allow_html=True
        )


def show_summary_table(results):
    st.subheader("Summary Table")
    import pandas as pd
    rows = [
        {
            "Location":    r.location,
            "MSSI Score":  r.mssi_score,
            "Normalized":  r.normalized_score,
            "Stress Level": r.stress_label,
            "Confidence":  f"{int(r.confidence*100)}%",
        }
        for r in results
    ]
    st.dataframe(pd.DataFrame(rows), use_container_width=True)


def show_gene_table(results):
    st.subheader("Gene Expression Details")
    import pandas as pd
    genes = list(results[0].gene_profile.keys())
    rows  = []
    for r in results:
        row = {"Location": r.location}
        for gene in genes:
            row[f"{gene} FC"] = r.gene_profile[gene]["fold_change"]
            row[f"{gene} SE"] = r.gene_profile[gene]["se"]
        rows.append(row)
    st.dataframe(pd.DataFrame(rows), use_container_width=True)


def show_plots(visualizer: MSSIVisualizer):
    st.subheader("Visualizations")
    tab1, tab2, tab3, tab4 = st.tabs([
        "MSSI Bar", "Heatmap", "Gene Profiles", "Confidence"
    ])
    with tab1:
        st.pyplot(visualizer.plot_mssi_bar(save=False))
    with tab2:
        st.pyplot(visualizer.plot_heatmap(save=False))
    with tab3:
        st.pyplot(visualizer.plot_gene_profiles(save=False))
    with tab4:
        st.pyplot(visualizer.plot_confidence(save=False))


def show_downloads(reporter: MSSIReporter, options: dict):
    st.subheader("Download Results")
    col1, col2 = st.columns(2)
    if options["export_pdf"]:
        pdf_path = reporter.export_pdf()
        with open(pdf_path, "rb") as f:
            col1.download_button(
                "Download PDF Report",
                data=f.read(),
                file_name="mssi_report.pdf",
                mime="application/pdf"
            )
    if options["export_json"]:
        json_path = reporter.export_json()
        with open(json_path, "rb") as f:
            col2.download_button(
                "Download JSON Results",
                data=f.read(),
                file_name="mssi_results.json",
                mime="application/json"
            )


def main():
    setup_page()

    uploaded = st.file_uploader(
        "Upload Excel File (.xlsx)",
        type=["xlsx", "xls"]
    )

    if uploaded is None:
        st.info("Please upload your Excel file to start analysis.")
        return

    with tempfile.NamedTemporaryFile(delete=False, suffix=".xlsx") as tmp:
        tmp.write(uploaded.read())
        tmp_path = tmp.name

    try:
        loader  = MSSILoader(tmp_path)
        data    = loader.load()
        genes   = loader.get_genes()
        options = sidebar_options(genes)

        with st.spinner("Running MSSI analysis..."):
            calculator = MSSICalculator(data, weights=options["weights"])
            stats      = calculator.get_full_stats()
            normalized = calculator.normalize_scores()
            classifier = MSSIClassifier(stats, normalized)
            results    = classifier.classify_all()

        show_metrics(results)
        st.divider()
        show_summary_table(results)
        st.divider()
        show_gene_table(results)

        if options["show_plots"]:
            st.divider()
            with tempfile.TemporaryDirectory() as tmpdir:
                visualizer = MSSIVisualizer(results, output_dir=tmpdir)
                reporter   = MSSIReporter(results,   output_dir=tmpdir)
                show_plots(visualizer)
                st.divider()
                show_downloads(reporter, options)

    except Exception as e:
        st.error(f"Error: {e}")
    finally:
        os.unlink(tmp_path)


if __name__ == "__main__":
    main()