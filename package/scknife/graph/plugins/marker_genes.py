import os
import pandas as pd
import networkx as nx
from scknife.utils.helpers import add_edge_with_source, normalize_celltype


canonical_markers = {
    "Epithelial": ["EPCAM","KRT8","KRT18","CLDN4","DSC2","ELF3","SFN","ITGB4","ITGB6","KRT6A","KRT7","LAD1","LAMA3","LAMB3","SLPI","ST14","TGFA","SCEL","PTGES","SH2D3A","AP1M2","MPZL2","EHF","G0S2","LSR","S100A14"],
    "Endothelial": ["PECAM1","ACKR1","ACVRL1","ANXA3","ART4","BMX","ERG","FLT1","GJA4","GPR4","HTR2B","KDR","KRT19","TACSTD2","PDE3A","RALA","MAPK12","SELE","SNTB2","TIE1","CLDN5","VWF","HYAL2","EIF2B2","TNFSF18","GDF3","IGF2BP3","ARHGEF15","MMRN1","CD93","CLEC1A","EMCN","SOX18","ROBO4","RASIP1","CXorf36","MMRN2","FAM124B","MYCT1","KANK3"],
    "B cells": ["CD79B","BLK","CD19","MS4A1","CD22","CD79A","PNOC","SNX2","MBD4","STAG3","PWP1","SP140","GGA2","STAP1","FCRL2","SMC6"],
    "T cells": ["CD3D","CD3E","CD3G"],
    "NK Cells": ["FASLG","CD247","GZMB","IL2RB","LIM2","PRF1","PTGDR","NCR1","NMUR1","GNLY","CD160","TBX21"],
    "CAFs": ["PDGFRB","PDGFRA","PDPN","DCN","COL1A1","COL1A2","FAP","THY1","CXCL12","MCAM","ACTA2","MYLK","MYL9","TAGLN","CALD1","CAV1","MEF2C","NOTCH3","RGS5"],
    "Plasma cells": ["JCHAIN","IGKC","IGHG1","IGHG2","IGHG3","IGHG4","AMPD1","TNFRSF17","CCNC","CD79A","FKBP2","GOLGA3","HSF4","KCNN3","MGAT2","PNOC","RGS13","SRP54","SSR4","HSP90B1","UBE2G1","RNF113A","MANF","RNGTT","EBAG9","PREB","MRPS31","SEC24A","FNDC3A","ATF6","MAST1","UBXN4","TBL2","SERP1","ALG5","GMPPA","YIPF1","LAX1","TMEM39A","LTB4R2","SPATS2","UBA5","ADM2","ZBP1","USP48","MTDH","IMP4"],
    "Myeloid cells": ["BMX","CA4","CEACAM3","FCGR3B","FPR2","CXCR1","CXCR2","MEFV","AATK","P2RY13","TREM1","VNN3","MMP25","BTNL8","TRPM6","CD1B","CD1E","CD86","FCER2","CCL17","ALDH1A2","RRP1B","KCNK13","AZU1","MS4A3","CLC","CPA3","CRHBP","CTSG","ELANE","MS4A2","FLT3","HDC","MPO","SERPINB10","PRG2","PRTN3","RNASE2","RNASE3","EPX","NAALADL1","HPGDS","ACP2","ATOX1","ATP6V0C","ATP6V1E1","C1QA","CD9","CD48","CLCN7","CCR1","CRYBB1","CYBB","FCER1G","FDX1","FOLR2","FPR3","FTL","HEXA","HEXB","KCNJ5","LAIR1","M6PR","MDH1","NUBP1","P2RX7","RAC1","CCL8","CCL22","STX4","TCEB1","TYROBP","VAMP8","MARCO","CD84","ATP6V0E1","ATP6V1F","CD163","CIAO1","LY86","ARHGEF11","HS3ST2","ARPC4","ATP6AP2","BAIAP2","PRDX3","ERP29","LILRB5","SDS","LILRB4","VSIG4","FKBP15","ZZZ3","SIGLEC7","SIGLEC9","ADAMDEC1","ORMDL2","COMMD9","SPG21","MS4A4A","UBE2D4","TRAPPC2L","TREM2","TMEM70","ARL8B","TMEM126B","SLAMF8","C12orf4","HAMP","GUF1","YIF1B","AIF1","APAF1","RHOA","RHOG","ARNT","ASGR2","C3AR1","CASP5","TNFRSF8","CD33","CSF1R","CSF3R","FCAR","FCN1","HCK","HK3","HRH2","HSPA6","KCNMB1","MAP3K3","MNDA","MYO1F","CFP","PHKG2","PTGIR","RARA","S100A12","UBE2D1","UPK3A","VASP","BEST1","LST1","NUP214","IQGAP1","ATP6V0D1","CD101","EIF4E2","CALCOCO2","LILRB2","TGOLN2","CAMKK2","LILRB1","CD300C","FGL2","LILRA1","LILRB3","LILRA2","TREX1","GABARAP","RPH3A","ACAP2","CLEC5A","IL17RA","FBXL5","CLEC4E","VENTX","PILRA","METTL9","TLR7","TLR8","RIN2","RHOT1","TMEM127","TMEM9B","RETN","RPGRIP1","CLEC7A","ISG15","IFITM1","IFITM2","IFITM3","LYZ","S100A8","S100A9","FCN1","LTF","BPI","DEFA1","DEFA3","DEFA4","DEFA5","DEFA6","CAMP"],
    "Epithelial Luminal Mature": ["FOXA1","MLPH","SPINK1","ESR1","PGR","GATA3","KRT18","KRT19","CDH1","TFF3"],
    "Epithelial Luminal Progenitor": ["LALBA","SLC34A2","CSN3","KIT","ELF5","S100A8","FOLR1","CLDN1"],
    "Epithelial cycling": ["MKI67","TOP2A","CDK1","AURKA","AURKB","CCNB1","CCNB2","CCNA2","BUB1","BUB1B","NEK2","PBK","KIF11","KIF20A","UBE2C","CDC20","KIF2C","CEP55","NDC80","DLGAP5","TPX2","PRC1","ASPM","CENPE","CENPF","KIF23","NUSAP1","HMMR","TACC3","ZWINT","RRM2","TYMS","PCNA","MCM2","MCM3","MCM5","MCM6","MCM7","CDC6","CDCA8","CDKN3"],
    "Basal_Myoepithelial": ["KRT5","KRT14","KRT17","OXTR","ACTA2","TAGLN","CNN1","MYLK","DST","MYL9"],
    "CD8 T cells": ["CD8A","CD8B","CD27","CTSW","CX3CR1","GZMH","GZMK","DNAJB1","KLRB1","LAIR2","LY9","PTGDR","PTPN4","RBL2","FBXW4","TSPAN32","KLRG1","CD160","RWDD3","IPCEF1","SIRPG","NAA16","FAM134C"],
    "CD4 T cells": ["CD2","CD4","CD28","CCR7","CCR8","ITK","TRAT1"],
    "T-Regulatory cells": ["CD5","CD28","CCR4","CCR8","CTLA4","IL2RA","IL10RA","FOXP3","ICOS"],
}

tnbc_biomarkers = [
    "KIF20A","DEDD","NCAPD2","PLK2","PPM1A","Cyclin E","FLI-1","ATIP3",
    "Telomerase","ErbB-2","ADAM12","PRMT5","HIC1","HOXB7","eIF4A","CPSF6",
    "TMPO-AS1","Circ-HER2","Circ_0000977","CircCAPG","PSMB5","miR-34a",
    "miR-138","miR-21","miR-211-5p","miR-101","miR-340","miR-877-5p",
    "miR-3613-3p","Hsp90","DNAJB4","RhoGDI-2","PON2","FSCN1","PDLIM2",
    "GRB7","CLU","PRKD3","EGFR","MET","PYK2","Glutathione","CHKA",
    "p70S6K","YBX1","PKM2","FZD7","NF-kB","Wnt/beta-catenin","LCN2","CD90",
    "CBP","FOXM1","B7-H4","LAG3","PD-1","BARD1","RAD51","PALB2",
    "Integrin","5-HT7 receptor","Sigma-2 receptor","EphA10","CDCP1",
    "PDGFRbeta","EP300","TENM4","EPHB6","MCT-1","IL-6","IL-6R","BCL6",
]


def load_marker_genes(celltypist_xlsx_path=None):
    try:
        import decoupler as dc
    except ImportError as exc:
        raise ImportError("Marker graph rebuilding requires pip install 'scknife[graphs]'; the frozen graph tutorial does not need it") from exc
    print("  [markers] Loading 4-source marker combination with NORMALIZATION...")

    dfs_to_combine = []

    print("    [1/4] PanglaoDB...")
    m = dc.get_resource(name="PanglaoDB", organism="human", license="academic")
    m = m[m["human"].astype(bool) & m["canonical_marker"].astype(bool) & (m["human_sensitivity"].astype(float) > 0.5)]
    print(f"      -> {len(m)} markers from PanglaoDB")
    dfs_to_combine.append(m[["genesymbol", "cell_type"]])

    print("    [2/4] Canonical markers...")
    can_df = pd.DataFrame([
        {"genesymbol": g, "cell_type": ct}
        for ct, gl in canonical_markers.items() for g in gl
    ])
    print(f"      -> {len(can_df)} markers from canonical set")
    dfs_to_combine.append(can_df)

    print("    [3/4] TNBC biomarkers...")
    tnbc_df = pd.DataFrame({"genesymbol": tnbc_biomarkers, "cell_type": "TNBC_marker"})
    print(f"      -> {len(tnbc_df)} TNBC biomarkers")
    dfs_to_combine.append(tnbc_df)

    print("    [4/4] CellTypist encyclopedia...")
    if celltypist_xlsx_path and os.path.exists(celltypist_xlsx_path):
        try:
            ct_xlsx = pd.read_excel(celltypist_xlsx_path)
            ct_list = []
            for _, row in ct_xlsx.iterrows():
                cell_type = str(row["High-hierarchy cell types"]) if pd.notna(row["High-hierarchy cell types"]) else None
                if not cell_type:
                    continue
                all_genes = []
                curated_str = str(row.get("Curated markers", "")) if pd.notna(row.get("Curated markers")) else ""
                if curated_str and curated_str != "nan":
                    genes = [g.strip() for g in str(curated_str).split(",") if g.strip()]
                    all_genes.extend(genes)
                top10_str = str(row.get("Top 10 important genes from the CellTypist model", "")) if pd.notna(row.get("Top 10 important genes from the CellTypist model")) else ""
                if top10_str and top10_str != "nan":
                    genes = [g.strip() for g in str(top10_str).split(",") if g.strip()]
                    all_genes.extend(genes)
                for gene in all_genes:
                    if gene:
                        ct_list.append({"genesymbol": gene, "cell_type": cell_type})
            if ct_list:
                ct_df = pd.DataFrame(ct_list).drop_duplicates()
                print(f"      -> {len(ct_df)} markers from CellTypist")
                dfs_to_combine.append(ct_df)
            else:
                print("      No markers parsed from CellTypist")
        except Exception as e:
            print(f"      Could not load CellTypist: {str(e)[:60]}")
    else:
        print(f"      CellTypist file not found or not provided")

    combined = pd.concat(dfs_to_combine, ignore_index=True).drop_duplicates()
    print(f"    Total (after dedup): {len(combined)} gene-celltype pairs")

    print("\n  [markers] Applying cell type normalization...")
    combined["cell_type"] = combined["cell_type"].apply(normalize_celltype)
    combined = combined.dropna(subset=["cell_type"])

    G = nx.Graph()
    for _, row in combined.iterrows():
        g = str(row["genesymbol"])
        c = str(row["cell_type"])
        if not G.has_node(g):
            G.add_node(g, kind="gene")
        if not G.has_node(c):
            G.add_node(c, kind="cell_type")
        add_edge_with_source(G, g, c, source="marker_graph", type="gene-celltype")

    n_genes = sum(1 for _, a in G.nodes(data=True) if a.get("kind") == "gene")
    n_ct = sum(1 for _, a in G.nodes(data=True) if a.get("kind") == "cell_type")
    print(f"\n    Graph built: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
    print(f"      (genes={n_genes}, cell_types={n_ct})")
    return {"graph": G, "node_kinds": {"gene": n_genes, "cell_type": n_ct}}
