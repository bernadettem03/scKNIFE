from scknife.graph.plugins.breast_ppi import load_breast_ppi
from scknife.graph.plugins.metabolic_tasks import load_metabolic_tasks
from scknife.graph.plugins.rxn_by_gene import load_rxn_by_gene
from scknife.graph.plugins.task_by_rxn import load_task_by_rxn
from scknife.graph.plugins.cellchat import load_cellchat
from scknife.graph.plugins.marker_genes import load_marker_genes
from scknife.graph.plugins.cell_state import load_cell_state
from scknife.graph.plugins.hallmark import load_hallmark
from scknife.graph.plugins.cytopus import load_cytopus_processes, load_cytopus_identities
from scknife.graph.plugins.enrichr import (
    load_enrichr_go,
    load_enrichr_reactome,
    load_enrichr_kegg,
    load_enrichr_trrust,
    load_enrichr_archs4,
    load_enrichr_tabula,
)
from scknife.graph.plugins.human_gem import load_human_gem
