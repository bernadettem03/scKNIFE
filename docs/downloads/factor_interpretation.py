import scknife

fit = scknife.tutorial.load_fit()
programs = scknife.fp.group_factors(fit, scknife.tutorial.load_graph(), n_groups=4, bootstraps=200, seed=42)
programs.selected.head()

scknife.fp.plot_usage(programs)

scknife.fp.plot_group_selection(programs)

scknife.fp.plot_consensus(programs)

scknife.fp.plot_resolutions(programs)

scknife.fp.plot_node_dictionary(programs, by="factor")

scknife.fp.plot_node_dictionary(programs, by="group")

scknife.fp.plot_composition(programs)

programs.entity_kinds, programs.groups, programs.ngroups

# Other kinds: "gene", "GO_term", "reaction", "task", "cellchat_entity"
scknife.fp.plot_entities(programs, kind="geneset")

scknife.fp.plot_group_activity(programs, group="NG")

scknife.fp.plot_group_activity(programs, group="CC")

scknife.fp.plot_associations(programs)

scknife.fp.plot_cell_profiles(programs)

scknife.fp.plot_coverage(programs)

scknife.fp.plot_constituents(programs)
