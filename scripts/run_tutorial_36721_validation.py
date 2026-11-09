import scknife

errors = scknife.optional.rank_projection(scknife.tutorial.load_fit(), scknife.tutorial.load_graph(), ranks=(8,16,24))
scknife.optional.plot_projection(errors)
