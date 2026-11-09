import scknife

fit = scknife.tutorial.load_fit()
means = scknife.optional.state_programs(fit, scknife.tutorial.load_activities(fit), scknife.tutorial.load_graph(), n=8)
scknife.optional.plot_states(means)
