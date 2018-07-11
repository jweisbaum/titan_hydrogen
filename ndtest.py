import networkx as nx
import ndlib.models.epidemics.SIRModel as sir
import ndlib.models.ModelConfig as mc
from bokeh.io import output_notebook, show
from ndlib.viz.bokeh.DiffusionTrend import DiffusionTrend
from ndlib.viz.bokeh.DiffusionPrevalence import DiffusionPrevalence
from ndlib.viz.bokeh.MultiPlot import MultiPlot
import ndlib.models.epidemics.SISModel as sis
import ndlib.models.epidemics.SIModel as si
import ndlib.models.epidemics.ThresholdModel as th

# Network Definition
g = nx.erdos_renyi_graph(1000, 0.1)

# Model Selection
model = sir.SIRModel(g)

# Model Configuration
config = mc.Configuration()
config.add_model_parameter('beta', 0.001)
config.add_model_parameter('gamma', 0.01)
config.add_model_parameter("percentage_infected", 0.05)
model.set_initial_status(config)

# Simulation
iterations = model.iteration_bunch(200)
trends = model.build_trends(iterations)

viz = DiffusionTrend(model, trends)
p = viz.plot(width=400, height=400)
show(p)
