import numpy as np
import matplotlib.pyplot as plt
import pyvcell.vcml as vc


biomodel = vc.load_vcml_url("https://raw.githubusercontent.com/virtualcell/pyvcell/refs/heads/main/examples/models/Tutorial_MultiApp_PDE.vcml")
model = biomodel.model
print(biomodel)
print(model.parameter_values)

model.set_parameter_value(name="r0.Kf", value=20.0)
print(model.parameter_values)


# make mesh smaller to make it run faster - was (100, 100, 36)
sim = biomodel.applications[0].simulations[0]
sim.mesh_size = (50, 50, 18)

# take N samples from log-normal distibution for r0.Kf
N=5
r0_Kf_values = np.random.lognormal(mean=0, sigma=1.0, size=N)

# run N simulations and store results
all_results = []
for val in r0_Kf_values:
    model.set_parameter_value("r0.Kf", val)
    print(f"running sim with r0.Kf={val}")
    all_results.append(vc.simulate(biomodel, sim.name))

for i in range(N):
    val = r0_Kf_values[i]
    result = all_results[i]
    print(f"r0.Kf={val}")
    # result.plotter.plot_concentrations() #doesn't have option to save plots
    t = result.times
    fig, ax = plt.subplots()
    ax.plot(t, result.concentrations.T)
    ax.set(xlabel="time (s)", ylabel="concentration", title="Concentration over time")

    y_labels = [c.label for c in result.channels if c.mean_values is not None]
    ax.legend(y_labels)
    ax.grid()
    plt.savefig(f"/scratch/xpz5km/cli_VCell/vcell_out/pyvcell_tutorial/{i}.pdf")
    plt.close()
