import gillespy2
import numpy as np
import matplotlib.pyplot as plt

n_trajectories = 1000
n_time_points = 1000*100
max_time = 1000

def create_soma(parameter_values=None):
    # First call the gillespy2.Model initializer.
    model = gillespy2.Model(name='Soma')

    # Define parameters for the rates of creation and dissociation.
    # n_gene_copies = gillespy2.Parameter(name='n_gene_copies', expression=1)
    gene_activation_rate_const = gillespy2.Parameter(name='gene_activation_rate_constant', expression=1)
    gene_deactivation_rate_const = gillespy2.Parameter(name='gene_deactivation_rate_constant', expression=0)
    transcription_rate_const = gillespy2.Parameter(name='transcription_rate_const', expression=(3*200/10000)*.001*3600 *1e3)
    mRNA_decay_rate_const = gillespy2.Parameter(name='mRNA_decay_rate_const', expression=1.2e-5*3600 * 1e1)

    translation_rate_const = gillespy2.Parameter(name='translation_rate_const', expression=0.021*3600 *0)
    protein_decay_rate_const = gillespy2.Parameter(name='protein_decay_rate_const', expression=1.21e-6*3600)

    model.add_parameter([gene_activation_rate_const, gene_deactivation_rate_const, transcription_rate_const, mRNA_decay_rate_const, translation_rate_const, protein_decay_rate_const])

    # Define variables for the molecular species
    inactive_genes = gillespy2.Species(name='inactive_genes', initial_value=1)
    active_genes = gillespy2.Species(name='active_genes', initial_value=0)
    mRNAs = gillespy2.Species(name='mRNAs',   initial_value=0)
    proteins = gillespy2.Species(name='proteins',   initial_value=0)

    model.add_species([inactive_genes, active_genes, mRNAs, proteins])

    # The list of reactants and products for a Reaction object are each a
    # Python dictionary in which the dictionary keys are Species objects
    # and the values are stoichiometries of the species in the reaction.
    gene_activation = gillespy2.Reaction(name="gene_activation", rate=gene_activation_rate_const, reactants={inactive_genes:1}, products={active_genes:1})
    gene_deactivation = gillespy2.Reaction(name="gene_deactivation", rate=gene_deactivation_rate_const, reactants={active_genes:1}, products={inactive_genes:1})

    mRNA_creation = gillespy2.Reaction(name="mRNA_creation", rate=transcription_rate_const, reactants={active_genes:1}, products={mRNAs:1, active_genes:1})
    mRNA_decay = gillespy2.Reaction(name="mRNA_decay", rate=mRNA_decay_rate_const, reactants={mRNAs:1}, products={})

    protein_creation = gillespy2.Reaction(name="protein_creation", rate=translation_rate_const, reactants={mRNAs:1}, products={proteins:1, mRNAs:1})
    protein_decay = gillespy2.Reaction(name="protein_decay", rate=protein_decay_rate_const, reactants={proteins:1}, products={})

    
    model.add_reaction([gene_activation, gene_deactivation, mRNA_creation, mRNA_decay, protein_creation, protein_decay])

    # Set the timespan for the simulation.
    tspan = gillespy2.TimeSpan.linspace(t=max_time, num_points=n_time_points)
    model.timespan(tspan)
    return model

model = create_soma()
results = model.run(number_of_trajectories=n_trajectories)

print('Done with Gillespie simulations')

gene_average = np.zeros(n_time_points)
mRNA_average = np.zeros(n_time_points)
prot_average = np.zeros(n_time_points)
for i in range(0, n_trajectories):
    gene_average += results[i]['active_genes']/n_trajectories
    mRNA_average += results[i]['mRNAs']/n_trajectories
    prot_average += results[i]['proteins']/n_trajectories

gene_var = np.zeros(n_time_points)
mRNA_var = np.zeros(n_time_points)
prot_var = np.zeros(n_time_points)
for i in range(0, n_trajectories):
    gene_var += (results[i]['active_genes'] - gene_average)**2/n_trajectories
    mRNA_var += (results[i]['mRNAs'] - mRNA_average)**2/n_trajectories
    prot_var += (results[i]['proteins'] - prot_average)**2/n_trajectories

fig_avrg, axs_avrg = plt.subplots(nrows=3, ncols=1, figsize=(10*1.6*1.1, 3.2*1.9*1.7))

axs_avrg[0].plot(results[0]['time'], gene_average)
axs_avrg[0].grid()
axs_avrg[1].plot(results[0]['time'], mRNA_average)
axs_avrg[1].grid()
axs_avrg[2].plot(results[0]['time'], prot_average)
axs_avrg[2].grid()

fig_var, axs_var = plt.subplots(nrows=3, ncols=1, figsize=(10*1.6*1.1, 3.2*1.9*1.7))

axs_var[0].plot(results[0]['time'], np.sqrt(gene_var))
axs_var[0].grid()
axs_var[1].plot(results[0]['time'], np.sqrt(mRNA_var))
axs_var[1].grid()
axs_var[2].plot(results[0]['time'], np.sqrt(prot_var))
axs_var[2].grid()

plt.show()
