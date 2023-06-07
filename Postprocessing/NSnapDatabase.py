###
# IMPORTS
###
import numpy as np
from os import path
import csv
import matplotlib.pyplot as plt
import NStaircase as sc
from io import TextIOWrapper
import zipfile
import matplotlib.animation as anim
import matplotlib as mpl
import datetime
import os
mpl.rcParams['animation.ffmpeg_path'] = r'C:\\Program Files (x86)\\FFMpeg\\bin\\ffmpeg.exe'

###
# MAIN CLASS
###

# Experiment List (ID=FileSpecifier explaining what differs from usual properties, sim_interval, stop_time):
class SnapDatabase:

    ###
    # DATA ACCESS
    ###

    # input: folder_path = 'Fig2/FigIIv1/', filespecifier = 'M1E3'
    #        zip = None (usual folder), zip file name
    def __init__(self, folder_path, file_specifier, zip):
        self.folder_path = folder_path
        self.file_specifier = file_specifier
        self.sim_interval = 1
        self.stop_time = 0
        self.l_prop = sc.LatProp()
        self.p_prop = []
        if zip is None:
            self.zip = False
        else:
            self.zip = True
            self.zip_file = zipfile.ZipFile(zip)
        self.read_preamble()

    # read preamble file
    # output: no output, changes the l_prop and p_prop of the whole class
    def read_preamble(self):
        snapshot_file = self.folder_path + "/OutSnapshot_"+self.file_specifier+".csv"
        if (not self.zip and path.exists(snapshot_file)) or (self.zip and (snapshot_file in self.zip_file.namelist())):
            if not self.zip:
                file = open(snapshot_file, 'r')
            else:
                file = TextIOWrapper(self.zip_file.open(snapshot_file, 'r'), 'utf-8')
            reader = list(csv.reader(file, skipinitialspace=True))
            if len(reader) != 0:
                self.sim_interval = float(reader[1][1])
                self.stop_time = float(reader[2][1])

                self.l_prop.PopulationNumber = int(reader[5][1])
                self.l_prop.GridBound = int(reader[6][1])
                self.l_prop.GenBound = int(reader[7][1])
                self.l_prop.CarryingCapacity = float(reader[8][1])
                self.l_prop.LD = int(reader[9][1])
                new_version = False
                if len(reader) > 36 and len(reader[35]) > 1:
                    if reader[35][0] == "CompetBelowStairs":   # new version
                        prop_num = (39 - 12 + 1) + 2  # number of property rows + extra two rows for indent and Population, ...
                        new_version = True
                else:
                    prop_num = (34 - 12 + 1) + 2  # number of property rows + extra two rows for indent and Population, ...
                for pop in range(self.l_prop.PopulationNumber):
                    self.p_prop.append(sc.PopProp())
                    self.p_prop[pop].InitCellsNum = int(reader[12+prop_num*pop][1])
                    self.p_prop[pop].InitCellsGen = int(reader[13+prop_num*pop][1])
                    self.p_prop[pop].InitCellsPos = int(reader[14+prop_num*pop][1])
                    self.p_prop[pop].BirthRate = float(reader[15+prop_num*pop][1])
                    self.p_prop[pop].StressBirthRate = float(reader[16+prop_num*pop][1])
                    self.p_prop[pop].ResistCost = float(reader[17+prop_num*pop][1])
                    self.p_prop[pop].DeathRate = float(reader[18+prop_num*pop][1])
                    self.p_prop[pop].StressDeathRate = float(reader[19+prop_num*pop][1])
                    self.p_prop[pop].MuteUp = float(reader[20+prop_num*pop][1])
                    self.p_prop[pop].MuteDown = float(reader[21+prop_num*pop][1])
                    self.p_prop[pop].StressMuteUp = float(reader[22+prop_num*pop][1])
                    self.p_prop[pop].StressMuteDown = float(reader[23+prop_num*pop][1])
                    self.p_prop[pop].Move = float(reader[24+prop_num*pop][1])
                    self.p_prop[pop].Chemotax = float(reader[25+prop_num*pop][1])
                    self.p_prop[pop].StressDepChemotax = bool(int(reader[26+prop_num*pop][1]))
                    self.p_prop[pop].Chemokin = float(reader[27+prop_num*pop][1])
                    self.p_prop[pop].SwitchUp = float(reader[28+prop_num*pop][1])
                    self.p_prop[pop].SwitchDown = float(reader[29+prop_num*pop][1])
                    self.p_prop[pop].ConsiderSwarm = bool(int(reader[30+prop_num*pop][1]))
                    self.p_prop[pop].SwarmingDensity = float(reader[31+prop_num*pop][1])
                    self.p_prop[pop].SwarmingMove = float(reader[32+prop_num*pop][1])
                    self.p_prop[pop].ConsiderHGT = bool(int(reader[33+prop_num*pop][1]))
                    self.p_prop[pop].HGTRate = float(reader[34+prop_num*pop][1])
                    if new_version:
                        self.p_prop[pop].CompetBelowStairs = float(reader[35+prop_num*pop][1])
                        self.p_prop[pop].ConsiderDensSwitch = bool(int(reader[36+prop_num*pop][1]))
                        self.p_prop[pop].DensSwitchBias = float(reader[37+prop_num*pop][1])
                        self.p_prop[pop].DensSwitchTot = float(reader[38+prop_num*pop][1])
                        self.p_prop[pop].DensSwitchDens = float(reader[39+prop_num*pop][1])

            else:
                print("Error! Cannot read the preamble, the snapshot_file is empty: " + snapshot_file)
                return None
        else:
            print("Error! Cannot read the preamble, snapshot_file does not exist: " + snapshot_file)
            return None

    # read snapshot file
    # output: [Time,GenSpaceTot[pop][gen][pos]]
    # input: sim_id(=0,...), task_id, population evolving mutation, LD, time_type (='T','D','PA') -> 'OutAdap_taskID_ID_pop_LDT/D.csv'
    def read_gen_space_tot(self, sim_index):
        snapshot_file = self.folder_path + "/OutSnapshot_"+self.file_specifier+"_"+str(sim_index)+".csv"
        GenSpaceTot = np.zeros((self.l_prop.PopulationNumber, self.l_prop.GenBound, self.l_prop.GridBound), dtype=float)
        if (not self.zip and path.exists(snapshot_file)) or (self.zip and (snapshot_file in self.zip_file.namelist())):
            if not self.zip:
                file = open(snapshot_file, 'r')
            else:
                file = TextIOWrapper(self.zip_file.open(snapshot_file, 'r'), 'utf-8')
            reader = list(csv.reader(file, skipinitialspace=True))
            if len(reader) != 0:
                Time = float(reader[1][1])
                for pop in range(self.l_prop.PopulationNumber):
                    for gen in range(self.l_prop.GenBound):
                        for pos in range(self.l_prop.GridBound):
                            GenSpaceTot[pop][gen][pos] = float(reader[pop*(self.l_prop.GenBound+1)+gen+3][pos])
                return [Time, GenSpaceTot]
            else:
                print("Warning! Tried to read gen_space_tot but the snapshot_file is empty: " + snapshot_file)
                return None
        else:
            print("Warning! Tried to read gen_space_tot but the snapshot_file does not exist: " + snapshot_file)
            return None

    # calculate gen_tot from gen_space_tot
    def gen_tot(self, gen_space_tot):
        gen_tot = np.zeros((self.l_prop.PopulationNumber, self.l_prop.GenBound))
        for pop in range(self.l_prop.PopulationNumber):
            for gen in range(self.l_prop.GenBound):
                for pos in range(self.l_prop.GridBound):
                    gen_tot[pop][gen] += gen_space_tot[pop][gen][pos]

    # calculate space_tot from gen_space_tot
    def space_tot(self, gen_space_tot):
        space_tot = np.zeros((self.l_prop.PopulationNumber, self.l_prop.GridBound))
        for pop in range(self.l_prop.PopulationNumber):
            for gen in range(self.l_prop.GenBound):
                for pos in range(self.l_prop.GridBound):
                    space_tot[pop][pos] += gen_space_tot[pop][gen][pos]
        return space_tot

    # calculate abs_gen_tot from gen_space_tot
    def abs_gen_tot(self, gen_space_tot):
        gen_tot = np.zeros(self.l_prop.GenBound)
        for pop in range(self.l_prop.PopulationNumber):
            for gen in range(self.l_prop.GenBound):
                for pos in range(self.l_prop.GridBound):
                    gen_tot[gen] += gen_space_tot[pop][gen][pos]
        return gen_tot

    # calculate abs_space_tot from gen_space_tot
    def abs_space_tot(self, gen_space_tot):
        space_tot = np.zeros(self.l_prop.GridBound)
        for pop in range(self.l_prop.PopulationNumber):
            for gen in range(self.l_prop.GenBound):
                for pos in range(self.l_prop.GridBound):
                    space_tot[pos] += gen_space_tot[pop][gen][pos]
        return space_tot

    # calculate abs_tot from gen_space_tot
    def abs_tot(self, gen_space_tot):
        abs_tot = 0
        for pop in range(self.l_prop.PopulationNumber):
            for gen in range(self.l_prop.GenBound):
                for pos in range(self.l_prop.GridBound):
                    abs_tot += gen_space_tot[pop][gen][pos]
        return abs_tot

    # sum over populations, [pos][gen]
    def space_gen_tot(self, gen_space_tot):
        space_gen_tot = np.zeros((self.l_prop.GridBound, self.l_prop.GenBound))
        for pop in range(self.l_prop.PopulationNumber):
            for gen in range(self.l_prop.GenBound):
                for pos in range(self.l_prop.GridBound):
                    space_gen_tot[pos][gen] += gen_space_tot[pop][gen][pos]
        return space_gen_tot

    # calculate horizontal gene flux in population pop
    def horizontal_gene_flux(self, pop, gen_space_tot):
        gene_flux = 0
        for gen in range(self.l_prop.GenBound-1):
            for pos in range(self.l_prop.GridBound):
                gen_donor = gen+1
                while gen_donor<self.l_prop.GenBound:
                    if self.p_prop[pop].ConsiderHGT:
                        gene_flux += self.p_prop[pop].HGTRate*gen_space_tot[pop][gen][pos]*gen_space_tot[pop][gen_donor][pos]
                    gen_donor += 1
        return gene_flux

    # calculate cumulative gene flux in population pop
    def cum_horizontal_gene_flux(self, times, fluxes, sim_index):
        cum_gene_flux = np.trapz(fluxes[:(sim_index+1)], times[:(sim_index+1)])
        return cum_gene_flux

    # finds the resistance number LD (=R) for a given GenSpaceTot
    def LD(self, gen_space_tot):
        pos = self.l_prop.GridBound - 1
        space_gen_tot = self.space_gen_tot(gen_space_tot)
        while pos >= 0:
            if pos == np.argmax(space_gen_tot[pos]):
                break
            pos -= 1
        if pos == -1:
            LD = np.argmax(self.abs_gen_tot(gen_space_tot)) + 1
        else:
            LD = pos+1
        return LD

    # finds the swarming number S for a given swarm_pop and GenSpaceTot
    def S(self, gen_space_tot, swarm_pop):
        if self.p_prop[swarm_pop].ConsiderSwarm:
            S = 0
            space_tot = self.space_tot(gen_space_tot)[swarm_pop]
            while S < self.l_prop.GridBound and space_tot[S] >= self.p_prop[swarm_pop].SwarmingDensity*0.99:
                S += 1
            return S
        else:
            return None

    # finds if there is a founder of the next state, if currently in resistance state R
    def next_founder(self, gen_space_tot, R):
        founder = False
        if R<min(self.l_prop.GenBound-1,self.l_prop.GridBound-1):
            for pop in range(self.l_prop.PopulationNumber):
                if gen_space_tot[pop][R+1][R+1] > 0:    # i.e., there is a founder and this founder managed to reproduce
                    founder = True
        return founder

    ####
    # PLOTTING FUNCTIONS
    ####

    # plots the staircase plot into given axis
    # input: sim_index, ax(axes to plot to), labels(=true,false), plot_pop (pop. to be plotted)
    def staircase_plot(self, sim_index, plot_pop, ax, labels):
        # get important data
        [Time, GenSpaceTot] = self.read_gen_space_tot(sim_index)
        # make the plot
        x_axis = np.asarray([0.5+pos for pos in range(self.l_prop.GridBound+1)])
        g_axis = np.asarray([0.5+gen for gen in range(self.l_prop.GenBound+1)])
        ax.pcolormesh(x_axis, g_axis, GenSpaceTot[plot_pop], vmin=0, vmax=self.l_prop.CarryingCapacity, cmap='Greys')
        for j in range(self.l_prop.GridBound-1):
            ax.hlines(j + 1.5, j + 1.5, j + 2.5, colors='black')
            ax.vlines(j + 1.5, j + 0.5, j + 1.5, colors='black')
        ax.set_xticks([1+pos for pos in range(self.l_prop.GridBound)])
        ax.set_yticks([1+pos for pos in range(self.l_prop.GridBound)])
        # make lables
        if not labels:
            ax.set_xticks([])
            ax.set_yticks([])
        else:
            ax.set_xlabel('space x')
            ax.set_ylabel('genotype g')
            ax.set_title('t = '+str(Time)+' h')

    # plots the stable wild-type plot into given axis
    # input: sim_index, ax(axes to plot to), labels(=true,false)
    def wild_type_plot(self, sim_index, plot_r_net, r_labels, pop_labels, x_label, special_x_tick, ax, ax_twin, analytics):
        # get important data
        if self.read_gen_space_tot(sim_index) is None:
            [Time,GenSpaceTot] = [np.inf, np.zeros((self.l_prop.PopulationNumber, self.l_prop.GridBound, self.l_prop.GridBound))]
        else:
            [Time, GenSpaceTot] = self.read_gen_space_tot(sim_index)
        LD = self.LD(GenSpaceTot)

        # make stable wild-type plot
        # colors of increasing shade [blue, yellow, green, orange], can later add: purple, grey from https://htmlcolorcodes.com/
        colors = [['#AED6F1', '#85C1E9', '#5DADE2', '#3498DB', '#2E86C1', '#2874A6', '#21618C', '#1B4F72'],
                  ['#F9E79F', '#F7DC6F', '#F4D03F', '#F1C40F', '#D4AC0D', '#B7950B', '#9A7D0A', '#7D6608'],
                  ['#ABEBC6', '#82E0AA', '#58D68D', '#2ECC71', '#28B463', '#239B56', '#1D8348', '#186A3B'],
                  ['#F5CBA7', '#F0B27A', '#EB984E', '#E67E22', '#CA6F1E', '#AF601A', '#935116', '#784212']]
        ax.tick_params(axis='y', labelcolor='black')
        if special_x_tick:
            ax.set_xticks([1, LD + 1, self.l_prop.GridBound])
            ax.set_xticklabels(["$1$", "$R+1$", "$L$"])
        else:
            ax.set_xticks([1+k for k in range(self.l_prop.GridBound)])
        ax.set_ylim([0, self.l_prop.CarryingCapacity])
        if self.p_prop[0].Chemotax != 0.5: # for chemotaxis
            ax.set_ylim([0, 7 * self.l_prop.CarryingCapacity / 2])
        ax.set_xlim([1, self.l_prop.GridBound])
        ax.set_yticks([])

        # plot profile if a single swarming population is present
        if self.l_prop.PopulationNumber == 1 and self.p_prop[0].ConsiderSwarm:
            S = self.S(GenSpaceTot, 0)
            if self.p_prop[0].Move > self.p_prop[0].SwarmingMove:
                colors[0], colors[1] = colors[1], colors[0]
            for g in range(self.l_prop.GenBound):
                g = self.l_prop.GenBound - 1 - g
                N_sim = np.zeros(2*self.l_prop.GridBound-1)
                for pos in range(self.l_prop.GridBound):
                    for g_i in range(g + 1):
                        N_sim[2*pos] += GenSpaceTot[0][g_i][pos]
                    if pos > 0:
                        N_sim[2*pos-1] = (N_sim[2*pos]+N_sim[2*(pos-1)])/2
                if g == self.l_prop.GenBound - 1:
                    if S == 0:
                        ax.fill_between([k/2 + 1 for k in range(2*self.l_prop.GridBound-1)], N_sim, color=colors[0][g], label='')
                    elif S == self.l_prop.GridBound:
                        ax.fill_between([k/2 + 1 for k in range(2*self.l_prop.GridBound-1)], N_sim, color=colors[1][g], label='')
                    else:
                        ax.fill_between([k/2 + 1 for k in range(2*S)], N_sim[:2*S], color=colors[1][g], label='')
                        ax.fill_between([k/2 + S + 1/2 for k in range(2*self.l_prop.GridBound-2*S)], N_sim[2*S-1:], color=colors[0][g], label='simulation')
                    if not analytics:
                        ax.plot([k/2 + 1 for k in range(2*self.l_prop.GridBound-1)], N_sim, color='black')
                else:
                    if S == 0:
                        ax.fill_between([k/2 + 1 for k in range(2*self.l_prop.GridBound-1)], N_sim, color=colors[0][g], label='')
                    elif S == self.l_prop.GridBound:
                        ax.fill_between([k/2 + 1 for k in range(2*self.l_prop.GridBound-1)], N_sim, color=colors[1][g], label='')
                    else:
                        ax.fill_between([k/2 + 1 for k in range(2*S)], N_sim[:2*S], color=colors[1][g], label='')
                        ax.fill_between([k/2 + S + 1/2 for k in range(2*self.l_prop.GridBound-2*S)], N_sim[2*S-1:], color=colors[0][g], label='')
        # plot a profile otherwise
        else:
            for plot_pop in range(self.l_prop.PopulationNumber):
                color = plot_pop
                if self.p_prop[0].ConsiderDensSwitch:
                    if self.p_prop[0].Move < self.p_prop[1].Move:
                        color = self.l_prop.PopulationNumber - 1 - plot_pop
                plot_pop = self.l_prop.PopulationNumber - 1 - plot_pop
                for g in range(self.l_prop.GenBound):
                    g = self.l_prop.GenBound - 1 - g
                    N_sim = np.zeros(self.l_prop.GridBound)
                    for x in range(self.l_prop.GridBound):
                        for g_i in range(g + 1):
                            for plot_pop_i in range(plot_pop + 1):
                                N_sim[x] += GenSpaceTot[plot_pop_i][g_i][x]
                    if g == self.l_prop.GenBound - 1 and plot_pop == self.l_prop.PopulationNumber - 1:
                        ax.fill_between([k + 1 for k in range(self.l_prop.GridBound)], N_sim, color=colors[color][g], label='simulation')
                        if not analytics:
                            ax.plot([k + 1 for k in range(self.l_prop.GridBound)], N_sim, color='black')
                    elif g == self.l_prop.GenBound - 1:
                        ax.fill_between([k + 1 for k in range(self.l_prop.GridBound)], N_sim, color=colors[color][g])
                        if not analytics:
                            ax.plot([k + 1 for k in range(self.l_prop.GridBound)],  N_sim, linestyle='--', color='black')
                    else:
                        ax.fill_between([k + 1 for k in range(self.l_prop.GridBound)], N_sim, color=colors[color][g], label='')
        # if net birth rate to be plotted
        if plot_r_net:
            plot_pop = 0  # assumes identical birth and death rates for all populations
            # prepare data for r_net plot
            r_net = np.zeros(self.l_prop.GridBound)
            N_tot = np.zeros(self.l_prop.GridBound)
            for pos in range(self.l_prop.GridBound):
                for pop in range(self.l_prop.PopulationNumber):
                    for gen in range(self.l_prop.GenBound):
                        N_tot[pos] += GenSpaceTot[pop][gen][pos]
                # birth rate
                if N_tot[pos] < self.l_prop.CarryingCapacity:
                    if pos < LD+1:
                        r_net[pos] += self.p_prop[plot_pop].BirthRate*(1-N_tot[pos]/self.l_prop.CarryingCapacity)*((1-self.p_prop[plot_pop].ResistCost)**(LD-pos))
                    else:
                        r_net[pos] += self.p_prop[plot_pop].StressBirthRate*(1-N_tot[pos]/self.l_prop.CarryingCapacity)
                # death rate:
                if pos < LD + 1:
                    r_net[pos] -= self.p_prop[plot_pop].DeathRate
                else:
                    r_net[pos] -= (self.p_prop[plot_pop].DeathRate+self.p_prop[plot_pop].StressDeathRate)
            # make the plot
            ax_twin.bar([k + 1 for k in range(self.l_prop.GridBound)], r_net, color='tab:red')
            ax_twin.tick_params(axis='y', labelcolor='tab:red')
            ax_twin.set_ylim([-1.5 * self.p_prop[plot_pop].DeathRate, self.p_prop[plot_pop].BirthRate-self.p_prop[plot_pop].DeathRate/2])
            ax_twin.set_yticks([])
            # make r_net labels
            if r_labels:
                ax_twin.set_ylabel('mutant net birth rate', color='tab:red')
                ax_twin.tick_params(axis='y', labelcolor='tab:red')
                yticks = [-self.p_prop[plot_pop].DeathRate, 0, self.p_prop[plot_pop].DeathRate]
                ylabels = ['$-\delta$', '$0$', '$\delta$']
                for i in range(int((self.p_prop[plot_pop].BirthRate-self.p_prop[plot_pop].DeathRate/2)/self.p_prop[plot_pop].DeathRate)):
                    if i > 0:
                        yticks.append((i + 1) * self.p_prop[plot_pop].DeathRate)
                        ylabels.append(str(i + 1) + "$\\delta$")
                ax_twin.set_yticks(yticks)
                ax_twin.set_yticklabels(ylabels)
        # make labels for the population profile graph
        if pop_labels:
            K = self.l_prop.CarryingCapacity
            ax.set_ylabel('$N_x$=$\\#$ cells', color='black')
            #ax.set_yticks([0, K / 2, K, 3 * K / 2, 2*K, 5*K/2, 3*K, 7*K/2])  # for chemotaxis
            #ax.set_yticklabels(["0", "K/2", "K", "3K/2", "2K", "5K/2", "3K", "7K/2"])  # for chemotaxis
            ax.set_yticks([0, K / 4, K / 2, 3 * K / 4, K])
            ax.set_yticklabels(["0", "K/4", "K/2", "3K/4", "K"])
            if self.p_prop[0].Chemotax != 0.5:  # for chemotaxis
                ax.set_yticks([0, K / 2, K, 3 * K / 2, 2 * K, 5 * K / 2, 3 * K, 7 * K / 2])
                ax.set_yticklabels(["0", "K/2", "K", "3K/2", "2K", "5K/2", "3K", "7K/2"])
        if x_label:
            ax.set_xlabel('space x')
        # make analytics
        if analytics:
            self.l_prop.LD = LD
            gen_bound = self.l_prop.GenBound
            self.l_prop.GenBound = 1
            staircase = sc.Staircase(self.l_prop, self.p_prop)
            self.l_prop.GenBound = gen_bound
            N_plot = np.zeros(self.l_prop.GridBound)
            for plot_pop in range(self.l_prop.PopulationNumber):
                for pos in range(self.l_prop.GridBound):
                    N_plot[pos] += staircase.N[plot_pop * self.l_prop.GridBound + pos]
                if plot_pop < self.l_prop.PopulationNumber - 1:
                    ax.plot([i + 1 for i in range(self.l_prop.GridBound)], N_plot, linestyle='--', color='black',
                            label='analytics')
                else:
                    ax.plot([i + 1 for i in range(self.l_prop.GridBound)], N_plot, color='black',
                            label='analytics')

    # produces an animation of staircase and wild_type plots
    # input: staircase (t/f), wildtype (t/f), max_time, fps speed (60 ususal), animation_name
    def animated_plot(self, staircase, wildtype, max_time, sim_delay, fps, titles, animation_name):
        # prepare figure properties
        staircase_num = int(staircase)*self.l_prop.PopulationNumber
        wildtype_num = int(wildtype)
        staircase_width = 5.5
        wildtype_width = 5.5
        width = staircase_width*staircase_num+wildtype_width*wildtype_num
        height = 5
        width_ratio = [staircase_width for i in range(staircase_num)]+[wildtype_width for i in range(wildtype_num)]
        title_colors = ['#2E86C1', '#D4AC0D', '#28B463', '#CA6F1E']
        # create figure
        fig, axes = plt.subplots(nrows=1, ncols=staircase_num + wildtype_num, figsize=(width, height),
                                gridspec_kw={'width_ratios': width_ratio})
        axs = []
        if staircase_num + wildtype_num == 1:
            axs.append(axes)
        else:
            for plot in range(staircase_num + wildtype_num):
                axs.append(axes[plot])
        if wildtype:
            ax_twin = axs[staircase_num].twinx()
            axs.append(ax_twin)
        # prepare folder
        date_time = datetime.datetime.now()
        folder_specifier = "_"+str(date_time.date().day)+"_"+str(date_time.date().month)+"_"+str(date_time.date().year)+"_"\
                           + str(date_time.time().hour) +"_"+str(date_time.time().minute)
        folder_name = "RawMovies\\"+animation_name+folder_specifier
        os.mkdir(folder_name)
        # draw the figure for each time and take a snapshot
        max_time = min(max_time, self.stop_time)
        for plot_index in range(int(max_time/sim_delay)):
            sim_index = int(sim_delay/self.sim_interval)*plot_index
            # plot staircase
            if staircase:
                for plot_pop in range(staircase_num):
                    self.staircase_plot(sim_index, self.l_prop.PopulationNumber - 1 - plot_pop, axs[plot_pop], True)
                    if staircase_num > 1:
                        axs[plot_pop].set_title(titles[plot_pop], color=title_colors[plot_pop])
            # plot wildtype
            if wildtype:
                self.wild_type_plot(sim_index, True, True, True, True, False, axs[staircase_num], axs[staircase_num+1], False)
            # read gen_space_tot
            if self.read_gen_space_tot(sim_index) is None:
                [Time, GenSpaceTot] = [np.inf, np.zeros((self.l_prop.PopulationNumber, self.l_prop.GridBound, self.l_prop.GridBound))]
            else:
                [Time, GenSpaceTot] = self.read_gen_space_tot(sim_index)
            # plot title
            if self.l_prop.PopulationNumber == 1 and self.p_prop[0].ConsiderSwarm:
                R = self.LD(GenSpaceTot)
                S = self.S(GenSpaceTot, 0)
                fig.suptitle('Time: '+"{:6.1f}".format(Time)+"h, R: "+str(R)+", S: "+str(S))
            else:
                R = self.LD(GenSpaceTot)
                fig.suptitle('Time: '+"{:6.1f}".format(Time)+"h, R: "+str(R))
            # take a snapshot and clear the canvas
            fig.savefig(folder_name+"\\"+animation_name+"-"+str(plot_index+1)+".png", dpi=300)  # save the figure to file
            for ax in axs:
                for artist in ax.lines + ax.collections + ax.patches:
                    artist.remove()
            print("Created frame: "+str(plot_index))

        # run ffmpeg in the command prompt
        # ffmpeg command such as H:\Můj disk\Summer Research 2021\Staircase model\Simulations\RawMovies\LowMotility_LowHGT_24_3_2022_10_41>ffmpeg -framerate 25 -i LowMotility_LowHGT-%d.png output.mp4
        os.system("h:; "
                  "cd H:\\Můj disk\\Summer Research 2021\\Staircase model\\Simulations\\"+folder_name+"; "
                  "ffmpeg -framerate "+str(fps)+" -i "+animation_name+"-%d.png "+animation_name+".mov")

    # produces a snapshot in the animated plot
    # input: staircase (t/f), wildtype (t/f), simulated time
    def snap_draw_plot(self, staircase, wildtype, sim_time, titles):
        # prepare figure properties
        staircase_num = int(staircase)*self.l_prop.PopulationNumber
        wildtype_num = int(wildtype)
        staircase_width = 3
        wildtype_width = 3
        width = staircase_width*staircase_num+wildtype_width*wildtype_num
        height = 3
        width_ratio = [staircase_width for i in range(staircase_num)]+[wildtype_width for i in range(wildtype_num)]
        title_colors = ['#2E86C1', '#D4AC0D', '#28B463', '#CA6F1E']
        # create figure
        fig, axes = plt.subplots(nrows=1, ncols=staircase_num + wildtype_num, figsize=(width, height),
                                gridspec_kw={'width_ratios': width_ratio}, constrained_layout=True)
        axs = []
        if staircase_num + wildtype_num == 1:
            axs.append(axes)
        else:
            for plot in range(staircase_num + wildtype_num):
                axs.append(axes[plot])
        if wildtype:
            ax_twin = axs[staircase_num].twinx()
            axs.append(ax_twin)
        # draw the figure for each time and take a snapshot
        sim_time = min(sim_time, self.stop_time)
        sim_index = int(sim_time/self.sim_interval)
        # plot staircase
        if staircase:
            for plot_pop in range(staircase_num):
                self.staircase_plot(sim_index, plot_pop, axs[plot_pop], True)
                if staircase_num > 1:
                    axs[plot_pop].set_title(titles[plot_pop], color=title_colors[plot_pop])
        # plot wildtype
        if wildtype:
            self.wild_type_plot(sim_index, True, True, True, True, False, axs[staircase_num], axs[staircase_num+1], False)
        # read gen_space_tot
        if self.read_gen_space_tot(sim_index) is None:
            [Time, GenSpaceTot] = [np.inf, np.zeros((self.l_prop.PopulationNumber, self.l_prop.GridBound, self.l_prop.GridBound))]
        else:
            [Time, GenSpaceTot] = self.read_gen_space_tot(sim_index)
        # plot title
        if self.l_prop.PopulationNumber == 1 and self.p_prop[0].ConsiderSwarm:
            R = self.LD(GenSpaceTot)
            S = self.S(GenSpaceTot, 0)
            fig.suptitle('Time: '+"{:6.1f}".format(Time)+", R: "+str(R)+", S: "+str(S))
        else:
            R = self.LD(GenSpaceTot)
            fig.suptitle('Time: '+"{:6.1f}".format(Time)+", R: "+str(R))
        plt.show()

    # plots the horizontal gene flux of pop as a function of time
    # input: (HGT population) pop, (plotting axes) ax, (maximal plotted time) = max_time, (name of line in legend) legend_name, lables (t/f)
    def gene_flux_plot(self, pop, max_time, ax, legend_name, labels):
        if max_time is None:
            sim_num = int(self.stop_time/self.sim_interval)
        else:
            sim_num = int(max_time / self.sim_interval)
        times = np.zeros(sim_num)
        fluxes = np.zeros(sim_num)
        for sim_index in range(sim_num):
            [time, gen_space_tot] = self.read_gen_space_tot(sim_index)
            times[sim_index] = time
            fluxes[sim_index] = self.horizontal_gene_flux(pop, gen_space_tot)
        ax.plot(times, fluxes, label=legend_name)
        if labels:
            ax.set_xlabel("time t (h)", fontsize=12)
            ax.set_ylabel("gene flux (1/h)", fontsize=12)

    # plots the cumulative horizontal gene flux of pop as a function of time
    # input: (HGT population) pop, (plotting axes) ax, (maximal plotted time) = max_time, (name of line in legend) legend_name, lables (t/f)
    def cum_gene_flux_plot(self, pop, max_time, ax, legend_name, labels):
        if max_time is None:
            sim_num = int(self.stop_time/self.sim_interval)
        else:
            sim_num = int(max_time / self.sim_interval)
        times = np.zeros(sim_num)
        fluxes = np.zeros(sim_num)
        cum_fluxes = np.zeros(sim_num)
        for sim_index in range(sim_num):
            [time, gen_space_tot] = self.read_gen_space_tot(sim_index)
            times[sim_index] = time
            fluxes[sim_index] = self.horizontal_gene_flux(pop, gen_space_tot)
        for sim_index in range(sim_num):
            cum_fluxes[sim_index] = self.cum_horizontal_gene_flux(times, fluxes, sim_index)
        ax.plot(times, cum_fluxes, label=legend_name)
        if labels:
            ax.set_xlabel("time t (h)", fontsize=12)
            ax.set_ylabel("cumulative gene flux (1/h)", fontsize=12)

    # plots the evolution of resistance (R) and swarming (S) as a function of time
    def resistance_swarming_plot(self, max_time, swarm_pop, sim_delay, ax, labels, show_tolerance):
        times_S = []
        times_R = []
        S_values = []
        R_values = []
        T_values = []
        times_T = []
        R = 0
        S = -1
        # go through snapshots and save the times when S or R changes
        max_time = min(max_time, self.stop_time)
        for plot_index in range(int(max_time / sim_delay)):
            sim_index = int(sim_delay / self.sim_interval) * plot_index
            [time, gen_space_tot] = self.read_gen_space_tot(sim_index)
            print("Current time: "+str(time))
            something_changed = False
            if self.S(gen_space_tot, swarm_pop)>S:
                times_S.append(time)
                S += 1
                S_values.append(S)
                something_changed = True
            if self.LD(gen_space_tot)>R:
                times_R.append(time)
                R += 1
                R_values.append(R)
                something_changed = True
            if something_changed:
                times_T.append(time)
                T_values.append(S-R)
        # make the plot
        ax.plot(times_S, S_values, marker='o', label="swarming S", color="tab:green")
        ax.plot(times_R, R_values, marker='o', label="resistance R", color="tab:red")
        if show_tolerance:
            ax.plot(times_T, T_values, marker='o', label="tolerance T=S-R")
        # set limits, etc...
        ax.set_ylim([0, self.l_prop.GridBound])
        ax.set_xlim([0, max_time])
        if labels:
            ax.legend()
            ax.set_xlabel("time t (h)")
            ax.set_ylabel("swarming and resistance levels")

###
# OTHER FUNCTIONS
###

# produces an animation of staircase and wild_type plots from mutliple simulations
# the code requires: SAME POPULATION NUMBER, SAME TITLES FOR RESPECTIVE POPULATIONS
# input: staircase (t/f), wildtype (t/f), max_time, fps speed (60 ususal), animation_name
def animated_plot(simulations, staircase, wildtype, max_time, sim_delay, titles, fig_titles, animation_name):
    # prepare figure properties
    sim_num = len(simulations)
    staircase_width = 5.5
    wildtype_width = 5.5
    staircase_num = []
    wildtype_num = int(wildtype)
    width = []
    for sim in range(sim_num):
        fig_titles[sim] = r"$\bf{" + fig_titles[sim].replace(" ", r"}$ $\bf{") + "}$"
        staircase_num_local = int(staircase) * simulations[sim].l_prop.PopulationNumber
        staircase_num.append(staircase_num_local)
        width.append(staircase_width * staircase_num_local + wildtype_width * wildtype_num)
    height = 5*sim_num
    title_colors = ['#2E86C1', '#D4AC0D', '#28B463', '#CA6F1E']
    # create figure
    fig = plt.figure(figsize=(max(width), height))
    subfigs = fig.subfigures(nrows=sim_num, ncols=1)
    if sim_num == 1:
        subfigs = [subfigs]
    axs = []
    for row, subfig in enumerate(subfigs):
        margin_width = (max(width)-width[row])/(staircase_num[row]+wildtype_width+1)
        if margin_width > 0:
            width_ratio = [margin_width]
            for i in range(staircase_num[row]):
                width_ratio.append(staircase_width*0.5)
                width_ratio.append(margin_width)
            for i in range(wildtype_num):
                width_ratio.append(wildtype_width*0.5)
                width_ratio.append(margin_width)
            axes = subfig.subplots(nrows=1, ncols=2*(staircase_num[row] + wildtype_num) + 1, gridspec_kw={'width_ratios': width_ratio})
            sub_axs = []
            for plot in range(staircase_num[row] + wildtype_num):
                sub_axs.append(axes[2*plot+1])
                axes[2*plot].grid(False)
                axes[2*plot].axis('off')
            if wildtype:
                ax_twin = sub_axs[staircase_num[row]].twinx()
                sub_axs.append(ax_twin)
            axs.append(sub_axs)
            axes[-1].grid(False)
            axes[-1].axis('off')
        else:
            width_ratio = [staircase_width for i in range(staircase_num[row])] + [wildtype_width for i in range(wildtype_num)]
            axes = subfig.subplots(nrows=1, ncols=staircase_num[row] + wildtype_num, gridspec_kw={'width_ratios': width_ratio})
            sub_axs = []
            if staircase_num[row] + wildtype_num == 1:
                sub_axs.append(axes)
            else:
                for plot in range(staircase_num[row] + wildtype_num):
                    sub_axs.append(axes[plot])
            if wildtype:
                ax_twin = sub_axs[staircase_num[row]].twinx()
                sub_axs.append(ax_twin)
            axs.append(sub_axs)
    # prepare folder
    date_time = datetime.datetime.now()
    folder_specifier = "_" + str(date_time.date().day) + "_" + str(date_time.date().month) + "_" + str(
        date_time.date().year) + "_" \
                       + str(date_time.time().hour) + "_" + str(date_time.time().minute)
    folder_name = "RawMovies\\" + animation_name + folder_specifier
    os.mkdir(folder_name)
    # draw the figure for each time and take a snapshot
    stop_times = [simulations[sim].stop_time for sim in range(sim_num)]
    max_time = min([max_time]+stop_times)
    #R = [0 for sim in range(sim_num)]
    for plot_index in range(int(max_time / sim_delay)):
        for sim in range(sim_num):
            sim_index = int(sim_delay / simulations[sim].sim_interval) * plot_index
            # plot staircase
            if staircase:
                for plot_pop in range(staircase_num[sim]):
                    if simulations[sim].p_prop[0].ConsiderDensSwitch:
                        simulations[sim].staircase_plot(sim_index, plot_pop, axs[sim][plot_pop], True)
                        if staircase_num[sim] > 1:
                            if simulations[sim].p_prop[0].Move > simulations[sim].p_prop[1].Move:
                                axs[sim][plot_pop].set_title(titles[plot_pop], color=title_colors[simulations[sim].l_prop.PopulationNumber - 1 - plot_pop])
                            else:
                                axs[sim][plot_pop].set_title(titles[plot_pop], color=title_colors[plot_pop])
                    else:
                        simulations[sim].staircase_plot(sim_index, simulations[sim].l_prop.PopulationNumber - 1 - plot_pop, axs[sim][plot_pop], True)
                        if staircase_num[sim] > 1:
                            axs[sim][plot_pop].set_title(titles[plot_pop], color=title_colors[plot_pop])
            # plot wildtype
            if wildtype:
                simulations[sim].wild_type_plot(sim_index, True, True, True, True, False, axs[sim][staircase_num[sim]], axs[sim][staircase_num[sim] + 1], False)
            # read gen_space_tot
            if simulations[sim].read_gen_space_tot(sim_index) is None:
                [Time, GenSpaceTot] = [np.inf, np.zeros((simulations[sim].l_prop.PopulationNumber, simulations[sim].l_prop.GridBound, simulations[sim].l_prop.GridBound))]
            else:
                [Time, GenSpaceTot] = simulations[sim].read_gen_space_tot(sim_index)
            # plot title
            R = simulations[sim].LD(GenSpaceTot)
            #if simulations[sim].next_founder(GenSpaceTot, R[sim]):
            #    R[sim] += 1
            subfigs[sim].suptitle(fig_titles[sim]+", Time: " + "{:6.1f}".format(Time) + "h, R: " + str(R))
        # take a snapshot and clear the canvas
        #fig.savefig(folder_name + "\\" + animation_name + "-" + str(plot_index + 1) + ".svg")  # save the figure to file
        fig.savefig(folder_name + "\\" + animation_name + "-" + str(plot_index + 1) + ".png", dpi=300)  # save the figure to file
        for sim in range(sim_num):
            for ax in axs[sim]:
                for artist in ax.lines + ax.collections + ax.patches:
                    artist.remove()
        print("Created frame: " + str(plot_index))

# produces a specific snapshot from the animation of staircase and wild_type plots from mutliple simulations
# the code requires: SAME POPULATION NUMBER, SAME TITLES FOR RESPECTIVE POPULATIONS
# input: staircase (t/f), wildtype (t/f), max_time, fps speed (60 ususal), animation_name
def snap_animated_plot(plot_index, simulations, staircase, wildtype, max_time, sim_delay, titles, fig_titles, animation_name):
    # prepare figure properties
    gen_space_tot = []
    sim_num = len(simulations)
    staircase_width = 5.5
    wildtype_width = 5.5
    staircase_num = []
    wildtype_num = int(wildtype)
    width = []
    for sim in range(sim_num):
        fig_titles[sim] = r"$\bf{" + fig_titles[sim].replace(" ", r"}$ $\bf{") + "}$"
        staircase_num_local = int(staircase) * simulations[sim].l_prop.PopulationNumber
        staircase_num.append(staircase_num_local)
        width.append(staircase_width * staircase_num_local + wildtype_width * wildtype_num)
    height = 5*sim_num
    title_colors = ['#2E86C1', '#D4AC0D', '#28B463', '#CA6F1E']
    # create figure
    fig = plt.figure(figsize=(max(width), height))
    subfigs = fig.subfigures(nrows=sim_num, ncols=1)
    if sim_num == 1:
        subfigs = [subfigs]
    axs = []
    for row, subfig in enumerate(subfigs):
        margin_width = (max(width)-width[row])/(staircase_num[row]+wildtype_width+1)
        if margin_width > 0:
            width_ratio = [margin_width]
            for i in range(staircase_num[row]):
                width_ratio.append(staircase_width*0.5)
                width_ratio.append(margin_width)
            for i in range(wildtype_num):
                width_ratio.append(wildtype_width*0.5)
                width_ratio.append(margin_width)
            axes = subfig.subplots(nrows=1, ncols=2*(staircase_num[row] + wildtype_num) + 1, gridspec_kw={'width_ratios': width_ratio})
            sub_axs = []
            for plot in range(staircase_num[row] + wildtype_num):
                sub_axs.append(axes[2*plot+1])
                axes[2*plot].grid(False)
                axes[2*plot].axis('off')
            if wildtype:
                ax_twin = sub_axs[staircase_num[row]].twinx()
                sub_axs.append(ax_twin)
            axs.append(sub_axs)
            axes[-1].grid(False)
            axes[-1].axis('off')
        else:
            width_ratio = [staircase_width for i in range(staircase_num[row])] + [wildtype_width for i in range(wildtype_num)]
            axes = subfig.subplots(nrows=1, ncols=staircase_num[row] + wildtype_num, gridspec_kw={'width_ratios': width_ratio})
            sub_axs = []
            if staircase_num[row] + wildtype_num == 1:
                sub_axs.append(axes)
            else:
                for plot in range(staircase_num[row] + wildtype_num):
                    sub_axs.append(axes[plot])
            if wildtype:
                ax_twin = sub_axs[staircase_num[row]].twinx()
                sub_axs.append(ax_twin)
            axs.append(sub_axs)
    # draw the figure for each time and take a snapshot
    stop_times = [simulations[sim].stop_time for sim in range(sim_num)]
    max_time = min([max_time]+stop_times)
    #R = [0 for sim in range(sim_num)]
    for sim in range(sim_num):
        sim_index = int(sim_delay / simulations[sim].sim_interval) * plot_index
        # plot staircase
        if staircase:
            for plot_pop in range(staircase_num[sim]):
                if simulations[sim].p_prop[0].ConsiderDensSwitch:
                    simulations[sim].staircase_plot(sim_index, plot_pop, axs[sim][plot_pop], True)
                    if staircase_num[sim] > 1:
                        if simulations[sim].p_prop[0].Move > simulations[sim].p_prop[1].Move:
                            axs[sim][plot_pop].set_title(titles[plot_pop], color=title_colors[simulations[sim].l_prop.PopulationNumber - 1 - plot_pop])
                        else:
                            axs[sim][plot_pop].set_title(titles[plot_pop], color=title_colors[plot_pop])
                else:
                    simulations[sim].staircase_plot(sim_index, simulations[sim].l_prop.PopulationNumber - 1 - plot_pop, axs[sim][plot_pop], True)
                    if staircase_num[sim] > 1:
                        axs[sim][plot_pop].set_title(titles[plot_pop], color=title_colors[plot_pop])
        # plot wildtype
        if wildtype:
            simulations[sim].wild_type_plot(sim_index, True, True, True, True, False, axs[sim][staircase_num[sim]], axs[sim][staircase_num[sim] + 1], False)
        # read gen_space_tot
        if simulations[sim].read_gen_space_tot(sim_index) is None:
            [Time, GenSpaceTot] = [np.inf, np.zeros((simulations[sim].l_prop.PopulationNumber, simulations[sim].l_prop.GridBound, simulations[sim].l_prop.GridBound))]
        else:
            [Time, GenSpaceTot] = simulations[sim].read_gen_space_tot(sim_index)
        # plot title
        R = simulations[sim].LD(GenSpaceTot)
        #if simulations[sim].next_founder(GenSpaceTot, R[sim]):
        #    R[sim] += 1
        subfigs[sim].suptitle(fig_titles[sim]+", Time: " + "{:6.1f}".format(Time) + "h, R: " + str(R))
        gen_space_tot.append(GenSpaceTot)
    # take a snapshot and clear the canvas
    return fig, gen_space_tot